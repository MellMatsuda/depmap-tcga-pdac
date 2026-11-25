# =============================================================================
# tcc_simulate.R — dados SIMULADOS com a mesma estrutura dos objetos reais.
# Serve SÓ para validar o código do TCCMELL.Rmd sem baixar DepMap/TCGA
# (params$simulate = TRUE). Nenhum número gerado aqui tem valor biológico.
# =============================================================================

simulate_tcc_data <- function(n_lines = 400, n_extra_expr = 30, n_tcga = 183,
                              n_genes = 1500, n_dep = 800, seed = 1) {
  set.seed(seed)
  lineages <- c("Pancreas", "Lung", "Breast", "Bowel", "Skin", "CNS/Brain",
                "Ovary/Fallopian Tube", "Esophagus/Stomach", "Kidney", "Lymphoid")
  n_all <- n_lines + n_extra_expr
  lin <- sample(lineages, n_all, TRUE, prob = c(.12, rep(.88 / 9, 9)))
  ids <- sprintf("ACH-%06d", seq_len(n_all))
  pt <- sprintf("PT-%05d", seq_len(n_all))
  share <- sample(n_all, 20); pt[share] <- pt[sample(n_all, 20)]
  subtype <- ifelse(lin == "Pancreas",
                    sample(c("Pancreatic Adenocarcinoma", "Pancreatic Neuroendocrine Tumor"),
                           n_all, TRUE, prob = c(.85, .15)), paste(lin, "tumor"))
  dm_model <- data.frame(
    ModelID = ids, PatientID = pt, CellLineName = paste0("LINE", seq_len(n_all)),
    OncotreeLineage = lin, OncotreePrimaryDisease = subtype, OncotreeSubtype = subtype,
    OncotreeCode = ifelse(subtype == "Pancreatic Adenocarcinoma", "PAAD", "OTHER"))

  # --- expressão: programas latentes + efeito de linhagem -------------------
  locus <- c("MTAP", "CDKN2A", "CDKN2B", "DMRTA1")
  named <- c(locus, "PRMT5", "MAT2A", "KRAS", "SMAD4", "TP53", "GATA6", "KRT5", "COL1A1",
             "ACTA2", "PTPRC", "PRSS1", "CPA1")
  genes <- c(named, sprintf("GENE%04d", seq_len(n_genes - length(named))))
  K <- 8
  L <- matrix(rnorm(length(lineages) * K, 0, 1.5), length(lineages), K, dimnames = list(lineages))
  z <- L[lin, ] + matrix(rnorm(n_all * K), n_all, K)
  Wg <- matrix(rnorm(n_genes * K, 0, 0.6), n_genes, K)
  mu_g <- runif(n_genes, 0.5, 7)
  Xl <- sweep(z %*% t(Wg), 2, mu_g, `+`) + matrix(rnorm(n_all * n_genes, 0, 0.5), n_all)
  dimnames(Xl) <- list(ids, genes)
  del_l <- runif(n_all) < ifelse(lin == "Pancreas", 0.3, 0.15)
  Xl[del_l, locus] <- rnorm(sum(del_l) * length(locus), 0.1, 0.1)
  Xl <- pmax(Xl, 0)

  # CN gene-level linear relativo (formato OmicsCNGeneWGS, 26Q3)
  cn <- matrix(rnorm(n_all * length(locus), 1, 0.15), n_all,
               dimnames = list(ids, locus))
  cn[del_l, ] <- abs(rnorm(sum(del_l) * length(locus), 0.05, 0.05))

  # --- dependências (Chronos-like) -------------------------------------------
  dep_genes <- c("PRMT5", "MAT2A", "KRAS", "WRN", "MYC", sprintf("DEP%04d", seq_len(n_dep - 5)))
  Bd <- matrix(rnorm(n_dep * K, 0, 0.12) * rbinom(n_dep * K, 1, 0.3), n_dep, K)
  base <- rnorm(n_dep, -0.15, 0.25)
  base[sample(n_dep, round(0.06 * n_dep))] <- -1.6       # pan-essenciais
  Y <- sweep(z %*% t(Bd), 2, base, `+`) + matrix(rnorm(n_all * n_dep, 0, 0.2), n_all)
  dimnames(Y) <- list(ids, dep_genes)
  Y[, "PRMT5"] <- -0.7 - 0.35 * del_l + rnorm(n_all, 0, 0.2)
  Y[, "MAT2A"] <- -0.6 - 0.30 * del_l + rnorm(n_all, 0, 0.2)
  Y[, "KRAS"]  <- -0.3 - 0.9 * (lin == "Pancreas") + rnorm(n_all, 0, 0.3)
  Y[sample(length(Y), 200)] <- NA
  keep_crispr <- seq_len(n_lines)                       # algumas linhas sem CRISPR

  # --- TCGA: tumores pancreáticos misturados com estroma ---------------------
  n_t <- n_tcga
  purity <- rbeta(n_t, 3, 4)
  zt <- matrix(L["Pancreas", ], n_t, K, byrow = TRUE) + matrix(rnorm(n_t * K), n_t, K)
  tum <- sweep(zt %*% t(Wg), 2, mu_g, `+`)
  del_t <- runif(n_t) < 0.25
  tum[del_t, genes %in% locus] <- 0.1
  stroma <- mu_g + rnorm(n_genes, 0, 1.2)
  stroma[genes %in% c("COL1A1", "ACTA2", "PTPRC", "PRSS1", "CPA1")] <- 10
  batch <- rnorm(n_genes, 0, 0.4)
  Xt <- purity * tum + (1 - purity) * matrix(stroma, n_t, n_genes, byrow = TRUE)
  Xt <- pmax(sweep(Xt, 2, batch, `+`) + matrix(rnorm(n_t * n_genes, 0, 0.4), n_t), 0)
  pids <- sprintf("TCGA-%s-%04d", sample(c("2J", "3A", "HZ", "IB", "US"), n_t, TRUE), seq_len(n_t))
  stype <- rep("Primary Tumor", n_t)
  stype[sample(n_t, 4)] <- "Solid Tissue Normal"
  vial <- ifelse(stype == "Primary Tumor", "01A", "11A")
  dupl <- sample(which(stype == "Primary Tumor"), 3)    # alíquotas repetidas
  bc <- paste0(pids, "-", vial, "-11R-A", sprintf("%03d", seq_len(n_t)), "-07")
  bc[dupl] <- paste0(pids[dupl - 1], "-01A-12R-A999-07")
  dimnames(Xt) <- list(bc, genes)
  gene_pos <- data.frame(symbol = genes, chr = sample(paste0("chr", 1:22), n_genes, TRUE),
                         start = sample.int(2e8, n_genes))
  gene_pos$chr[gene_pos$symbol %in% locus] <- "chr9"
  gene_pos$start[gene_pos$symbol %in% locus] <- c(21802000, 21967000, 22002000, 22447000)
  gene_pos$end <- gene_pos$start + 2e4
  tc <- list(expr = Xt, gene_pos = gene_pos, n_dup = 0L, n_genes_total = n_genes,
             gene_types = data.frame(symbol = genes, gene_type = "protein_coding"),
             coldata = data.frame(barcode = bc, sample_type = stype))

  # --- tabela tipo Raphael 2017 (150 amostras curadas) -----------------------
  s16 <- substr(bc, 1, 16)
  cur <- unique(s16[stype == "Primary Tumor"])
  cur <- sample(cur, min(150, length(cur)))
  m <- match(cur, s16)
  stage <- sample(c("stage ia.ib.iia_n0", "stage iib_n1", "stage iii", "stage iv"),
                  length(cur), TRUE, prob = c(.25, .67, .04, .04))
  hz <- 0.0015 * exp(0.6 * (stage != "stage ia.ib.iia_n0") + 0.4 * zt[m, 1])
  t_ev <- rexp(length(cur), hz); t_c <- runif(length(cur), 100, 2500)
  raph <- data.frame(
    `Tumor Sample ID` = cur,
    `ABSOLUTE Purity` = round(purity[m], 2),
    `Purity Class (high or low)` = ifelse(purity[m] > 0.33, "high", "low"),
    `mRNA Moffitt clusters (All 150 Samples) 1basal  2classical` =
      as.character(ifelse(zt[m, 2] > 0.5, 1, 2)),
    `KRAS Mutated (1 or 0)` = as.character(rbinom(length(cur), 1, 0.93)),
    `CDKN2A Status` = ifelse(del_t[m], "Deleted", NA),
    `Days to death` = ifelse(t_ev < t_c, round(t_ev), NA),
    `Follow up days` = round(pmin(t_ev, t_c)),
    `Censored 1=yes 0=no` = as.character(as.integer(t_ev >= t_c)),
    `AJCC pathologic tumor stage` = stage,
    `Histological type by RHH` = sample(c("Ductal adenocarcinoma", "Adenosquamous", "Other"),
                                        length(cur), TRUE, prob = c(.95, .02, .03)),
    `Age at initial pathologic diagnosis` = as.character(round(rnorm(length(cur), 65, 10))),
    Gender = sample(c("male", "female"), length(cur), TRUE),
    check.names = FALSE)
  tcga_cnv_mtap <- data.frame(sample = s16, min_copy_number = ifelse(del_t, 0, 2))

  list(dm_expr = Xl, dm_crispr = Y[keep_crispr, ], dm_model = dm_model,
       dm_cn = cn, tc = tc, raph = raph, clin = NULL, tcga_cnv_mtap = tcga_cnv_mtap)
}
