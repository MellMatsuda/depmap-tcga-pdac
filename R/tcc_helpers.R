# =============================================================================
# tcc_helpers.R — funções auxiliares do TCCMELL.Rmd (pipeline v2, sem scRNA)
# Cada função é pura (sem efeitos colaterais além de mensagens) para que o Rmd
# possa ajustá-las DENTRO de cada fold de treino, sem vazamento.
# =============================================================================

# ---- utilitários -------------------------------------------------------------

#' Avalia `expr` e guarda em RDS; na próxima execução lê do disco.
cached <- function(name, expr, use_cache = TRUE, dir = "results/cache") {
  f <- file.path(dir, paste0(name, ".rds"))
  if (use_cache && file.exists(f)) {
    message("[cache] lendo ", f)
    return(readRDS(f))
  }
  val <- force(expr)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(val, f)
  val
}

#' Primeira coluna existente entre os candidatos (nomes mudam entre releases).
pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) df[[hit[1]]] else rep(NA, nrow(df))
}

#' Correlação coluna a coluna entre duas matrizes de mesma dimensão.
#' Sem NA (o chamador garante). Coluna constante -> NA.
colwise_cor <- function(a, b, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  a <- as.matrix(a); b <- as.matrix(b)
  stopifnot(all(dim(a) == dim(b)))
  if (method == "spearman") {
    a <- apply(a, 2, rank); b <- apply(b, 2, rank)
    a <- matrix(a, ncol = ncol(b)); b <- matrix(b, ncol = ncol(a))
  }
  a <- sweep(a, 2, colMeans(a)); b <- sweep(b, 2, colMeans(b))
  den <- sqrt(colSums(a^2) * colSums(b^2))
  r <- colSums(a * b) / den
  r[!is.finite(r)] <- NA
  r
}

#' Como colwise_cor, mas tolera NA (casos completos por coluna). Devolve vetor nomeado.
colwise_cor_na <- function(a, b, method = "pearson") {
  a <- as.matrix(a); b <- as.matrix(b)
  stopifnot(all(dim(a) == dim(b)))
  r <- vapply(seq_len(ncol(a)), function(j) {
    ok <- !is.na(a[, j]) & !is.na(b[, j])
    if (sum(ok) < 3) return(NA_real_)
    suppressWarnings(stats::cor(a[ok, j], b[ok, j], method = method))
  }, numeric(1))
  stats::setNames(r, colnames(a))
}

# ---- 00. aquisição DepMap ----------------------------------------------------

#' Baixa arquivos de um release do DepMap pela listagem pública do portal.
#' TODO(verificar): o endpoint "api/download/files" e as colunas
#' (release, filename, url) foram usados por versões anteriores do portal;
#' se falhar, baixar manualmente em https://depmap.org/portal/data_page/ e
#' salvar os CSV em data/depmap/.
depmap_download <- function(release, filenames, dest_dir,
                            index_url = "https://depmap.org/portal/api/download/files") {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  idx <- tryCatch(data.table::fread(index_url, showProgress = FALSE),
                  error = function(e) NULL)
  if (is.null(idx) || !all(c("release", "filename", "url") %in% names(idx))) {
    warning("Não consegui ler a listagem do DepMap. Baixe manualmente: ",
            paste(filenames, collapse = ", "), " -> ", dest_dir)
    return(invisible(FALSE))
  }
  for (fn in filenames) {
    dest <- file.path(dest_dir, fn)
    if (file.exists(dest)) next
    row <- idx[idx$release == release & idx$filename == fn, ]
    if (!nrow(row)) { warning("Não achei ", fn, " em ", release); next }
    message("Baixando ", fn, " (", release, ")")
    utils::download.file(row$url[1], dest, mode = "wb")
  }
  invisible(TRUE)
}

#' Colunas de metadados conhecidas nos CSV do DepMap (formatos 23Q4–26Q3).
DEPMAP_META_COLS <- c("V1", "ProfileID", "SequencingID", "ModelID", "ModelConditionID",
                      "is_default_entry", "IsDefaultEntryForModel", "IsDefaultEntryForMC",
                      "is_default_entry_for_mc")

#' Lê uma matriz DepMap -> matriz modelos (ACH-) x símbolos HGNC.
#' Formatos aceitos:
#'  (a) antigo: 1ª coluna = ModelID, colunas de gene "SIMBOLO (ENTREZ)" (CRISPRGeneEffect);
#'  (b) 24Q4: SequencingID, ModelID, IsDefaultEntryForModel, ... + genes;
#'  (c) 26Q3: uma linha por profile — ProfileID, is_default_entry, ModelID + genes
#'      (OmicsExpressionTPMLogp1HumanProteinCodingGenes, OmicsCNGeneWGS).
#' Com default_only = TRUE fica só o profile default de cada modelo. Colunas de gene
#' = todas as que não são metadados; aceita "SIMBOLO (ENTREZ)" ou só "SIMBOLO".
read_depmap_matrix <- function(path, default_only = TRUE) {
  dt <- data.table::fread(path, showProgress = FALSE)
  nm <- names(dt)
  if (nm[1] == "") { data.table::setnames(dt, 1, "V1"); nm <- names(dt) }
  meta <- intersect(nm, DEPMAP_META_COLS)
  gene_cols <- setdiff(nm, meta)
  if (!length(gene_cols)) stop("Nenhuma coluna de gene em ", path)
  n_rows_raw <- nrow(dt)
  def_col <- intersect(c("is_default_entry", "IsDefaultEntryForModel"), nm)
  if (default_only && length(def_col)) {
    v <- toupper(as.character(dt[[def_col[1]]]))
    dt <- dt[v %in% c("TRUE", "YES", "1", "T"), ]
  }
  id <- if ("ModelID" %in% nm) dt[["ModelID"]] else dt[[1]]
  if (!all(grepl("^ACH-", id[!is.na(id)])))
    stop("IDs de linha não são ModelID (ACH-) em ", path,
         ". Se houver só ProfileID, mapeie via OmicsProfiles.csv.")
  m <- as.matrix(dt[, gene_cols, with = FALSE])
  storage.mode(m) <- "double"
  has_entrez <- grepl(" \\(\\d+\\)$", gene_cols)
  sym <- sub(" \\(\\d+\\)$", "", gene_cols)
  entrez <- ifelse(has_entrez, sub("^.* \\((\\d+)\\)$", "\\1", gene_cols), NA_character_)
  rownames(m) <- id
  colnames(m) <- sym
  # HGNC 1:1 — símbolos duplicados (raros) são descartados, não somados
  dup <- sym %in% sym[duplicated(sym)]
  keep_row <- !is.na(id) & !duplicated(id)
  m <- m[keep_row, !dup, drop = FALSE]
  attr(m, "entrez") <- stats::setNames(entrez[!dup], sym[!dup])
  attr(m, "n_dup_symbols") <- sum(dup)
  attr(m, "n_rows_raw") <- n_rows_raw
  attr(m, "n_dup_models_dropped") <- sum(duplicated(id[!is.na(id)]))
  attr(m, "default_filter") <- if (default_only && length(def_col)) def_col[1] else "nenhum"
  m
}

# ---- 00/01. TCGA -------------------------------------------------------------

#' SummarizedExperiment do GDC (STAR - Counts) -> matriz amostras x símbolos
#' em log2(TPM+1). Adaptado de tcga_march1.Rmd (l. 149-160): protein-coding,
#' assay tpm_unstrand, duplicatas de símbolo somadas (decisão declarada).
prep_tcga_expr <- function(se) {
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  gene_types <- unique(data.frame(symbol = rd$gene_name, gene_type = rd$gene_type))
  n_total <- nrow(rd)
  keep <- rd$gene_type == "protein_coding" & !is.na(rd$gene_name) & rd$gene_name != ""
  se <- se[keep, ]
  rd <- rd[keep, ]
  tpm <- SummarizedExperiment::assay(se, "tpm_unstrand")
  sym <- rd$gene_name
  n_dup <- sum(duplicated(sym))
  tpm <- rowsum(tpm, group = sym)
  rr <- as.data.frame(SummarizedExperiment::rowRanges(se))
  pos <- data.frame(symbol = sym, chr = as.character(rr$seqnames),
                    start = rr$start, end = rr$end)
  pos <- pos[!duplicated(pos$symbol), ]
  list(expr = t(log2(tpm + 1)), gene_pos = pos, n_dup = n_dup, gene_types = gene_types,
       n_genes_total = n_total,
       coldata = as.data.frame(SummarizedExperiment::colData(se)))
}

#' Genes a até `flank` pb de um gene âncora (ex.: locus 9p21 em torno de MTAP).
genes_near <- function(gene_pos, anchor = "MTAP", flank = 1e6) {
  a <- gene_pos[gene_pos$symbol == anchor, ]
  if (!nrow(a)) return(character())
  hit <- gene_pos$chr == a$chr[1] &
    gene_pos$end >= a$start[1] - flank & gene_pos$start <= a$end[1] + flank
  gene_pos$symbol[hit]
}

# ---- 05. alinhamento tipo Celligner ------------------------------------------

#' PCA por SVD; devolve centro, rotação (p x k) e scores (n x k).
pca_fit <- function(X, k) {
  k <- min(k, nrow(X) - 1, ncol(X))
  center <- colMeans(X)
  Xc <- sweep(X, 2, center)
  s <- svd(Xc, nu = 0, nv = k)
  list(center = center, rotation = s$v, scores = Xc %*% s$v, sdev = s$d / sqrt(nrow(X) - 1))
}

#' Direções de PCA contrastivo: autovetores de Cov(tumor) - Cov(linhagem).
#' Como a matriz tem posto <= n_t + n_l, resolve no espaço das amostras:
#' M = [A; B] = U D V',  Cov_dif = V (D U' S U D) V',  S = diag(+1 tumor, -1 linha).
cpca_directions <- function(Tc, Lc, k = 4) {
  A <- Tc / sqrt(nrow(Tc) - 1)
  B <- Lc / sqrt(nrow(Lc) - 1)
  s <- svd(rbind(A, B))
  keep <- s$d > 1e-8 * s$d[1]
  U <- s$u[, keep, drop = FALSE]; D <- s$d[keep]; V <- s$v[, keep, drop = FALSE]
  sgn <- c(rep(1, nrow(A)), rep(-1, nrow(B)))
  C <- (D * t(U)) %*% (sgn * U)
  C <- sweep(C, 2, D, `*`)
  e <- eigen((C + t(C)) / 2, symmetric = TRUE)
  W <- V %*% e$vectors[, seq_len(k), drop = FALSE]
  list(W = W, values = e$values[seq_len(k)])
}

#' Correção MNN (Haghverdi 2018, como usada no Celligner): leva os TUMORES para
#' perto das linhagens. Pares mútuos são achados no espaço de PCs; os vetores de
#' correção (espaço de genes) são suavizados por kernel gaussiano entre tumores.
mnn_correction <- function(Lr, Tr, Sl, St, k_t2l = 50, k_l2t = 5, sigma_frac = 0.5) {
  nn_t2l <- FNN::get.knnx(Sl, St, k = min(k_t2l, nrow(Sl)))$nn.index
  nn_l2t <- FNN::get.knnx(St, Sl, k = min(k_l2t, nrow(St)))$nn.index
  pairs <- data.frame(i = rep(seq_len(nrow(St)), ncol(nn_t2l)), j = as.vector(nn_t2l))
  mutual <- mapply(function(i, j) i %in% nn_l2t[j, ], pairs$i, pairs$j)
  pairs <- pairs[mutual, , drop = FALSE]
  corr <- matrix(0, nrow(Tr), ncol(Tr), dimnames = dimnames(Tr))
  if (!nrow(pairs)) {
    warning("MNN: nenhum par mútuo encontrado; tumores não corrigidos.")
    return(list(correction = corr, n_pairs = 0, n_tumors_paired = 0))
  }
  ti <- sort(unique(pairs$i))
  vec <- t(vapply(ti, function(i) {
    js <- pairs$j[pairs$i == i]
    colMeans(Lr[js, , drop = FALSE]) - Tr[i, ]
  }, numeric(ncol(Tr))))
  d_all <- as.matrix(stats::dist(St))
  sigma <- sigma_frac * stats::median(d_all[upper.tri(d_all)])
  Wk <- exp(-d_all[, ti, drop = FALSE]^2 / (2 * sigma^2))
  rs <- rowSums(Wk)
  zero <- rs < 1e-12
  if (any(zero)) {  # sem vizinho com par: usa o tumor pareado mais próximo
    nearest <- apply(d_all[zero, ti, drop = FALSE], 1, which.min)
    Wk[zero, ] <- 0
    Wk[cbind(which(zero), nearest)] <- 1
    rs[zero] <- 1
  }
  corr[] <- (Wk %*% vec) / rs
  list(correction = corr, n_pairs = nrow(pairs), n_tumors_paired = length(ti))
}

#' Ajusta o alinhamento usando SÓ as linhagens de treino + TCGA (sem Y).
#' Passos: genes de maior variância (X apenas) -> centralização por fonte ->
#' remoção dos cPCs -> MNN nos tumores -> PCA combinada (features).
#' `contrast_idx`: índices das linhagens usadas no contraste da cPCA (NULL =
#' todas; variante de diagnóstico = só PDAC). `return_matrices`: devolve as
#' matrizes intermediárias e os scores de cada estágio (para exploração).
fit_alignment <- function(Xl, Xt, n_genes = 5000, n_cpc = 4, n_pc = 50,
                          n_pc_mnn = 50, k_t2l = 50, k_l2t = 5, sigma_frac = 0.5,
                          exclude_genes = character(), contrast_idx = NULL,
                          anchor_idx = NULL, return_matrices = FALSE,
                          offset_target = NULL) {
  # `offset_target`: vetor (nomeado por gene) com o centróide-alvo dos tumores,
  # alternativa a `anchor_idx` quando as linhagens-âncora não estão em Xl
  # (ex.: PDAC do teste no leave-PDAC-out, transdutivo, só X).
  # `anchor_idx`: linhagens cujo centróide recebe o centróide dos tumores.
  # Motivo: com TCGA só PAAD, centralizar os tumores na própria média apaga o
  # sinal pancreático médio (o Celligner original centraliza um TCGA
  # pan-câncer). Com âncora = linhagens PDAC de treino, o centróide tumoral é
  # levado ao centróide PDAC das linhagens. NULL = Celligner puro.
  genes <- setdiff(intersect(colnames(Xl), colnames(Xt)), exclude_genes)
  v <- matrixStats::colVars(rbind(Xl[, genes, drop = FALSE], Xt[, genes, drop = FALSE]))
  genes <- genes[order(v, decreasing = TRUE)][seq_len(min(n_genes, length(genes)))]
  Xl <- Xl[, genes, drop = FALSE]; Xt <- Xt[, genes, drop = FALSE]
  mu_l <- colMeans(Xl); mu_t <- colMeans(Xt)
  Lc <- sweep(Xl, 2, mu_l); Tc <- sweep(Xt, 2, mu_t)
  Lcon <- if (is.null(contrast_idx)) Lc else
    sweep(Xl[contrast_idx, , drop = FALSE], 2, colMeans(Xl[contrast_idx, , drop = FALSE]))
  cp <- if (n_cpc > 0) cpca_directions(Tc, Lcon, n_cpc) else list(W = matrix(0, ncol(Lc), 0), values = numeric())
  W <- cp$W
  offset <- if (!is.null(offset_target)) offset_target[genes] - mu_l else
    if (is.null(anchor_idx)) rep(0, ncol(Lc)) else
    colMeans(Xl[anchor_idx, , drop = FALSE]) - mu_l
  Tc_cpc <- Tc                                  # scores de cPC sem o offset
  Tc <- sweep(Tc, 2, offset, `+`)
  Lr <- Lc - (Lc %*% W) %*% t(W)
  Tr <- Tc - (Tc %*% W) %*% t(W)
  nl <- nrow(Lr)
  pc0 <- pca_fit(rbind(Lr, Tr), n_pc_mnn)
  mnn <- mnn_correction(Lr, Tr, pc0$scores[seq_len(nl), , drop = FALSE],
                        pc0$scores[-seq_len(nl), , drop = FALSE],
                        k_t2l, k_l2t, sigma_frac)
  Ta <- Tr + mnn$correction
  pc <- pca_fit(rbind(Lr, Ta), n_pc)
  colnames(pc$scores) <- paste0("PC", seq_len(ncol(pc$scores)))
  out <- list(genes = genes, mu_l = mu_l, mu_t = mu_t, W = W, cpc_values = cp$values,
              pc_center = pc$center, rotation = pc$rotation, pc_sdev = pc$sdev,
              feat_lines = pc$scores[seq_len(nl), , drop = FALSE],
              feat_tcga = pc$scores[-seq_len(nl), , drop = FALSE],
              n_pairs = mnn$n_pairs, n_tumors_paired = mnn$n_tumors_paired)
  if (return_matrices) {
    out$cpc_scores_tcga <- Tc_cpc %*% W      # projeção dos tumores nos cPCs
    out$cpc_scores_lines <- Lc %*% W
    out$scores_cpca_only <- pc0$scores       # estágio: centralizado + sem cPCs
    out$dist_nn_before <- FNN::get.knnx(pc0$scores[seq_len(nl), , drop = FALSE],
                                        pc0$scores[-seq_len(nl), , drop = FALSE], k = 1)$nn.dist[, 1]
    out$dist_nn_after <- FNN::get.knnx(out$feat_lines, out$feat_tcga, k = 1)$nn.dist[, 1]
  }
  out
}

#' Projeta linhagens novas (held-out) com os parâmetros do treino.
#' Linhagens não recebem correção MNN (no Celligner só os tumores são movidos).
transform_lines <- function(fit, Xnew) {
  Xc <- sweep(Xnew[, fit$genes, drop = FALSE], 2, fit$mu_l)
  Xr <- Xc - (Xc %*% fit$W) %*% t(fit$W)
  out <- sweep(Xr, 2, fit$pc_center) %*% fit$rotation
  colnames(out) <- paste0("PC", seq_len(ncol(out)))
  out
}

#' Checagem do alinhamento: para cada tumor, fração dos k vizinhos (entre
#' linhagens) que são da linhagem-alvo, e fração de tumores entre os k vizinhos
#' no espaço combinado (mistura de fonte; menor = mais misturado).
alignment_metrics <- function(F_lines, F_tcga, is_pancreas, is_pdac, k = 10) {
  nn <- FNN::get.knnx(F_lines, F_tcga, k = k)$nn.index
  frac_panc <- rowMeans(matrix(is_pancreas[nn], nrow = nrow(nn)))
  frac_pdac <- rowMeans(matrix(is_pdac[nn], nrow = nrow(nn)))
  all_f <- rbind(F_lines, F_tcga)
  src_t <- c(rep(FALSE, nrow(F_lines)), rep(TRUE, nrow(F_tcga)))
  nn2 <- FNN::get.knnx(all_f, F_tcga, k = k + 1)$nn.index[, -1, drop = FALSE]
  frac_tumor <- rowMeans(matrix(src_t[nn2], nrow = nrow(nn2)))
  c(viz_pancreas = mean(frac_panc), acaso_pancreas = mean(is_pancreas),
    viz_pdac = mean(frac_pdac), acaso_pdac = mean(is_pdac),
    viz_mesma_fonte_tumor = mean(frac_tumor),
    acaso_fonte = nrow(F_tcga) / nrow(all_f))
}

#' Escore simples de marcadores: média dos z-scores (dentro da coorte) dos
#' genes presentes. Não é ESTIMATE; é um proxy declarado.
marker_score <- function(X, genes) {
  g <- intersect(genes, colnames(X))
  if (!length(g)) return(rep(NA_real_, nrow(X)))
  z <- scale(X[, g, drop = FALSE])
  z[!is.finite(z)] <- NA
  rowMeans(z, na.rm = TRUE)
}

# ---- 06. alvo Y (dentro do fold) ---------------------------------------------

#' Lê CRISPRInferredCommonEssentials.csv (uma coluna "Essentials", "SIMBOLO (ENTREZ)").
read_common_essentials <- function(path) {
  v <- data.table::fread(path, header = TRUE)[[1]]
  unique(sub(" \\(\\d+\\)$", "", as.character(v)))
}

#' Seleciona genes-alvo usando SÓ as linhagens de treino.
#' `pan_genes`: lista externa de pan-essenciais (ex.: CRISPRInferredCommonEssentials);
#' se NULL, usa a regra interna (fração de linhagens de treino com efeito < pan_thr).
#' Genes com mais de `max_na` de NA no treino ficam fora (o mgaussian exige Y completo).
select_targets <- function(Ytr, n_top = 1000, pan_thr = -0.5, pan_frac = 0.9,
                           anchors = character(), pan_genes = NULL, max_na = 0.05) {
  pan <- if (is.null(pan_genes)) {
    colnames(Ytr)[colMeans(Ytr < pan_thr, na.rm = TRUE) >= pan_frac]
  } else intersect(pan_genes, colnames(Ytr))
  na_ok <- colMeans(is.na(Ytr)) <= max_na
  sds <- matrixStats::colSds(Ytr, na.rm = TRUE)
  names(sds) <- colnames(Ytr)
  cand <- setdiff(colnames(Ytr)[na_ok & is.finite(sds)], pan)
  top <- names(sort(sds[cand], decreasing = TRUE))[seq_len(min(n_top, length(cand)))]
  union(top, intersect(anchors, colnames(Ytr)[na_ok]))
}

#' Imputa NA pela média da coluna (usar só no Y de TREINO; a avaliação ignora NA).
impute_col_mean <- function(Y) {
  mu <- colMeans(Y, na.rm = TRUE)
  idx <- which(is.na(Y), arr.ind = TRUE)
  if (nrow(idx)) Y[idx] <- mu[idx[, 2]]
  Y
}

# ---- 07. features dentro do fold ---------------------------------------------

#' Refaz o alinhamento com as linhagens de treino + TCGA (sem Y) e projeta as
#' linhagens de teste. `anchor` = "pdac" usa o centróide das linhagens PDAC do
#' TREINO; "pdac_test" usa as PDAC do teste (transdutivo, só X; sensibilidade
#' do leave-PDAC-out); "none" = Celligner puro.
fold_features <- function(Xl, Xt, tr, te, is_pdac, cfg, exclude_genes = character(),
                          anchor = c("pdac", "none", "pdac_test")) {
  anchor <- match.arg(anchor)
  anchor_idx <- NULL; offset_target <- NULL
  if (anchor == "pdac") {
    anchor_idx <- which(is_pdac[tr])
    if (!length(anchor_idx)) anchor_idx <- NULL
  } else if (anchor == "pdac_test") {
    # centróide das PDAC do TESTE (só X, sem Y): sensibilidade transdutiva declarada
    offset_target <- colMeans(Xl[te[is_pdac[te]], , drop = FALSE])
  }
  fit <- fit_alignment(Xl[tr, , drop = FALSE], Xt, n_genes = cfg$n_genes_align,
                       n_cpc = cfg$n_cpc, n_pc = cfg$n_pc, k_t2l = cfg$k_t2l,
                       k_l2t = cfg$k_l2t, exclude_genes = exclude_genes,
                       anchor_idx = anchor_idx, offset_target = offset_target)
  Fte <- transform_lines(fit, Xl[te, , drop = FALSE])
  rownames(fit$feat_lines) <- rownames(Xl)[tr]
  rownames(Fte) <- rownames(Xl)[te]
  rownames(fit$feat_tcga) <- rownames(Xt)
  list(Ftr = fit$feat_lines, Fte = Fte, Ftcga = fit$feat_tcga, tr = tr, te = te,
       genes = fit$genes, n_pairs = fit$n_pairs)
}

# ---- 09. folds agrupados -----------------------------------------------------

#' K folds agrupados (todas as linhagens de um mesmo paciente no mesmo fold),
#' estratificados por um rótulo lógico (ex.: PDAC) para espalhá-lo entre folds.
make_group_folds <- function(groups, K = 5, strata = NULL, seed = 1) {
  set.seed(seed)
  groups <- as.character(groups)
  ug <- unique(groups)
  gs <- if (is.null(strata)) rep("all", length(ug)) else
    as.character(tapply(strata, groups, any)[ug])
  fold_g <- stats::setNames(integer(length(ug)), ug)
  for (s in unique(gs)) {
    g <- ug[gs == s]
    g <- g[sample.int(length(g))]
    fold_g[g] <- rep_len(sample.int(K), length(g))
  }
  unname(fold_g[groups])
}

# ---- 10. modelos -------------------------------------------------------------

#' k-NN multi-saída: média de Y dos k vizinhos de treino no espaço de features.
knn_predict <- function(Ftr, Ytr, Fte, k = 20) {
  idx <- FNN::get.knnx(Ftr, Fte, k = min(k, nrow(Ftr)))$nn.index
  out <- t(apply(idx, 1, function(ii) colMeans(Ytr[ii, , drop = FALSE], na.rm = TRUE)))
  out <- matrix(out, nrow = nrow(Fte), dimnames = list(rownames(Fte), colnames(Ytr)))
  out
}

#' Elastic net por gene (glmnet), lambda escolhido por CV interna agrupada.
#' glmnet padroniza X com média/DP do treino internamente (sem vazamento).
en_fit_predict <- function(Ftr, Ytr, Fte, foldid, alpha = 0.5, keep_fits = FALSE) {
  pred <- matrix(NA_real_, nrow(Fte), ncol(Ytr), dimnames = list(rownames(Fte), colnames(Ytr)))
  fits <- if (keep_fits) vector("list", ncol(Ytr)) else NULL
  nnz <- integer(ncol(Ytr))
  for (g in seq_len(ncol(Ytr))) {
    y <- Ytr[, g]
    ok <- !is.na(y)
    fid <- as.integer(factor(foldid[ok]))
    if (sum(ok) < 20 || length(unique(fid)) < 3) next
    cv <- glmnet::cv.glmnet(Ftr[ok, , drop = FALSE], y[ok], alpha = alpha,
                            foldid = fid, standardize = TRUE)
    pred[, g] <- as.numeric(stats::predict(cv, newx = Fte, s = "lambda.min"))
    nnz[g] <- cv$nzero[cv$index["min", 1]]
    if (keep_fits) fits[[g]] <- cv
  }
  if (keep_fits) names(fits) <- colnames(Ytr)
  list(pred = pred, nzero = stats::setNames(nnz, colnames(Ytr)), fits = fits)
}

#' Elastic net multi-resposta (glmnet family = "mgaussian"): um único ajuste
#' para todos os genes de Y. A penalidade é de grupo por feature (cada feature
#' entra ou sai para TODOS os genes ao mesmo tempo); declarar no texto.
#' Lambda por CV interna agrupada (foldid). X padronizado pelo glmnet com
#' média/DP do treino; Y centrado internamente (intercepto por gene).
mgaussian_fit_predict <- function(Ftr, Ytr, Fte, foldid, alpha = 0.5) {
  cv <- glmnet::cv.glmnet(Ftr, Ytr, family = "mgaussian", alpha = alpha,
                          foldid = foldid, standardize = TRUE)
  pred <- stats::predict(cv, newx = Fte, s = "lambda.min")
  pred <- matrix(pred[, , 1], nrow = nrow(Fte), dimnames = list(rownames(Fte), colnames(Ytr)))
  list(pred = pred, lambda = cv$lambda.min,
       nzero = unname(cv$nzero[which(cv$lambda == cv$lambda.min)]))
}

#' Roda os modelos de um fold a partir das features `ff` (saída de fold_features).
#' Modelos: trivial (média de treino), knn, en_A (mgaussian só com PCs; sem
#' 9p21/CNV) e en_B (PCs + `extra`, ex.: status MTAP/CDKN2A).
#' `perm_seed` != NULL embaralha as linhas de Y no treino (nulo por permutação);
#' X e alinhamento não mudam, porque não dependem de Y.
run_fold_models <- function(ff, Y, cfg, groups, pan_genes = NULL,
                            models = c("trivial", "knn", "en_A", "en_B"),
                            perm_seed = NULL, extra = NULL, targets = NULL) {
  Ytr <- Y[ff$tr, , drop = FALSE]
  if (is.null(targets))
    targets <- select_targets(Ytr, cfg$n_genes_y, anchors = cfg$anchors, pan_genes = pan_genes)
  Ytr <- impute_col_mean(Ytr[, targets, drop = FALSE])
  if (!is.null(perm_seed)) {
    set.seed(perm_seed)
    Ytr <- Ytr[sample.int(nrow(Ytr)), , drop = FALSE]
  }
  out <- list(targets = targets, te = ff$te, info = list())
  if ("trivial" %in% models) {
    mu <- colMeans(Ytr)
    out$trivial <- matrix(mu, length(ff$te), length(mu), byrow = TRUE,
                          dimnames = list(rownames(ff$Fte), targets))
  }
  if ("knn" %in% models) out$knn <- knn_predict(ff$Ftr, Ytr, ff$Fte, cfg$knn_k)
  inner <- make_group_folds(groups[ff$tr], K = cfg$inner_k, seed = cfg$seed)
  if ("en_A" %in% models) {
    m <- mgaussian_fit_predict(ff$Ftr, Ytr, ff$Fte, inner, cfg$en_alpha)
    out$en_A <- m$pred; out$info$en_A <- m[c("lambda", "nzero")]
  }
  if ("en_B" %in% models && !is.null(extra)) {
    m <- mgaussian_fit_predict(cbind(ff$Ftr, extra[ff$tr, , drop = FALSE]), Ytr,
                               cbind(ff$Fte, extra[ff$te, , drop = FALSE]), inner, cfg$en_alpha)
    out$en_B <- m$pred; out$info$en_B <- m[c("lambda", "nzero")]
  }
  out
}

#' Junta as predições held-out de vários folds numa matriz linhas x genes.
assemble_preds <- function(fold_res, model, row_ids) {
  genes <- unique(unlist(lapply(fold_res, function(f) colnames(f[[model]]))))
  P <- matrix(NA_real_, length(row_ids), length(genes), dimnames = list(row_ids, genes))
  for (f in fold_res) {
    if (is.null(f[[model]])) next
    P[rownames(f[[model]]), colnames(f[[model]])] <- f[[model]]
  }
  P
}

# ---- 12. avaliação -----------------------------------------------------------

#' r por gene (Pearson/Spearman) nas linhagens `idx`, com IC bootstrap
#' (reamostragem de linhagens) e RMSE. Só genes completos em `idx`.
eval_genes <- function(P, Y, idx, genes = NULL, n_boot = 1000, seed = 1) {
  if (is.null(genes)) genes <- colnames(P)
  genes <- intersect(genes, colnames(Y))
  ids <- rownames(P)[idx]
  Pm <- P[ids, genes, drop = FALSE]; Ym <- Y[ids, genes, drop = FALSE]
  ok <- colSums(is.na(Pm)) == 0 & colSums(is.na(Ym)) == 0
  Pm <- Pm[, ok, drop = FALSE]; Ym <- Ym[, ok, drop = FALSE]
  r_p <- colwise_cor(Pm, Ym, "pearson")
  r_s <- colwise_cor(Pm, Ym, "spearman")
  set.seed(seed)
  boot <- matrix(NA_real_, n_boot, ncol(Pm))
  for (b in seq_len(n_boot)) {
    ii <- sample.int(nrow(Pm), replace = TRUE)
    boot[b, ] <- suppressWarnings(colwise_cor(Pm[ii, , drop = FALSE], Ym[ii, , drop = FALSE]))
  }
  data.frame(gene = colnames(Pm), n = nrow(Pm), r_pearson = r_p, r_spearman = r_s,
             ci_low = apply(boot, 2, stats::quantile, 0.025, na.rm = TRUE),
             ci_high = apply(boot, 2, stats::quantile, 0.975, na.rm = TRUE),
             rmse = sqrt(colMeans((Pm - Ym)^2)),
             mse = colMeans((Pm - Ym)^2),
             row.names = NULL)
}

#' p empírico unilateral (r observado maior que o nulo agrupado entre genes).
perm_pvalue <- function(r_obs, r_null) {
  r_null <- r_null[is.finite(r_null)]
  vapply(r_obs, function(r) (1 + sum(r_null >= r)) / (1 + length(r_null)), numeric(1))
}

# ---- 14. validação nos pacientes ---------------------------------------------

#' Efeito de MTAP-del na dependência predita, ajustado por pureza, com IC (lm)
#' e p por permutação do rótulo de MTAP.
mtap_test <- function(yhat, mtap_del, purity, n_perm = 10000, seed = 1) {
  d <- stats::na.omit(data.frame(y = yhat, del = as.integer(mtap_del), pur = purity))
  fit <- stats::lm(y ~ del + pur, data = d)
  est <- unname(stats::coef(fit)["del"])
  ci <- stats::confint(fit)["del", ]
  X <- cbind(1, d$del, d$pur)
  set.seed(seed)
  null <- replicate(n_perm, {
    X[, 2] <- sample(d$del)
    stats::lm.fit(X, d$y)$coefficients[2]
  })
  data.frame(n = nrow(d), n_del = sum(d$del), efeito = est,
             ic_low = unname(ci[1]), ic_high = unname(ci[2]),
             dif_medias = mean(d$y[d$del == 1]) - mean(d$y[d$del == 0]),
             p_perm = (1 + sum(abs(null) >= abs(est))) / (1 + n_perm))
}

#' Efeito padronizado de MTAP-del (ajustado por pureza) para TODOS os genes
#' preditos: distribuição de referência para o par "falso" gene-marcador.
mtap_effect_all_genes <- function(Yhat, mtap_del, purity) {
  ok <- !is.na(mtap_del) & !is.na(purity)
  Ys <- scale(Yhat[ok, , drop = FALSE])
  X <- cbind(1, as.integer(mtap_del[ok]), purity[ok])
  B <- qr.coef(qr(X), Ys)
  stats::setNames(B[2, ], colnames(Yhat))
}
