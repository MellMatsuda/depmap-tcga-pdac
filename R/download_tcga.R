# Baixa e prepara TCGA-PAAD (STAR counts, clínica Biotab, CNV gene-level) para
# data/tcga/, com o mesmo código do chunk 'tcga-download' do TCCMELL.Rmd.
# Uso (a partir da raiz do projeto): Rscript R/download_tcga.R [cnv_workflow]
args <- commandArgs(trailingOnly = TRUE)
cnv_workflow <- if (length(args)) args[1] else "ASCAT3"
suppressPackageStartupMessages(library(TCGAbiolinks))
tcga_dir <- "data/tcga"
gdc_dir <- file.path(tcga_dir, "GDCdata")
dir.create(gdc_dir, recursive = TRUE, showWarnings = FALSE)
f_se <- file.path(tcga_dir, "TCGA_PAAD_STAR_SE.rds")
f_cl <- file.path(tcga_dir, "TCGA_PAAD_clinical_biotab.rds")
f_cn <- file.path(tcga_dir, "TCGA_PAAD_CNV_gene.rds")

if (!file.exists(f_se)) {
  q <- GDCquery(project = "TCGA-PAAD", data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
  GDCdownload(q, method = "api", files.per.chunk = 10, directory = gdc_dir)
  saveRDS(GDCprepare(q, directory = gdc_dir), f_se)
  message("OK: ", f_se)
}
if (!file.exists(f_cl)) {
  q <- GDCquery(project = "TCGA-PAAD", data.category = "Clinical",
                data.type = "Clinical Supplement", data.format = "BCR Biotab")
  GDCdownload(q, directory = gdc_dir)
  saveRDS(GDCprepare(q, directory = gdc_dir)$clinical_patient_paad, f_cl)
  message("OK: ", f_cl)
}
if (!file.exists(f_cn)) {
  cnv <- tryCatch({
    q <- GDCquery(project = "TCGA-PAAD", data.category = "Copy Number Variation",
                  data.type = "Gene Level Copy Number", workflow.type = cnv_workflow)
    GDCdownload(q, method = "api", files.per.chunk = 20, directory = gdc_dir)
    GDCprepare(q, directory = gdc_dir)
  }, error = function(e) { message("CNV falhou (", cnv_workflow, "): ", conditionMessage(e)); NULL })
  if (!is.null(cnv)) { saveRDS(cnv, f_cn); message("OK: ", f_cn) }
}
message("Concluído.")
