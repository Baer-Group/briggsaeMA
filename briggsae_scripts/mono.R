library(rstudioapi)
library(readxl)
library(data.table)
setwd(dirname(getActiveDocumentContext()$path))
cat("Working directory:", getwd(), "\n")

# Load mono repeat BED (0-based start, 1-based end - standard BED)
bed <- fread("briggsae_direct_perfect_mono_len5.bed",
             col.names = c("CHROM", "start0", "end1","repeat_unit"),
             header = FALSE)
bed[, start1 := start0 + 1]  # convert to 1-based for VCF matching
setkey(bed, CHROM, start1, end1)

classify_mono <- function(mutations) {
  dt <- as.data.table(mutations)
  dt[, pos_end := POS]          # point interval: POS to POS
  setkey(dt, CHROM, POS, pos_end)
  
  hits <- foverlaps(dt, bed,
                    by.x = c("CHROM", "POS", "pos_end"),
                    by.y = c("CHROM", "start1", "end1"),
                    type = "within", nomatch = NA)
  
  dt$mono_repeat <- !is.na(hits$start1)
  dt[, pos_end := NULL]
  dt
}

# Process both files, all sheets
sheets <- c("3x", "5x", "7x", "10x")
files  <- c("SR_final_mutation_list.xlsx", "SR_LR_final_mutation_list.xlsx")

for (f in files) {
  out_name <- sub(".xlsx", "_mono_classified.xlsx", f)
  wb_out   <- openxlsx::createWorkbook()
  
  for (s in sheets) {
    muts      <- read_excel(f, sheet = s)
    classified <- classify_mono(muts)
    
    openxlsx::addWorksheet(wb_out, s)
    openxlsx::writeData(wb_out, s, classified)
  }
  openxlsx::saveWorkbook(wb_out, out_name, overwrite = TRUE)
  cat("Saved:", out_name, "\n")
}
