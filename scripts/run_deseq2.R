suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
  stop(
    "Usage: Rscript run_deseq2.R ",
    "<count.csv> <coldata.csv> <output_dir>"
  )
}

counts_file <- args[[1]]
coldata_file <- args[[2]]
output_dir <- args[[3]]

control_level <- "complete"
alpha <- 0.05
ape_method <- "nbinomCR"


read_counts <- function(path) {
  if (!file.exists(path)) {
    stop("Counts file not found: ", path)
  }

  x <- read.csv(
    path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  if (!nrow(x) || !ncol(x)) {
    stop("Counts file is empty: ", path)
  }

  rownames(x) <- trimws(rownames(x))
  colnames(x) <- trimws(colnames(x))

  if (any(!nzchar(rownames(x)))) {
    stop("Counts contain empty gene IDs")
  }

  if (any(!nzchar(colnames(x)))) {
    stop("Counts contain empty sample IDs")
  }

  if (anyDuplicated(rownames(x))) {
    stop("Counts contain duplicated gene IDs")
  }

  if (anyDuplicated(colnames(x))) {
    stop("Counts contain duplicated sample IDs")
  }

  x <- as.matrix(x)
  suppressWarnings(storage.mode(x) <- "numeric")

  if (anyNA(x)) {
    stop("Counts contain missing or non-numeric values")
  }

  if (any(x < 0)) {
    stop("Counts contain negative values")
  }

  if (any(x != round(x))) {
    stop("Counts must be unnormalized integer values")
  }

  storage.mode(x) <- "integer"

  x
}


read_coldata <- function(path) {
  if (!file.exists(path)) {
    stop("Sample table not found: ", path)
  }

  x <- read.csv(
    path,
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    comment.char = "",
    colClasses = "character",
    stringsAsFactors = FALSE,
    blank.lines.skip = TRUE
  )

  if (nrow(x)) {
    nonempty <- apply(
      x,
      1L,
      function(values) {
        any(!is.na(values) & nzchar(trimws(values)))
      }
    )

    x <- x[nonempty, , drop = FALSE]
  }

  if (!nrow(x)) {
    stop("Sample table contains no data rows: ", path)
  }

  required <- c(
    "sample_id",
    "condition",
    "patient"
  )

  missing <- setdiff(
    required,
    colnames(x)
  )

  if (length(missing)) {
    stop(
      "Sample table is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  x[required] <- lapply(
    x[required],
    trimws
  )

  if (
    anyNA(x[required]) ||
    any(x[required] == "")
  ) {
    stop(
      "sample_id, condition and patient must not be empty"
    )
  }

  if (anyDuplicated(x$sample_id)) {
    stop(
      "Sample table contains duplicated sample_id values"
    )
  }

  rownames(x) <- x$sample_id

  x
}


write_matrix <- function(x, path) {
  write.table(
    as.data.frame(x),
    file = path,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}


format_results <- function(res) {
  x <- as.data.frame(res)

  data.frame(
    ensembl_id = rownames(x),
    x,
    row.names = NULL,
    check.names = FALSE
  )
}


counts <- read_counts(counts_file)
coldata <- read_coldata(coldata_file)

missing_in_coldata <- setdiff(
  colnames(counts),
  rownames(coldata)
)

missing_in_counts <- setdiff(
  rownames(coldata),
  colnames(counts)
)

if (
  length(missing_in_coldata) ||
  length(missing_in_counts)
) {
  stop(
    "Sample IDs do not match.\n",
    "Missing in coldata: ",
    paste(
      missing_in_coldata,
      collapse = ", "
    ),
    "\n",
    "Missing in counts: ",
    paste(
      missing_in_counts,
      collapse = ", "
    )
  )
}

coldata <- coldata[
  colnames(counts),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(counts),
    rownames(coldata)
  )
)

if (!control_level %in% coldata$condition) {
  stop(
    "Control level not found: ",
    control_level
  )
}

coldata$patient <- factor(
  coldata$patient
)

condition_values <- coldata$condition

condition_levels <- levels(
  relevel(
    factor(condition_values),
    ref = control_level
  )
)

if (nlevels(coldata$patient) != 3L) {
  stop(
    "Expected exactly three patients"
  )
}

if (length(condition_levels) != 5L) {
  stop(
    "Expected exactly five conditions including the control"
  )
}

condition_codes <- gsub(
  "[^A-Za-z0-9_.]+",
  "_",
  condition_levels
)

condition_codes <- gsub(
  "^_+|_+$",
  "",
  condition_codes
)

empty_codes <- !nzchar(
  condition_codes
)

condition_codes[empty_codes] <- paste0(
  "condition_",
  which(empty_codes)
)

condition_codes <- make.unique(
  condition_codes,
  sep = "_"
)

names(condition_codes) <- condition_levels

coldata$condition <- factor(
  unname(
    condition_codes[
      condition_values
    ]
  ),
  levels = unname(
    condition_codes
  )
)

sample_layout <- table(
  patient = coldata$patient,
  condition = condition_values
)

if (any(sample_layout != 1L)) {
  stop(
    "Expected exactly one sample for every patient-condition combination"
  )
}

design_matrix <- model.matrix(
  ~ patient + condition,
  data = coldata
)

if (
  qr(design_matrix)$rank <
  ncol(design_matrix)
) {
  stop(
    "The design matrix is not full rank"
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "Counts:   ",
  normalizePath(counts_file),
  "\n",
  sep = ""
)

cat(
  "Coldata:  ",
  normalizePath(coldata_file),
  "\n",
  sep = ""
)

cat(
  "Output:   ",
  normalizePath(output_dir),
  "\n\n",
  sep = ""
)

print(sample_layout)


dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ patient + condition
)

dds <- dds[
  rowSums(
    counts(dds)
  ) > 0,
]

dds <- DESeq(dds)


norm_counts <- counts(
  dds,
  normalized = TRUE
)

vsd <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(vsd)


write_matrix(
  norm_counts,
  file.path(
    output_dir,
    "counts_normalized_deseq2.tsv"
  )
)

write_matrix(
  vst_matrix,
  file.path(
    output_dir,
    "vst_assay.tsv"
  )
)


if (ncol(dds) <= 60L) {
  rld <- rlog(
    dds,
    blind = FALSE
  )

  write_matrix(
    assay(rld),
    file.path(
      output_dir,
      "rlog_assay.tsv"
    )
  )
}


conditions <- setdiff(
  condition_levels,
  control_level
)

control_code <- unname(
  condition_codes[[control_level]]
)


for (condition_name in conditions) {

  condition_code <- unname(
    condition_codes[[condition_name]]
  )

  comparison <- paste0(
    condition_code,
    "_vs_",
    control_code
  )

  coef_name <- paste0(
    "condition_",
    condition_code,
    "_vs_",
    control_code
  )

  if (
    !coef_name %in%
    resultsNames(dds)
  ) {
    stop(
      "Coefficient not found: ",
      coef_name
    )
  }

  res <- results(
    dds,
    contrast = c(
      "condition",
      condition_code,
      control_code
    ),
    alpha = alpha
  )

  res_shrunk <- lfcShrink(
    dds,
    coef = coef_name,
    res = res,
    type = "apeglm",
    apeMethod = ape_method
  )

  write.table(
    format_results(res),
    file = file.path(
      output_dir,
      paste0(
        comparison,
        "_raw.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    format_results(
      res_shrunk
    ),
    file = file.path(
      output_dir,
      paste0(
        comparison,
        "_lfcShrink_apeglm.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(
    condition_name,
    " vs ",
    control_level,
    ": ",
    sum(
      res$padj < alpha,
      na.rm = TRUE
    ),
    " genes with padj < ",
    alpha,
    "\n",
    sep = ""
  )
}


writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)
