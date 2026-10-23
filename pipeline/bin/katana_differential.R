#!/usr/bin/env Rscript
# Katana upstream DEG adapter. Statistical computations remain native tximport
# and DESeq2; this adapter validates identity/design and writes complete tables.
# Method references: tximport vignette and DESeq2 vignette.

fail <- function(...) stop(paste0(...), call. = FALSE)

parse_cli <- function(args) {
    opt <- list(stage = "all", design = "~ condition", outdir = ".",
                `quant-root` = ".", alpha = "0.05", `sf-type` = "ratio",
                `fit-type` = "parametric")
    allowed <- c(names(opt), "samples", "contrasts", "gtf", "gene-rds",
                 "quant-manifest", "help")
    i <- 1L
    while (i <= length(args)) {
        key <- sub("^--", "", args[[i]])
        if (!startsWith(args[[i]], "--") || !key %in% allowed)
            fail("Unknown option: ", args[[i]])
        if (key == "help") { opt$help <- TRUE; i <- i + 1L; next }
        if (i == length(args) || startsWith(args[[i + 1L]], "--"))
            fail("Missing value for --", key)
        opt[[key]] <- args[[i + 1L]]
        i <- i + 2L
    }
    if (!opt$stage %in% c("import", "normalize", "contrast", "differential", "all", "validate"))
        fail("--stage must be import, normalize, contrast, differential, all or validate")
    opt$alpha <- suppressWarnings(as.numeric(opt$alpha))
    if (length(opt$alpha) != 1L || !is.finite(opt$alpha) || opt$alpha <= 0 || opt$alpha >= 1)
        fail("--alpha must be between 0 and 1")
    if (!opt[["sf-type"]] %in% c("ratio", "poscounts"))
        fail("--sf-type must be ratio or poscounts; no automatic normalization fallback")
    if (!opt[["fit-type"]] %in% c("parametric", "local", "mean"))
        fail("--fit-type must be parametric, local or mean")
    opt
}

read_table <- function(path) {
    if (is.null(path) || !file.exists(path)) fail("Required input does not exist: ", path)
    sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
    read.table(path, sep = sep, header = TRUE, quote = "\"", comment.char = "",
               stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character",
               na.strings = character())
}

require_columns <- function(x, columns, label) {
    missing <- setdiff(columns, names(x))
    if (length(missing)) fail(label, " missing columns: ", paste(missing, collapse = ", "))
    for (column in columns)
        if (any(!nzchar(trimws(x[[column]])))) fail(label, " has blank values in ", column)
}

read_samples <- function(path) {
    x <- read_table(path)
    require_columns(x, c("sample", "condition"), "samplesheet")
    if (!nrow(x)) fail("samplesheet is empty")
    if (anyDuplicated(x$sample)) fail("Duplicate sample IDs; aggregate technical partitions before analysis")
    if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", x$sample)))
        fail("Sample IDs must be safe path components (letters, numbers, dot, underscore or hyphen)")
    rownames(x) <- x$sample
    x
}

read_contrasts <- function(path) {
    x <- read_table(path)
    require_columns(x, c("contrast", "treatment", "control"), "contrastsheet")
    if (!nrow(x)) fail("contrastsheet is empty")
    if (anyDuplicated(x$contrast)) fail("Duplicate contrast IDs")
    if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", x$contrast)))
        fail("Contrast IDs must be safe path components")
    if (any(x$treatment == x$control)) fail("Treatment and control must differ")
    x
}

validate_design <- function(samples, treatment, control, design) {
    if (!all(c(treatment, control) %in% samples$condition))
        fail("Contrast levels absent from samplesheet: ", treatment, " vs ", control)
    selected <- samples[samples$condition %in% c(control, treatment), , drop = FALSE]
    model <- tryCatch(as.formula(design), error = function(e) fail("Invalid design: ", design))
    terms <- terms(model)
    labels <- attr(terms, "term.labels")
    if (length(model) != 2L || attr(terms, "intercept") != 1L ||
        !"condition" %in% labels || any(!labels %in% c("condition", "batch", "subject")) ||
        any(attr(terms, "order") != 1L))
        fail("Supported DEG design is an additive intercept model using condition, batch and/or subject; got ", design)
    require_columns(selected, labels, "selected samplesheet")
    for (v in labels) selected[[v]] <- factor(selected[[v]])
    selected$condition <- factor(as.character(selected$condition), levels = c(control, treatment))
    for (v in setdiff(labels, "condition"))
        if (nlevels(selected[[v]]) < 2L)
            fail("Design variable ", v, " has only one level in this contrast; declare a valid reduced design explicitly")

    # Subject identity, when supplied, defines independent biological units.
    if ("subject" %in% names(selected) && any(nzchar(as.character(selected$subject))) && !all(nzchar(as.character(selected$subject))))
        fail("Subject identity is only partially specified in this contrast")
    has_subject <- "subject" %in% names(selected) && all(nzchar(as.character(selected$subject)))
    units <- if (has_subject) as.character(selected$subject) else selected$sample
    for (level in c(control, treatment))
        if (length(unique(units[selected$condition == level])) < 2L)
            fail("At least two independent biological replicates per condition are required: ", level)
    if (has_subject) {
        if (anyDuplicated(paste(units, selected$condition, sep = "\r")))
            fail("Repeated subject within one condition would create pseudoreplication; combine technical replicates first")
        crossed <- intersect(units[selected$condition == control], units[selected$condition == treatment])
        if (length(crossed) && !"subject" %in% labels)
            fail("Repeated subjects across conditions require an explicit paired design (~ subject + condition)")
        if ("subject" %in% labels && !setequal(units[selected$condition == control], units[selected$condition == treatment]))
            fail("Paired DEG design requires complete treatment/control subject pairs for this contrast")
    }
    mat <- model.matrix(model, selected)
    if (qr(mat)$rank != ncol(mat))
        fail("Inestimable/confounded DEG design for ", treatment, " vs ", control, ": ", design)
    if (nrow(mat) <= ncol(mat)) fail("DEG model has no residual degrees of freedom: ", design)
    list(samples = selected, formula = model, matrix = mat, design = design,
         biological_units = units, paired = "subject" %in% labels)
}

validate_contrasts <- function(samples, contrasts, design) {
    lapply(seq_len(nrow(contrasts)), function(i) {
        declared <- if ("design" %in% names(contrasts) && nzchar(contrasts$design[[i]])) contrasts$design[[i]] else design
        tryCatch(validate_design(samples, contrasts$treatment[[i]], contrasts$control[[i]], declared),
                 error = function(e) fail("Contrast ", contrasts$contrast[[i]], ": ", conditionMessage(e)))
    })
}

gtf_attribute <- function(attributes, name) {
    pattern <- paste0("(?:^|;[[:space:]]*)", name, "[[:space:]]+\"([^\"]*)\"")
    match <- regexec(pattern, attributes, perl = TRUE)
    values <- regmatches(attributes, match)
    vapply(values, function(x) if (length(x) >= 2L) x[[2L]] else NA_character_, character(1L))
}

read_gtf_mapping <- function(path) {
    if (is.null(path) || !file.exists(path)) fail("GTF does not exist: ", path)
    con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con))
    tx_chunks <- list(); gene_chunks <- list(); index <- 0L
    repeat {
        lines <- readLines(con, n = 50000L, warn = FALSE)
        if (!length(lines)) break
        lines <- lines[!startsWith(lines, "#") & nzchar(lines)]
        if (!length(lines)) next
        fields <- strsplit(lines, "\t", fixed = TRUE)
        if (any(lengths(fields) != 9L)) fail("Malformed GTF: expected nine tab-separated fields")
        table <- do.call(rbind, fields)
        table <- table[table[, 3] %in% c("gene", "transcript", "exon"), , drop = FALSE]
        if (!nrow(table)) next
        gid <- gtf_attribute(table[, 9], "gene_id")
        tid <- gtf_attribute(table[, 9], "transcript_id")
        gname <- gtf_attribute(table[, 9], "gene_name")
        biotype <- gtf_attribute(table[, 9], "gene_type")
        alternate <- gtf_attribute(table[, 9], "gene_biotype")
        biotype[is.na(biotype)] <- alternate[is.na(biotype)]
        valid <- !is.na(gid) & nzchar(gid)
        index <- index + 1L
        tx_chunks[[index]] <- unique(data.frame(transcript_id = tid[valid & !is.na(tid) & nzchar(tid)],
            gene_id = gid[valid & !is.na(tid) & nzchar(tid)], stringsAsFactors = FALSE))
        gene_chunks[[index]] <- unique(data.frame(gene_id = gid[valid], gene_name = gname[valid],
            gene_biotype = biotype[valid], stringsAsFactors = FALSE))
    }
    if (!length(tx_chunks)) fail("No gene/transcript mapping found in GTF")
    tx2gene <- unique(do.call(rbind, tx_chunks))
    if (!nrow(tx2gene)) fail("No transcript_id/gene_id mapping found in GTF")
    if (anyDuplicated(tx2gene$transcript_id)) fail("GTF maps a transcript ID to multiple gene IDs")
    genes <- unique(do.call(rbind, gene_chunks))
    ids <- unique(genes$gene_id)
    annotation <- data.frame(gene_id = ids, gene_name = NA_character_, gene_biotype = NA_character_,
                             stringsAsFactors = FALSE)
    for (column in c("gene_name", "gene_biotype")) {
        nonmissing <- genes[!is.na(genes[[column]]) & nzchar(genes[[column]]), c("gene_id", column), drop = FALSE]
        nonmissing <- unique(nonmissing)
        if (anyDuplicated(nonmissing$gene_id)) fail("Conflicting GTF ", column, " values for one gene_id")
        annotation[[column]] <- nonmissing[[column]][match(ids, nonmissing$gene_id)]
    }
    list(tx2gene = tx2gene, annotation = annotation)
}

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE,
                                          row.names = FALSE, na = "NA")
write_matrix <- function(x, path) write_tsv(data.frame(gene_id = rownames(x), x, check.names = FALSE), path)
write_kv <- function(values, path) write_tsv(data.frame(key = names(values), value = unlist(values),
                                                       stringsAsFactors = FALSE), path)
write_versions <- function(packages, path, process) {
    values <- c(R = paste(R.version$major, R.version$minor, sep = "."),
                setNames(vapply(packages, function(p) as.character(packageVersion(p)), character(1L)), packages))
    writeLines(c(paste0(process, ":"), paste0("  ", names(values), ": \"", values, "\"")), path)
}
require_package <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) fail("Missing R dependency: ", package)
}

quant_files <- function(opt, samples) {
    if (!is.null(opt[["quant-manifest"]])) {
        manifest <- read_table(opt[["quant-manifest"]])
        require_columns(manifest, c("sample", "quant"), "quant manifest")
        if (anyDuplicated(manifest$sample) || !setequal(manifest$sample, samples$sample))
            fail("Quant manifest must match samplesheet exactly once per sample")
        files <- manifest$quant[match(samples$sample, manifest$sample)]
    } else files <- file.path(opt[["quant-root"]], samples$sample, "quant.sf")
    if (any(!file.exists(files))) fail("Missing quant.sf for: ", paste(samples$sample[!file.exists(files)], collapse = ", "))
    names(files) <- samples$sample
    files
}

import_gene_quant <- function(opt, samples) {
    require_package("tximport")
    files <- quant_files(opt, samples)
    reference <- read_gtf_mapping(opt$gtf)
    transcript_ids <- NULL
    for (sample in names(files)) {
        q <- read_table(files[[sample]])
        require_columns(q, c("Name", "Length", "EffectiveLength", "TPM", "NumReads"), paste0("quant.sf ", sample))
        if (anyDuplicated(q$Name)) fail("Duplicate quantified transcript IDs in ", sample)
        for (column in c("Length", "EffectiveLength", "TPM", "NumReads")) {
            vals <- suppressWarnings(as.numeric(q[[column]]))
            if (any(!is.finite(vals)) || any(vals < 0)) fail("Invalid ", column, " in ", sample)
        }
        if (is.null(transcript_ids)) transcript_ids <- q$Name
        else if (!identical(transcript_ids, q$Name)) fail("Quantified transcript IDs/order differ across samples: ", sample)
        if (sum(as.numeric(q$NumReads)) <= 0) fail("No quantified fragments in ", sample)
    }
    missing <- setdiff(transcript_ids, reference$tx2gene$transcript_id)
    coverage <- data.frame(metric = c("quantified_transcripts", "mapped_transcripts", "unmapped_transcripts", "coverage"),
        value = c(length(transcript_ids), length(transcript_ids) - length(missing), length(missing),
                  (length(transcript_ids) - length(missing)) / length(transcript_ids)))
    write_tsv(coverage, file.path(opt$outdir, "tx2gene_coverage.tsv"))
    if (length(missing)) {
        write_tsv(data.frame(transcript_id = missing), file.path(opt$outdir, "unmapped_transcripts.tsv"))
        fail("Quantified transcript IDs missing from matching GTF: ", length(missing), "; no fabricated gene mappings or silent ID stripping")
    }
    mapping <- reference$tx2gene[match(transcript_ids, reference$tx2gene$transcript_id), , drop = FALSE]
    txi <- tximport::tximport(files, type = "salmon", tx2gene = mapping,
                             countsFromAbundance = "no", ignoreTxVersion = FALSE,
                             ignoreAfterBar = FALSE, dropInfReps = TRUE)
    if (!identical(colnames(txi$counts), samples$sample)) fail("tximport sample order changed unexpectedly")
    if (any(!is.finite(txi$length)) || any(txi$length <= 0)) fail("Gene effective lengths must be finite and positive")
    if (any(!is.finite(txi$counts)) || any(txi$counts < 0)) fail("Invalid estimated gene counts")
    saveRDS(txi, file.path(opt$outdir, "gene_txi.rds"), version = 2)
    write_matrix(txi$counts, file.path(opt$outdir, "gene.estimated_counts.tsv"))
    write_matrix(txi$abundance, file.path(opt$outdir, "gene.tpm.tsv"))
    write_matrix(txi$length, file.path(opt$outdir, "gene.effective_lengths.tsv"))
    write_tsv(mapping, file.path(opt$outdir, "tx2gene.tsv"))
    annotation <- reference$annotation[match(rownames(txi$counts), reference$annotation$gene_id), , drop = FALSE]
    write_tsv(annotation, file.path(opt$outdir, "annotation.tsv"))
    write_kv(list(stage = "import", counts_representation = "Salmon estimated fragment counts (potentially fractional)",
                  countsFromAbundance = "no", transcript_ids = "exact; no version stripping",
                  gene_level_artifacts_created = "true", samples = paste(samples$sample, collapse = ","),
                  gtf = basename(opt$gtf)), file.path(opt$outdir, "import_metadata.tsv"))
    write_versions("tximport", file.path(opt$outdir, "statistics.import.versions.yml"), "KATANA_TXIMPORT")
    capture.output(sessionInfo(), file = file.path(opt$outdir, "import.R_sessionInfo.log"))
    txi
}

validate_txi <- function(txi, samples) {
    if (!is.list(txi) || !all(c("counts", "abundance", "length", "countsFromAbundance") %in% names(txi)))
        fail("gene-rds is not a tximport list")
    if (!identical(txi$countsFromAbundance, "no")) fail("DEG expects countsFromAbundance=no with effective-length offsets")
    if (!identical(colnames(txi$counts), samples$sample)) fail("gene-rds sample IDs/order must match samplesheet exactly")
    for (name in c("counts", "abundance", "length")) {
        x <- txi[[name]]
        if (!is.matrix(x) || !identical(dimnames(x), dimnames(txi$counts)) || any(!is.finite(x)) || any(x < 0))
            fail("Invalid tximport ", name, " matrix")
    }
    if (any(txi$length <= 0)) fail("Gene effective lengths must be positive")
    if (is.null(rownames(txi$counts)) || anyDuplicated(rownames(txi$counts))) fail("Gene IDs must be unique")
}

run_normalization <- function(opt, samples, txi) {
    require_package("DESeq2")
    validate_txi(txi, samples)
    # This is a cohort-wide expression matrix, not a differential test. Native
    # DESeq2 normalization uses average transcript lengths supplied by tximport.
    all_dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = samples, design = ~ 1)
    all_dds <- DESeq2::estimateSizeFactors(all_dds, type = opt[["sf-type"]])
    write_matrix(DESeq2::counts(all_dds, normalized = TRUE), file.path(opt$outdir, "all.normalised_counts.tsv"))
    write_matrix(DESeq2::normalizationFactors(all_dds), file.path(opt$outdir, "all.normalization_factors.tsv"))
    write_matrix(DESeq2::counts(all_dds), file.path(opt$outdir, "gene.deseq2_integer_counts.tsv"))
    write_versions("DESeq2", file.path(opt$outdir, "statistics.normalization.versions.yml"), "KATANA_NORMALIZE")
}

run_differential <- function(opt, samples, contrasts, models, txi) {
    require_package("DESeq2")
    validate_txi(txi, samples)
    per_contrast <- identical(opt$stage, "contrast")
    if (per_contrast && nrow(contrasts) != 1L) fail("Exactly one comparison required for --stage contrast")
    if (!per_contrast) run_normalization(opt, samples, txi)
    records <- list()
    for (i in seq_len(nrow(contrasts))) {
        record <- contrasts[i, , drop = FALSE]
        model <- models[[i]]
        ids <- model$samples$sample
        subset_txi <- txi
        for (name in c("counts", "abundance", "length")) subset_txi[[name]] <- txi[[name]][, ids, drop = FALSE]
        dds <- DESeq2::DESeqDataSetFromTximport(subset_txi, colData = model$samples, design = model$formula)
        dds <- DESeq2::DESeq(dds, test = "Wald", fitType = opt[["fit-type"]],
                              sfType = opt[["sf-type"]], minReplicatesForReplace = Inf,
                              parallel = FALSE)
        res <- DESeq2::results(dds, contrast = c("condition", record$treatment, record$control),
                               alpha = opt$alpha, independentFiltering = TRUE, pAdjustMethod = "BH")
        table <- data.frame(gene_id = rownames(res), as.data.frame(res), check.names = FALSE)
        table$contrast <- record$contrast
        table$treatment <- record$treatment
        table$control <- record$control
        table$direction <- "log2(treatment/control)"
        table$test_status <- ifelse(table$baseMean == 0, "all_zero",
                                   ifelse(is.na(table$pvalue), "not_tested_or_outlier",
                                          ifelse(is.na(table$padj), "independent_filter", "tested")))
        prefix <- file.path(opt$outdir, record$contrast)
        write_tsv(table, paste0(prefix, ".deseq2.results.tsv"))
        write_matrix(DESeq2::counts(dds, normalized = TRUE), paste0(prefix, ".normalised_counts.tsv"))
        saveRDS(dds, paste0(prefix, ".dds.rds"), version = 2)
        writeLines(c(model$design, paste0("effect: log2(", record$treatment, "/", record$control, ")"),
                     paste0("samples: ", paste(ids, collapse = ",")),
                     paste0("coefficients: ", paste(DESeq2::resultsNames(dds), collapse = ","))),
                   paste0(prefix, ".deseq2.model.txt"))
        records[[i]] <- data.frame(contrast = record$contrast, treatment = record$treatment, control = record$control,
            design = model$design, paired = model$paired, sample_count = length(ids),
            control_biological_replicates = length(unique(model$biological_units[model$samples$condition == record$control])),
            treatment_biological_replicates = length(unique(model$biological_units[model$samples$condition == record$treatment])),
            gene_count = nrow(table), tested_pvalues = sum(!is.na(table$pvalue)), reported_padj = sum(!is.na(table$padj)),
            effect_direction = "log2(treatment/control)", stringsAsFactors = FALSE)
    }
    suffix <- if (per_contrast) paste0(".", contrasts$contrast[[1]]) else ""
    write_tsv(do.call(rbind, records), file.path(opt$outdir, paste0("models", suffix, ".tsv")))
    write_kv(list(stage = "differential", test = "native DESeq2 Wald", alpha = opt$alpha,
                  independentFiltering = "true", pAdjustMethod = "BH", lfc_shrinkage = "none; native MLE retained",
                  automatic_count_replacement = "disabled (minReplicatesForReplace=Inf)",
                  countsFromAbundance = "no", effective_length_normalization = "DESeqDataSetFromTximport",
                  all_normalised_counts = "cohort-wide ~1 normalization; not a differential test",
                  fit_type_requested = opt[["fit-type"]], sf_type = opt[["sf-type"]],
                  contrast_scope = "each declared treatment/control subset", genes_filtered_before_model = "none",
                  samples = paste(samples$sample, collapse = ",")), file.path(opt$outdir, paste0("differential_metadata", suffix, ".tsv")))
    write_versions("DESeq2", file.path(opt$outdir, paste0("statistics.differential", suffix, ".versions.yml")), "KATANA_DESEQ2")
    capture.output(sessionInfo(), file = file.path(opt$outdir, paste0("differential", suffix, ".R_sessionInfo.log")))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
    opt <- parse_cli(args)
    if (isTRUE(opt$help)) {
        cat("Usage: katana_differential.R --stage import|normalize|contrast|differential|all|validate --samples samples.csv\n",
            "Import: --gtf annotation.gtf --quant-root DIR [or --quant-manifest sample_quant.tsv]\n",
            "Differential: --gene-rds gene_txi.rds --contrasts contrasts.csv --design '~ condition'\n",
            "All combines import and differential. Outputs: --outdir DIR\n",
            "Optional: --alpha 0.05 --fit-type parametric --sf-type ratio\n", sep = "")
        return(invisible(NULL))
    }
    samples <- read_samples(opt$samples)
    contrasts <- models <- NULL
    if (!opt$stage %in% c("import", "normalize")) {
        contrasts <- read_contrasts(opt$contrasts)
        models <- validate_contrasts(samples, contrasts, opt$design)
    }
    if (opt$stage == "validate") {
        cat("VALID_DESIGN: ", length(models), " declared contrasts\n", sep = "")
        return(invisible(models))
    }
    if (!dir.exists(opt$outdir) && !dir.create(opt$outdir, recursive = TRUE)) fail("Cannot create output directory")
    txi <- if (opt$stage %in% c("import", "all")) import_gene_quant(opt, samples) else {
        if (is.null(opt[["gene-rds"]]) || !file.exists(opt[["gene-rds"]])) fail("Missing --gene-rds")
        readRDS(opt[["gene-rds"]])
    }
    if (opt$stage == "normalize") run_normalization(opt, samples, txi)
    if (opt$stage %in% c("contrast", "differential", "all")) run_differential(opt, samples, contrasts, models, txi)
    invisible(NULL)
}

if (sys.nframe() == 0L) main()
