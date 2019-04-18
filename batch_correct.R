library(argparser)
p = arg_parser("Normalize training and test set")
p = add_argument(p, "train", help = "training matrix")
p = add_argument(p, "test", help = "testing matrix")
p = add_argument(p, "prefix", help = "prefix to write files")

argv = parse_args(p)
argv = parse_args(p, c('data/raw/cells.rds', 'data/raw/bortezomib.rds', 'debug'))

cell = readRDS(argv$train)
clin = readRDS(argv$test)


cell = cell[, !apply(is.na(cell), 2L, any)]
clin = clin[, !apply(is.na(clin), 2L, any)]

genes = intersect(colnames(cell), colnames(clin))
stopifnot(length(genes) > 0L)
cell = cell[, genes]
clin = clin[, genes]

suppressPackageStartupMessages(library(sva))

batches = rep(c('cell', 'clin'), times = c(nrow(cell), nrow(clin)))
batches = as.factor(batches)

cell_clin = t(rbind(cell, clin))
cell_clin = ComBat(cell_clin, batches)
cell_clin = t(cell_clin)

cell = cell_clin[batches == 'cell', ]
clin = cell_clin[batches == 'clin', ]

saveRDS(cell, paste0(argv$prefix, "cells.rds"))
saveRDS(clin, paste0(argv$prefix, "clinical.rds"))
