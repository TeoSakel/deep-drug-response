#!/usr/bin/env Rscript
# inputs    -------------------------------------------------------------------
library(argparser)
library(pROC)

p = arg_parser("Run training pipeline")
p = add_argument(p, "drug", help = "Drug to learn")
p = add_argument(p, "matrix", help = "Matrix to train on (RDS file)")
p = add_argument(p, "hyperparams", help = "json file with hyper parameters")
p = add_argument(p, "prefix", help = "prefix to write files")
p = add_argument(p, "--max_models", default = 15L,
                 help = "Max number of models to compute")
p = add_argument(p, "--FSmethod", default = "mad",
                 help = "Feature Selection method to use {mad|wilcox}")
p = add_argument(p, "--FSthres", default = 1,
                 help = "Threshold for feature selection")
p = add_argument(p, "--learner", default = "DL",
                 help = "Algorithm to use for learning (DL, RF, GLM)")
p = add_argument(p, "--test", help = "Path to clinical data")
p = add_argument(p, "--kfold", default = 5L, help = "Number of CV-splits")
p = add_argument(p, "--max_mem", default = "64G", help = "H2O parameter")
p = add_argument(p, "--cores", default = "Slurm",
                 help = "Method to determine cores with future::availableCores")

argv = parse_args(p)
# argv = parse_args(p, c("Bortezomib", "bortezomib_train.rds", "hyperparam_grid_dl.json",
#                        "debug", "--max_models", "2", "--varcut", "50"))

message("Running with the following arguments")
str(argv)

drug = argv$drug
mat_file = argv$matrix
lift_null_file = argv$lift
hyperparam_file = argv$hyperparams
prefix = argv$prefix
Nmodels = argv$max_models
FSmethod = match.arg(argv$FSmethod, c("mad", "wilcox"))
FSthres = argv$FSthres
learner = match.arg(argv$learner, c("DL", "RF", "GLM"))
kfold = argv$kfold
max_mem = argv$max_mem
test_file = argv$test
cores_method = argv$cores

Ncores = future::availableCores(method = cores_method)
options(mc.cores = Ncores)
cat('*****************\n')
cat(learner,' : ' ,mat_file,'\n')
cat('*****************\n')

# Library    -------------------------------------------------------------------

suppressPackageStartupMessages({
    library(tidyverse)
    library(parallel)
    library(glmnet)
    library(glmnetUtils)
    library(h2o)
	library(pROC)
    library(caret)
})

all_equal = function(x, y) {
    if (length(y) == 1)
        y = rep(y, length(x))
    !is.character(all.equal(x, y, check.attributes = FALSE))
}

is_scaled = function(mat) {
    # result of `stat::scale`
    if (!is.null(attr(mat, "scaled:center")) )
        return(TRUE)
    # zero-mean
    cmeans = colMeans(mat, na.rm = TRUE)
    if (!all_equal(cmeans, 0))
        return(FALSE)
    # unit variance
    csd = apply(mat, 2, sd)
    if (!all_equal(csd, 1))
        return(FALSE)

    return(TRUE)
}

rbind_data_frame = function(...)
    # simple wrapper
    rbind.data.frame(..., make.row.names = FALSE, stringsAsFactors = FALSE)

delist_col = function(x, collapse = '|') {
    if (!is.list(x))
        return(x)
    vapply(x, paste, '', sep = '', collapse = collapse)
}

delist_df = function(x, collapse = '|')
    as.data.frame(lapply(x, delist_col, collapse = collapse))

FS_top = function(mat, top) {
    top = floor(top * length(genes))
    ix = apply(mat[, genes], 2, mad, na.rm = TRUE)
    ix = order(ix, decreasing = TRUE)[1:top]
    genes = genes[ix]
    if (!is.na(genes))
    c(tissues,genes)
	else
	tissues
}

FS_wilcox = function(mat, thres) {
    pvals = vapply(genes, function(gene) {
        frml = as.formula(sprintf('`%s` ~ %s', gene, drug))
        suppressWarnings(wilcox.test(frml, data = mat)$p.value)
    }, 0, USE.NAMES = FALSE)
    pvals = p.adjust(pvals, method = "fdr")
    genes = genes[pvals < thres]
    if (length(genes) > 0) c(tissues, genes) else tissues
}

expand_hyperparams_h2o = function(hyper_params) {
    hyper_params = split(hyper_params, seq(nrow(hyper_params)))
    hyper_params = lapply(hyper_params, as.list)
    names(hyper_params) = NULL
    for (i in seq(length(hyper_params))) {
        for (j in seq(length(hyper_params[[i]]))) {
            if (is.list(hyper_params[[i]][[j]]))
                hyper_params[[i]][[j]] = unlist(hyper_params[[i]][[j]],
                                                recursive = FALSE,
                                                use.names = FALSE)
        }
    }
    hyper_params
}

validate_h2o = function(model, newdata) {
    y_hat = h2o.performance(model, newdata)
    data.frame(MSE = h2o.mse(y_hat), AUC = h2o.auc(y_hat))
}

validate_glmnet = function(model, newdata, s = 'lambda.min') {
  idrug = which(colnames(newdata) == drug)
  #//y_hat = predict(model, as(newdata[, -idrug],'dgCMatrix'), s = s)
  y_hat = predict(model, newdata[, -idrug], s = s)
  dy = y_hat - newdata[, idrug]
  d1=roc(as.numeric(newdata[,idrug]),as.numeric(y_hat))
  data.frame(lambda = model[[s]],  # hacky
             MSE = mean(dy * dy), AUC =d1$auc)
}

# Preprocessing -----------------------------------------------------------------

message(paste0('Sourcing ', hyperparam_file, "..."))
hyper_params <- jsonlite::fromJSON(hyperparam_file,
                                   simplifyVector = TRUE,
                                   simplifyDataFrame = FALSE,
                                   simplifyMatrix = FALSE)

if (FSmethod == "mad") {
    FS = FS_top
} else if (FSmethod == "wilcox") {
    FS = FS_wilcox
} else {
    stop("Unknown FSmethod")
}

message(paste0('Reading ', mat_file, '...'))
#mat = as.matrix(readRDS(mat_file))
mat = readRDS(mat_file)
mat$recurrence=as.factor(mat$recurrence) ##hua zhou 2020.03.13
tissues = grep("^Tissue_", colnames(mat), value = TRUE)
if (any(grepl("_GE", colnames(mat)))) {
    genes = grep("_GE$", colnames(mat), value = TRUE)
} else {
    genes = setdiff(colnames(mat), c(drug, tissues))
}

message('Subseting matrix...')
mat = mat[!is.na(mat[, drug]),  ]

# k-fold cross-validation
set.seed(123)
mat = mat[sample(nrow(mat)), ]
fold_id = cut(1:nrow(mat), breaks = kfold, labels = FALSE)
flds <- createFolds(mat[,1], k = kfold, list = TRUE, returnTrain = FALSE)
for (i in 1:kfold) {
  fold_id[flds[[i]]]=i
}
fold_size = sum(fold_id != 1L)

# message(sprintf('Selecting genes based on top %.1f%% variance cutoff', FSthres))
features = mclapply(1:kfold, function(i) FS(mat[fold_id != i, ], FSthres))
features[[kfold + 1L]] = FS(mat, FSthres)
cat(features[[1]])

message('Setting up learning algorithm...')
if (learner == "DL") {
    h2o.init(nthreads = Ncores, max_mem_size = max_mem)
    convert = as.h2o
    learn = h2o.deeplearning
    validate = validate_h2o
    predict = h2o.predict
    hyper_params_df = lapply(hyper_params, expand.grid,
                             KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    hyper_params_df = Reduce(rbind_data_frame, hyper_params_df)
    Nmodels = min(Nmodels, nrow(hyper_params_df))
    hyper_params_df = hyper_params_df[sample(nrow(hyper_params_df), Nmodels), ]
    hyper_params = expand_hyperparams_h2o(hyper_params_df)
    hyper_params_df = cbind(data.frame(model_id = 1:Nmodels), hyper_params_df)
    saveModel = function(x, nam) h2o.saveModel(x, paste0(nam, ".h2o"))
} else if (learner == "RF") {
    h2o.init(nthreads = Ncores, max_mem_size = max_mem)
    convert = as.h2o
    learn = h2o.randomForest
    validate = validate_h2o
    predict = h2o.predict
    hyper_params_df = lapply(hyper_params, expand.grid,
                             KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    hyper_params_df = Reduce(rbind_data_frame, hyper_params_df)
    Nmodels = min(Nmodels, nrow(hyper_params_df))
    hyper_params_df = hyper_params_df[sample(nrow(hyper_params_df), Nmodels), ]
    hyper_params = expand_hyperparams_h2o(hyper_params_df)
    hyper_params_df = cbind(data.frame(model_id = 1:Nmodels), hyper_params_df)
    saveModel = function(x, nam) h2o.saveModel(x, paste0(nam, ".h2o"))
} else {
    # GLM
    validate = validate_glmnet
    convert = function(x, ...) {
		x[,1]=(x[,1]==1)+0
			as.matrix(x)}
# learn = function(training_frame, y, x, alpha) {
#        cv.glmnet(training_frame[, x], training_frame[, y],
#                  parallel = TRUE, alpha = alpha, family = "gaussian")
#    }
    learn = function(training_frame, y, x, alpha) {
        cv.glmnet(training_frame[, x], as.vector(training_frame[, y]),
                  parallel = TRUE, alpha = alpha, family = "binomial") }
    hyper_params = hyper_params
    hyper_params_df = data.frame(model_id = 1:length(hyper_params),
                                 alpha = hyper_params)
    saveModel = function(x, nam) saveRDS(x, paste0(nam, ".rds"))
    Nmodels = length(hyper_params)
}

message('Running cross-validation')
cvSummary = vector("list", kfold)
cvSummary2 = vector("list", kfold)
for (i in seq_len(kfold)) {
    kmat = mat[, c(drug, features[[i]]), drop = FALSE]
    ix = which(fold_id == i)
    train_hex = convert(kmat[-ix, , drop = FALSE])
    valid_hex = convert(kmat[ix, , drop = FALSE])
 if (length(features[[i]])>1 || learner!='GLM') {
    kmodel = lapply(hyper_params, function(params) {
                    x = features[[i]]
                    if ('mtries' %in% names(params)) {# special case
                        params[['mtries']] = round(params[['mtries']] * length(x))
                    }
                    do.call(learn, c(list(y = drug, x = x,
                                          training_frame = train_hex),
                                     params))
                  })
    kvalid = lapply(kmodel, validate, newdata = valid_hex)
    kvalid = do.call(rbind, kvalid)
    tmp_df = cbind(data.frame(fold_id = rep(i, Nmodels)),
                   hyper_params_df,
                   kvalid)
    rownames(tmp_df) = NULL
    cvSummary[[i]] = tmp_df } else {
    kmodel=glm(as.formula(paste0(drug,'~',features[[i]])),family='binomial',data=data.frame(train_hex))
    cauc=auc(valid_hex[,1],predict(kmodel,newdata=data.frame(valid_hex)))

    kvalid = data.frame(model_id=1,MSE=0,AUC=cauc)
    tmp_df = cbind(data.frame(fold_id = rep(i, 1)),
                   hyper_params_df[1,],
                   kvalid)
    rownames(tmp_df) = NULL
    cvSummary[[i]] = kvalid }

}

# Write results -----------------------------------------------------------
cvSummary = do.call(rbind, cvSummary)
cvSummary = delist_df(cvSummary)
cvSummary = cvSummary[order(cvSummary[["AUC"]],decreasing = T), ]
write_csv(cvSummary, paste0(prefix, "_cv.csv"))

cvSummary = cvSummary %>%
  group_by_at(vars(model_id)) %>%
  summarize(MSE_sd = sd(MSE, na.rm = TRUE),
            MSE = mean(MSE, na.rm = TRUE),
            AUC_sd = sd(AUC, na.rm = TRUE),
            AUC = mean(AUC, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(AUC)
write_csv(cvSummary, paste0(prefix, "_cvSummary.csv"))



# Train best model on all of the data
best_model = cvSummary$model_id[nrow(cvSummary)]
train_hex = convert(mat)
params = hyper_params[[best_model]]
x = features[[kfold + 1L]]
if ('mtries' %in% names(params))
    params[['mtries']] = floor(params[['mtries']] * length(x))
  if (length(x)>1 || learner != 'GLM')
model = do.call(learn, c(list(y = drug, x = x, training_frame = train_hex), params)) else
model=glm(as.formula(paste0(drug,'~',x)),family='binomial',data=data.frame(train_hex))
saveModel(model, prefix)




best_model = cvSummary$model_id[nrow(cvSummary)]
###start getting predictions
validate_h2o = function(model, newdata) {
if (learner!='GLM') {
  predict=h2o.predict
  y_hat = predict(model, newdata)
  #data.frame(MSE = h2o.mse(y_hat), AUC=h2o.auc(y_hat))
  y_hat=as.data.frame(y_hat)
  data.frame(y=as.data.frame(newdata[,1]),ypred=y_hat$predict,yprob=y_hat$p1)}
  else {
    y_hat = predict(model, newdata[,-1],s='lambda.min',type='response')
    #data.frame(MSE = h2o.mse(y_hat), AUC=h2o.auc(y_hat))
    y_hat=as.data.frame(y_hat)
    data.frame(y=newdata[,1],ypred=y_hat[,1],yprob=y_hat[,1])
  }

}

validate = validate_h2o
message('Running cross-validation')
cvSummary = vector("list", kfold)
subjs=c()
for (i in seq_len(kfold)) {
  kmat = mat[, c(drug, features[[i]]), drop = FALSE]
  ix = which(fold_id == i)
  train_hex = convert(kmat[-ix, , drop = FALSE])
  valid_hex = convert(kmat[ix, , drop = FALSE])
  subjs=c(subjs,rownames(kmat[ix, , drop = FALSE]))
if (length(features[[i]])>1) {
  kmodel = lapply(hyper_params, function(params) {
    x = features[[i]]
    if ('mtries' %in% names(params)) {# special case
      params[['mtries']] = round(params[['mtries']] * length(x))
    }
    do.call(learn, c(list(y = drug, x = x,
                          training_frame = train_hex),
                     params))
  })
  kvalid =validate(kmodel[[best_model]],valid_hex)
  cvSummary[[i]] = kvalid }  else {
 kmodel=glm(as.formula(paste0(drug,'~',features[[i]])),family='binomial',data=data.frame(train_hex))
    cauc=predict(kmodel,newdata=data.frame(valid_hex))
    cvSummary[[i]] = data.frame(recurrence=valid_hex[,1],pred=cauc)
}
}

cvSummary = do.call(rbind, cvSummary)
cvSummary2=data.frame(subjs=subjs,cvSummary)

write_csv(cvSummary2, paste0(prefix, "_pred_prob.csv"))



if (!is.na(test_file)) {
    test = readRDS(test_file)
    patient_id = rownames(test)
    test = convert(test[, x, drop = FALSE])
    y = as.vector(predict(model, test))
    df = tibble(patient_id = patient_id, ic50 = y)
    write_csv(df, paste0(prefix, "_clinicalPrediction.csv"))
}

if (learner %in% c("DL", "RF"))
    h2o.shutdown(FALSE)



