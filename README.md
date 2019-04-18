# Deep Drug Response

This repository contains the code to train the drug-response models used in
*"A machine learning framework for mining mechanistic insights and predicting
response to therapy in cancer"* by *Vougas, Sakellaropoulos et al.*

The main script `tune_drug.R` is run independently for each drug. The man page
of the script is:

```
Run training pipeline

positional arguments:
  drug            Drug to learn
  matrix          Matrix to train on (RDS file)
  hyperparams     json file with hyper parameters
  prefix          prefix to write files

flags:
  -h, --help                     show this help message and exit

optional arguments:
  -x, --opts OPTS                RDS file containing argument values
  -m, --max_models MAX_MODELS    Max number of models to compute [default: 15]
  -v, --varcut VARCUT            Percent of genes to keep based on variance. [default: 100]
  -l, --learner LEARNER          Algorithm to use for learning (DL, RF, GLM) [default: DL]
  -t, --test TEST                Path to clinical data
  -k, --kfold KFOLD              Number of CV-splits [default: 5]
  --max_mem MAX_MEM              H2O parameter [default: 64G]
  -c, --cores CORES              Method to determine cores with future::availableCores [default: Slurm]
```

`drug` must be one of the columns of `matrix` the rest are treated as genes.
The hyperparameters used in the paper are specifed in `hyperparams/`.
All of them are specified as lists in order to be able to expand the grid.

The script `batch_correct.R` is used to correct for batch effects between the training and
test dataset if required.
