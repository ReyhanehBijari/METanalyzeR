This README file provides an overview of the datasets stored in the `Data` directory of the METanalyzeR repository.

## 1. Soybean Yield Data from Multi-Environment Trials

### Description
This dataset represents the yield outcomes from two multi-environment experiments conducted as part of a soybean breeding program in the United States.

### Format
The dataset is structured as a data frame comprising 2000 observations with details across 5 variables:
- **Experiment**: Identifier for the experiment, with two distinct levels.
- **Genotype**: Categorizes the 46 soybean varieties. It's important to note that Genotype G1 in Experiment 1 is distinct from G2 in Experiment 2.
- **Environment**: Identifier for the 44 different environments where the trials were conducted.
- **Replication (rep)**: Indicates the replication factor, with two levels.
- **Yield**: Measures the yield in bushels per acre.

### Details
Breeding programs have complex decision-making processes, often stretching over several years from initial experimental stages to the commercialization of new genotypes. This dataset represents soybean data from two late-stage experiments in such a breeding program.

## 2. shafii.rapeseed: Multi-environment Trial of Rapeseed

This dataset is part of the "agridat" package on CRAN, which is a comprehensive resource offering datasets from various agricultural experiments. More information can be found on the [shafii.rapeseed dataset page](https://rdrr.io/cran/agridat/man/shafii.rapeseed.html).
