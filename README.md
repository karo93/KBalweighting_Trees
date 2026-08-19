# KBal Weighting for Tree-Based Methods

Replication materials for Kernel Balancing in Tree-Based Methods, which integrates kernel balancing (KBal) into the causal forest framework (Wager and Athey, 2018) and the X-learner framework (Künzel et al., 2019) for estimating heterogeneous treatment effects.

# Overview

Reliable estimation of conditional average treatment effects (CATEs) requires sufficient overlap between treated and control units. When overlap is weak, conventional propensity-score-based weighting methods may become unstable due to extreme weights, finite-sample bias, and limited covariate balance.

This repository investigates kernel balancing (KBal) as an alternative weighting strategy for tree-based causal machine learning methods. Rather than relying on estimated propensity scores, KBal constructs weights that balance treated and control groups in a transformed feature space. These weights are incorporated into two popular CATE estimators:

* KBal-weighted Causal Forests, where kernel balancing weights are used as sample weights during forest estimation.
* KBal-weighted X-Learner (XRF), where kernel balancing replaces the propensity score in the treatment-effect aggregation step.

The repository contains the simulation study evaluating estimator performance under varying degrees of overlap, as well as a semi-synthetic application based on the IHDP benchmark dataset.

## Repository structure

```text
KBalweighting_Trees/
├── README.md
├── LICENSE
├── .gitignore
├── simulation/
│   ├── data/                  simulation datasets
│   ├── raw_results/           simulation output (generated)
│   ├── aggregated_results/    summary statistics (generated)
│   └── *.R                    simulation scripts
├── application/
│   ├── data/                  IHDP benchmark datasets
│   ├── raw_results/           application output (generated)
│   └── *.R                    empirical application scripts
└── output/
    ├── sim_plots/
    ├── sim_tables/
    ├── application_plots/
    └── application_tables/
```
