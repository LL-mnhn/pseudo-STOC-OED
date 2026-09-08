# pseudo-STOC-OED

## Overview
This repository serves as preparation for works on **Optimal Experimental Designs (OEDs)** for the BIODICAPT project. It focuses on the use of **STOC**: a large dataset of bird sightings, to setup the algorithms and workflow.

> [!TIP]
> A quick analysis of our results is available [**here**](https://ll-mnhn.github.io/pseudo-STOC-OED/report/report.html).

## Context
### BIODICAPT
BIODICAPT is a French initiative that aims at monitoring biodiversity of agricultural lands on a large (national) scale through the use of various recording devices.

BIODICAPT started in early-2026, the first results (data extraction of species distributions) will not be available until late-2026 or early-2027. In the meantime, we use **STOC** to prepare our workflow.

### STOC
The French Breeding Bird Survey (FBBS, or STOC in french) is a standardised multi-species abundance dataset, based on a citizen science program, which has been running since 1989.

It contains bird sightings and explanatory/environmental variables. Since it is already pre-processed, we can work directly with it, out of the box.


## Project description
### Structure
```
├── report/*                # A report summarising results
├── scripts/                # Core scripts to run
├── R/                      # Functions
├── data/
│   ├── config              # Configuration files
│   ├── raw_data            # Original datasets (hidden)
│   └── preprocessed_data   # Preprocessed datasets
├── outputs/
│   ├── figures             # Visualizations and plots
│   └── results             # Data outputs
├── renv/                   # Information about R environment (packages)
├── DESCRIPTION             # Standard DESCRIPTION file for R packages
└── README.md               # This file
```

### Progress tracker
- [X] Anonymize agricultural plots locations
- [X] Load STOC dataset
- [X] Pre-processing
- [X] Training of HMSC models
- [X] Implementation of custom cost function
- [X] OED using HMSC on STOC
- [X] Add an automatic report
- [X] Does adding samples has an effect ?
- [ ] Explain results using simulated data ?
- [ ] Add other models (GJAM, RF, DNN,...)


## Getting Started
### Reading results
Open the report in your browser (available [here](https://ll-mnhn.github.io/pseudo-STOC-OED/report/report.html)). Or, alternatively, read the PDF equivalent (available [here](https://github.com/LL-mnhn/pseudo-STOC-OED/blob/main/report/report.pdf)).

### How to use scripts
1\. Clone this repository on your machine
```bash
cd /your/local/folder
git clone https://github.com/LL-mnhn/pseudo-STOC-OED.git
```

2\. Install dependencies

Open `pseudo-STOC-OED` as a new session in [Rstudio](https://docs.posit.co/ide/user/) or [Positron](https://positron.posit.co/welcome.html). Use R 4.6.1 (version used during development) 

Install `renv` if not already installed on your machine. Then run:
```R
install.packages("renv")
renv::restore()
```

3\. Once your environment is ready, local scripts can be run. E.g.

```R
source("scripts/4-summary.R")
```

### Usage Notes
Important results and figures are already saved in the `outputs` folder.

When running on an external machine: 
- `0-verify_datasets.R` can not be run (raw dataset are only accessible by authors)
- `1-pre_processing.R` can only be run partially (it processes raw datasets, which are unavailable; but it can show and save basic figures).
- Other files in `scripts` can be run without restrictions.

## Contact
For inquiries, please contact Loïc Lehnhoff (UMR CESCO - MNHN) at <loic.lehnhoff@mnhn.fr> 

*This work is supervised by Nicolas Parisey (UMR IGEPP - INRAE) and Karine Princé (UMR CESCO - MNHN)*
