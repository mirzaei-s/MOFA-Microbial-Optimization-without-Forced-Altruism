# MOFA: Microbial Optimization without Forced Altruism 
A nonlinear programming method for analyzing microbial community-level growth rate without enforced altruism.
This repository contains a Python implementation of MOFA (Microbial Optimization without Forced Altruism) applied to microbial communities.

## Overview

The goal of this project is to model and optimize the growth of microbial communities by merging individual species' metabolic networks and applying nonlinear optimization to predict community behavior. MOFA allows the study of microbial interactions without imposing artificial constraints on altruistic behavior.

### DV&MM Folder

This folder contains the MOFA implementation for **Desulfovibrio vulgaris** and **Methanococcus maripaludis**, as used in the article. 

**Gut microbiome species:**
- For gut microbiome species, the Python script `MOFA-Gut-Microbiome.py` and the species’ metabolic network files are located in the `model` folder. Make sure to update the file paths in the code to match your local or working directory, and specify the correct model names and biomass IDs for each species. (We have already added this in the code—just uncomment the relevant lines and select the diet for each species.)

- All results of MOFA for this community are saved in `OurResult.xlsx`.



