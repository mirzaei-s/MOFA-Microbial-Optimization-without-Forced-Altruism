# MOFA: Microbial Optimization without Forced Altruism 
A nonlinear programming method for analyzing microbial community-level growth rate without enforced altruism.
This repository contains a Python implementation of MOFA (Microbial Optimization without Forced Altruism) applied to microbial communities.

## Overview

The goal of this project is to model and optimize the growth of microbial communities by merging individual species' metabolic networks and applying nonlinear optimization to predict community behavior. MOFA allows the study of microbial interactions without imposing artificial constraints on altruistic behavior.

#

### DV&MM Folder

This folder contains the MOFA implementation for **Desulfovibrio vulgaris** and **Methanococcus maripaludis**, as used in the article. 

#

### Gut Microbiome Species Folder
- For gut microbiome species, the Python script `MOFA-Gut-Microbiome.py` and the species’ metabolic network files are located in the `model` folder. Make sure to update the file paths in the code to match your local or working directory, and specify the correct model names and biomass IDs for each species. (We have already added this in the code—just uncomment the relevant lines and select the diet for each species.)

- All results of MOFA for this community are saved in `OurResult.xlsx`.

#
**Toy Model:**  
   Use the metabolic network files `organism_Q_sbml.xml` and `organism_P_sbml.xml` with the notebook `MOFA_toymodel.ipynb`.

#

### MOFA Folder
- Contains `MOFA.py`, a version of the MOFA algorithm that users can apply to their own models. A `ReadMe.pdf` file is also included, providing step-by-step instructions for users who want to apply MOFA to their own models.

#

### MOFA(matlab)
- This folder contains the MOFA version in which the model-merging function is implemented in MATLAB. To use this implementation, you need to upload the `.S` file, which is an Excel output representing the merged metabolic network of the community model returned from MATLAB to the `.ipynb` notebook.
