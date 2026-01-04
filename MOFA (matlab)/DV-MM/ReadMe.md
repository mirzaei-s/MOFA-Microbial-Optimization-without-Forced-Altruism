### Step 1: Construct the community metabolic network

First, use the **Community_Growth_Rate** MATLAB code to construct a merged metabolic network for the community.  
Specify the names of your species’ metabolic network models, for example:

```matlab
models_name = {'MM','DV'};
```
The S file is the Excel output representing the merged metabolic network of the community model.

**Desulfovibrio vulgaris & Methanococcus maripaludis:**  
- For *Desulfovibrio vulgaris* and *Methanococcus maripaludis* species, see the `DV-MM` folder and run `MOFA_for_DV&MM.ipynb`. Upload the file `S_DV_MM.xlsx` to the path `/content/S_DV_MM.xlsx` in the Jupyter notebook environment before execution.
