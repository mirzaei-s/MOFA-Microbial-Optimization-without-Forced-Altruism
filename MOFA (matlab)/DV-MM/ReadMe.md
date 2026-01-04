### Step 1: Construct the community metabolic network

First, use the **Community_Growth_Rate** MATLAB code to construct a merged metabolic network for the community.  
Specify the names of your species’ metabolic network models, for example:

```matlab
models_name = {'MM','DV'};
```
The S file is the Excel output representing the merged metabolic network of the community model.

**Desulfovibrio vulgaris & Methanococcus maripaludis:**  
- For *Desulfovibrio vulgaris* and *Methanococcus maripaludis* species, see the `DV-MM` folder and run `MOFA_for_DV&MM.ipynb`. Upload the file `S_DV_MM.xlsx` to the path `/content/S_DV_MM.xlsx` in the Jupyter notebook environment before execution.


**NECom Implementation**

1. Use the `Create_Necom_Model.m` function to generate the stoichiometric matrix for each species, ensuring compatibility with the NECom framework.

2. **Toy Model:**  
   Use the Excel files `S_organism_P_sbml.xls` and `S_organism_Q_sbml.xls` with the notebook `NoFA_Opt_&_NECom_(Toymodel).ipynb`.

3. **Desulfovibrio vulgaris & Methanococcus maripaludis:**  
   Use the Excel files `S_DV.xls` and `S_MM.xls` with the notebook `NECom_DV&MM.ipynb`.
