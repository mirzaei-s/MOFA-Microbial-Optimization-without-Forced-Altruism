# MOFA: Microbial Optimization without Forced Altruism 
A nonlinear programming method for analyzing microbial community-level growth rate without enforced altruism.

### Step 1: Construct the community metabolic network

First, use the **Community_Growth_Rate** MATLAB code to construct a merged metabolic network for the community.  
Specify the names of your species’ metabolic network models, for example:

```matlab
models_name = {'organism_Q_sbml', 'organism_P_sbml'};
```
The S file is the Excel output representing the merged metabolic network of the community model.

### Step 2: Run the Python Notebook

Run the provided `.ipynb` file and use the **S file** (Excel output) as input to compute the community growth rate.  

**Toy Model:**
- For the toy model, the notebook is located in the `toymodel` folder: run `NoFA_Opt_&_NECom_(Toymodel).ipynb` with the file `S_toymodel.xlsx`.
  
**Desulfovibrio vulgaris & Methanococcus maripaludis:**  
- For *Desulfovibrio vulgaris* and *Methanococcus maripaludis* species, see the `DV-MM` folder and run `MOFA_for_DV&MM.ipynb` with the file `S_DV_MM.xlsx`.

**Gut microbiome species:**
- For gut microbiome species, the Python code `MOFA-Gut-Microbiome.py` uses `S_human.xlsx` located in the zipped folder. All results for this community can be found in `OurResult.xlsx`.


  
**NECom Implementation**

1. Use the `Create_Necom_Model.m` function to generate the stoichiometric matrix for each species, ensuring compatibility with the NECom framework.

2. **Toy Model:**  
   Use the Excel files `S_organism_P_sbml.xls` and `S_organism_Q_sbml.xls` with the notebook `NoFA_Opt_&_NECom_(Toymodel).ipynb`.

3. **Desulfovibrio vulgaris & Methanococcus maripaludis:**  
   Use the Excel files `S_DV.xls` and `S_MM.xls` with the notebook `NECom_DV&MM.ipynb`.
