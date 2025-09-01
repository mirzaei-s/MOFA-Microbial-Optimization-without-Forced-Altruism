# MOFA: Microbial Optimization without Forced Altruism 
 Microbial Community Growth Rate


### Step 1: Construct the community metabolic network

First, use the **Community_Growth_Rate** MATLAB code to construct a merged metabolic network for the community.  
Specify the names of your species’ metabolic network models, for example:

```matlab
models_name = {'organism_Q_sbml', 'organism_P_sbml'};
```
The S file is the Excel output representing the merged metabolic network of the community model.

### Step 2: Run the Python Notebook

Run the provided `.ipynb` file and use the **S file** (Excel output) as input to compute the community growth rate.  

- For the toy model, the notebook is located in the `toymodel` folder:  
  `NoFA_Opt_&_NECom_(Toymodel).ipynb` with the file `S_toymodel.xlsx`.
- For *Desulfovibrio vulgaris* and *Methanococcus maripaludis* species, see the `DV-MM` folder and run `MOFA_for_DV&MM.ipynb` with the file `S_DV_MM.xlsx`.
- For gut microbiome species, the Python code `MOFA-Gut-Microbiome.py` uses `S_human.xlsx` located in the zipped folder.


  
For NECom implementation, use the `Create_Necom_Model.m` function to generate the stoichiometric matrix for each species, ensuring compatibility with the NECom framework.
 use the excel file to run NECom_DV&MM.ipynb


For each pair of species in our study, we placed the corresponding S (stoichiometric) Excel file in their respective folders, updated the file path in the read_excel function, and then ran the method’s function.
