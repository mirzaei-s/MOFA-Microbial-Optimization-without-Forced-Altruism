# Microbial Community Growth Rate

MOFA: Microbial Optimization without Forced Altruism

### Step 1: Construct the community metabolic network

First, use the **Community_Growth_Rate** MATLAB code to construct a merged metabolic network for the community.  
Specify the names of your species’ metabolic network models, for example:

```matlab
models_name = {'organism_Q_sbml', 'organism_P_sbml'};
```
The S file is the Excel output representing the merged metabolic network of the community model.

### Step 2: Run the Python Notebook

Run the provided `.ipynb` file and use the **S file** (Excel output) as input to compute the community growth rate.  

- For the toy model, the notebook is located in `toymodel` file:  `NoFA_Opt_&_NECom_(Toymodel).ipynb`
- For *Desulfovibrio vulgaris* and *Methanococcus maripaludis* species, see the `DV-MM` file.
- For gut microbiome species, the Python code is in the `gut microbiome` folder.
- 
For NECom implementation, use the CreateNECom function to generate the stoichiometric matrix for each species, ensuring compatibility with the NECom framework.

For each pair of species in our study, we placed the corresponding S (stoichiometric) Excel file in their respective folders, updated the file path in the read_excel function, and then ran the method’s function.
