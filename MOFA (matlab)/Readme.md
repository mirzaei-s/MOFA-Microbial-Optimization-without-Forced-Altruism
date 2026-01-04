# MOFA: Microbial Optimization without Forced Altruism 
A nonlinear programming method for analyzing microbial community-level growth rate without enforced altruism.

### Step 1: Construct the community metabolic network

First, use the **Community_Growth_Rate** MATLAB code to construct a merged metabolic network for the community.  
Specify the names of your species’ metabolic network models, for example:

```matlab
models_name = {'organism_Q_sbml', 'organism_P_sbml'};
```
The S file is the Excel output representing the merged metabolic network of the community model.

The Create_Community_Model.m function is used to merge individual species into a community model.

### Step 2: Run the Python Notebook

Run the provided `.ipynb` file and use the **S file** (Excel output) as input to compute the community growth rate.  

**NECom Implementation**

1. Use the `Create_Necom_Model.m` function to generate the stoichiometric matrix for each species, ensuring compatibility with the NECom framework.
