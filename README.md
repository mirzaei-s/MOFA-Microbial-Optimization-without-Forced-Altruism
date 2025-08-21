# Microbial Community Growth Rate

MOFA: Microbial Optimization without Forced Altruism

First, use Community_Growth_Rate to construct a merged metabolic network for the community. The S file, which is the Excel output of the community model, is then used in Python to compute the community growth rate.

For NECom implementation, use the CreateNECom function to generate the stoichiometric matrix for each species, ensuring compatibility with the NECom framework.

For each pair of species in our study, we placed the corresponding S (stoichiometric) Excel file in their respective folders, updated the file path in the read_excel function, and then ran the method’s function.
