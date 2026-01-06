# MOFA: Microbial Optimization without Forced Altruism (DV&MM)

## Notebook: DV&MM_MOFA.ipynb

This Jupyter notebook provides a full workflow for:

1. **Reading metabolic networks**: Import models in `.smbl` or `.mat` formats.
2. **Merging species networks**: The function `add_model_to_community()` combines individual species' metabolic networks to create a comprehensive community model.
3. **Setting the diet**: Apply specific environmental conditions or nutrient uptake rates to the community.
4. **Running MOFA optimization**: Solve the nonlinear optimization problem to determine community growth rates and flux distributions, with constraints such as lactate (Lac) and hydrogen (H₂) set to experimental values.

## Functions

- `add_model_to_community(model, communitymodel, species_id)`: Merges a species model into the community model.

## Requirements

- Python 3.x
- COBRApy
- NumPy, SciPy, pandas
- Jupyter Notebook
- **IPOPT solver** (required for nonlinear optimization)
- Optional: MATLAB `.mat` files reader (`scipy.io`) if using `.mat` metabolic models
