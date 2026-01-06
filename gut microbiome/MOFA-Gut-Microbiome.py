import os
from pathlib import Path
from cobra.io import load_matlab_model,read_sbml_model
from cobra import Model, Reaction, Metabolite
from cobra.util.array import create_stoichiometric_matrix
from cobra.util.solver import linear_reaction_coefficients
import numpy as np
import pandas as pd
import re
import copy
import pyomo.environ as pyo
from pyomo.environ import *
import time

def merge_models_block_diagonal(model1, model2, prefix2="m2_"):
    """
    Merge two COBRA models with NO shared metabolites or reactions.
    Reactions and metabolites of model2 are appended after model1.
    """


    m1 = model1.copy()
    m2 = model2.copy()

    # ---------- rename metabolites of model2 ----------
    for met in m2.metabolites:
        met.id =  met.id + prefix2

    # ---------- rename reactions of model2 ----------
    for rxn in m2.reactions:
        rxn.id =  rxn.id + prefix2

    # ---------- add metabolites ----------
    m1.add_metabolites(m2.metabolites)

    # ---------- add reactions ----------
    m1.add_reactions(m2.reactions)
    return m1


def add_model_to_community(model, communitymodel, species_id):

  EX_1 = [i for i, rxn in enumerate(model.reactions) if 'EX_' in rxn.id]
  EX_2 = [i for i, rxn in enumerate(communitymodel.reactions) if 'EX_' in rxn.id]

  i1, _ = np.where(create_stoichiometric_matrix(model)[:, EX_1] != 0)
# community model
  i2, _ = np.where(create_stoichiometric_matrix(communitymodel)[:, EX_2] != 0)


  mets_1 = {
    met.id
    for rxn in model.reactions if 'EX_' in rxn.id
    for met in rxn.metabolites
  }
  mets_2 = {
    met.id.split('_species')[0]
    for rxn in communitymodel.reactions if 'EX_' in rxn.id
    for met in rxn.metabolites
   }

  shared_mets = list(set(mets_1) & set(mets_2))
  for i_loop in range(len(shared_mets)):
     met1_index = next(
        i for i, met in enumerate(model.metabolites)
        if met.id == shared_mets[i_loop]
     )

    # reactions involving this metabolite
     temp = np.where(create_stoichiometric_matrix(model)[met1_index, :] != 0)[0]
    # exchange reactions involving this metabolite
     rxn1_index = [
        j for j in temp
        if 'EX_' in model.reactions[j].id
     ]

     if model.reactions[rxn1_index[0]].lower_bound < 0 and model.reactions[rxn1_index[0]].upper_bound > 0:

      # Check if the reaction exports the metabolite (negative stoichiometric coefficient)
      if model.reactions[rxn1_index[0]].get_coefficient(model.metabolites[met1_index])<0:

          new_rxn = Reaction('EX_uptake_' + model.reactions[rxn1_index[0]].id.split('EX_')[1])
          new_rxn.lower_bound = 0
          new_rxn.upper_bound = -1 * model.reactions[rxn1_index[0]].lower_bound


          # Set the stoichiometry for the metabolite
          rxn = model.reactions[rxn1_index[0]]  # get the reaction object
          met = model.metabolites[met1_index]   # get the metabolite object
          stoich_coeff = rxn.get_coefficient(met)  # negative if metabolite is consumed

          new_rxn.add_metabolites({met: -1 * stoich_coeff})

          # Add new reaction to the model
          model.add_reactions([new_rxn])
          new_rxn.objective_coefficient = 0
          # Rename the original reaction to EX_export_...
          model.reactions[rxn1_index[0]].id = 'EX_export_' + model.reactions[rxn1_index[0]].id.split('EX_')[1]
          model.reactions[rxn1_index[0]].lower_bound = 0


      else:

          new_rxn = Reaction('EX_export_' + model.reactions[rxn1_index[0]].id.split('EX_')[1])
          new_rxn.lower_bound = 0
          new_rxn.upper_bound = -1 * model.reactions[rxn1_index[0]].lower_bound


          # Set the stoichiometry for the metabolite
          rxn = model.reactions[rxn1_index[0]]  # get the reaction object
          met = model.metabolites[met1_index]   # get the metabolite object
          stoich_coeff = rxn.get_coefficient(met)  # negative if metabolite is consumed

          new_rxn.add_metabolites({met: -1 * stoich_coeff})
          # Add new reaction to the model
          model.add_reactions([new_rxn])
          # Rename the original reaction to EX_export_...
          new_rxn.objective_coefficient = 0


     elif  model.reactions[rxn1_index[0]].get_coefficient(model.metabolites[met1_index])<0:

      model.reactions[rxn1_index[0]].id = 'EX_export_' + model.reactions[rxn1_index[0]].id.split('EX_')[1]
     else:
      model.reactions[rxn1_index[0]].id = 'EX_uptake_' + model.reactions[rxn1_index[0]].id.split('EX_')[1]

 #--------------------------------------------------------------------------------------------------
#                                             second model
 #---------------------------------------------------------------------------------------------------


     mets_base = [m.id.split('_species')[0] for m in communitymodel.metabolites]

     # find metabolite index
     met2_index_arr = np.where(np.array(mets_base) == shared_mets[i_loop])[0]
     met2_index = met2_index_arr[0]

     # find reactions involving this metabolite
     S_comm = create_stoichiometric_matrix(communitymodel)
     temp = np.where(S_comm[met2_index, :] != 0)[0]

    # exchange reactions involving this metabolite
     rxn2_index = [
        j for j in temp
        if 'EX_' in communitymodel.reactions[j].id
     ]
     rxn_ids = np.array([rxn.id for rxn in communitymodel.reactions])
     if not np.any(np.char.find(rxn_ids[rxn2_index], '_[Env]') != -1):
      if communitymodel.reactions[rxn2_index[0]].lower_bound < 0 and communitymodel.reactions[rxn2_index[0]].upper_bound > 0:

      # Check if the reaction exports the metabolite (negative stoichiometric coefficient)
       if communitymodel.reactions[rxn2_index[0]].get_coefficient(communitymodel.metabolites[met2_index])<0:

          rxn2_origin=communitymodel.reactions[rxn2_index[0]].id.split('EX_')[1]
          oldlb=communitymodel.reactions[rxn2_index[0]].lower_bound
           # Rename the original reaction to EX_export_...
          communitymodel.reactions[rxn2_index[0]].id = 'EX_export_' +rxn2_origin
          communitymodel.reactions[rxn2_index[0]].lower_bound = 0

          new_rxn = Reaction('EX_uptake_' + rxn2_origin)
          new_rxn.lower_bound = 0
          new_rxn.upper_bound = -1 * oldlb

        # Set the stoichiometry for the metabolite
          rxn = communitymodel.reactions[rxn2_index[0]]  # get the reaction object
          met = communitymodel.metabolites[met2_index]   # get the metabolite object
          stoich_coeff = rxn.get_coefficient(met)  # negative if metabolite is consumed
          new_rxn.add_metabolites({met: -1 * stoich_coeff})
          # Add new reaction to the model
          communitymodel.add_reactions([new_rxn])
          new_rxn.objective_coefficient = 0


       else:
          rxn2_origin=communitymodel.reactions[rxn2_index[0]].id.split('EX_')[1]
           # Rename the original reaction to EX_export_...
          communitymodel.reactions[rxn2_index[0]].id = 'EX_uptake_' + rxn2_origin
          oldlb=communitymodel.reactions[rxn2_index[0]].lower_bound;
          communitymodel.reactions[rxn2_index[0]].lower_bound = 0

          new_rxn = Reaction('EX_export_' + rxn2_origin)
          new_rxn.lower_bound = 0
          new_rxn.upper_bound = -1 * oldlb

          # Set the stoichiometry for the metabolite
          rxn = communitymodel.reactions[rxn2_index[0]]  # get the reaction object
          met = communitymodel.metabolites[met2_index]   # get the metabolite object
          stoich_coeff = rxn.get_coefficient(met)  # negative if metabolite is consumed

          new_rxn.add_metabolites({met: -1 * stoich_coeff})
          # Add new reaction to the model
          communitymodel.add_reactions([new_rxn])
          new_rxn.objective_coefficient = 0


      elif  communitymodel.reactions[rxn2_index[0]].get_coefficient(communitymodel.metabolites[met2_index])<0:
       communitymodel.reactions[rxn2_index[0]].id = 'EX_export_' + communitymodel.reactions[rxn2_index[0]].id.split('EX_')[1]
      else:
       communitymodel.reactions[rxn2_index[0]].id = 'EX_uptake_' + communitymodel.reactions[rxn2_index[0]].id.split('EX_')[1]


################################################################################
##                      Merge models
################################################################################


  communitymodel = merge_models_block_diagonal(communitymodel, model, prefix2="_species2")

  for i_loop in range(len(shared_mets)):
    mets_base = [m.id.split('_species')[0] for m in communitymodel.metabolites]
    met_index_arr = np.where(np.array(mets_base) == shared_mets[i_loop])[0]
    met_index = met_index_arr
    temp = np.where(np.any(create_stoichiometric_matrix(communitymodel)[met_index, :] != 0, axis=0))[0]
    rxn_index = [
    j for j in temp
    if 'EX_' in communitymodel.reactions[j].id
    ]

    if not any('_[Env]' in communitymodel.reactions[i].id for i in rxn_index):

      env_met_id = f"{shared_mets[i_loop]}"
      if env_met_id not in [m.id for m in communitymodel.metabolites]:
        env_met = Metabolite(
        id=env_met_id,
        name=env_met_id,
        compartment="e"   # یا هر compartment دلخواه
        )
        communitymodel.add_metabolites([env_met])

        rxn_id = f"EX_export_{shared_mets[i_loop]}_[Env]"
        new_rxn = Reaction(rxn_id)
        new_rxn.lower_bound = 0
        new_rxn.upper_bound = 1000
        new_rxn.add_metabolites({ env_met: -1})
        communitymodel.add_reactions([new_rxn])

        rxn_id = f"EX_uptake_{shared_mets[i_loop]}_[Env]"
        new_rxn = Reaction(rxn_id)
        new_rxn.lower_bound = 0
        new_rxn.upper_bound = 1000
        new_rxn.add_metabolites({ env_met: 1})
        communitymodel.add_reactions([new_rxn])

        for idx in rxn_index:
          rxn = communitymodel.reactions[idx]
          rxn.id = 'EXCom_' + rxn.id.split('EX_')[-1]

          jj = np.where(create_stoichiometric_matrix(communitymodel)[:, idx] != 0)[0]

          met = communitymodel.metabolites[jj[0]]   # get the metabolite object

          stoich_coeff = rxn.get_coefficient(met)  # negative if metabolite is consumed
          rxn.add_metabolites({env_met:  -1 * stoich_coeff})

          if '_species' not in met.id:
            suffix = rxn.id.split('_species', 1)[1]
            met.id = met.id + '_species' + suffix



    else:
      rxn = [j for j in rxn_index if 'species' in communitymodel.reactions[j].id]
      for j in rxn:
       old_id = communitymodel.reactions[j].id
       if 'EX_' in old_id:
         communitymodel.reactions[j].id = 'EXCom_' + old_id.split('EX_')[1]

      ii= np.where(create_stoichiometric_matrix(communitymodel)[:,  rxn[0] ] != 0)[0]

      rxn_ids = np.array([rxn.id for rxn in communitymodel.reactions])
      temp = np.where(np.char.find(rxn_ids[rxn_index], "_[Env]") != -1)[0]

      for j_idx in rxn:
       rxn_obj = communitymodel.reactions[j_idx]


       main_met = communitymodel.metabolites[ii]
       env_met = communitymodel.metabolites[jj[0]]

       coeff = rxn_obj.get_coefficient(main_met)
       rxn_obj.add_metabolites({env_met: -coeff})

  pattern = re.compile(r'^EX_.*?(?<!_\[Env\])$')

  Ex_index = [i for i, rxn in enumerate(communitymodel.reactions)
            if pattern.search(rxn.id)]

  for i in range(len(Ex_index)):
    idx = Ex_index[i]
    rxn = communitymodel.reactions[idx]
    if rxn.lower_bound < 0 and rxn.upper_bound > 0:
      S = create_stoichiometric_matrix(communitymodel, array_type='dense')
      met_index = np.where(S[:, Ex_index[i]] != 0)[0]

      met, coeff = next(iter(rxn.metabolites.items()))
      name_rxn = rxn.id.split('EX_', 1)[1]


      if coeff < 0:

          # rename export reaction
          rxn.id = 'EXCom_export_' + name_rxn
          old_lb=rxn.lower_bound
          rxn.lower_bound = 0

          uptake = Reaction(id='EXCom_uptake_' + name_rxn)
          uptake.add_metabolites({met: coeff*-1})
          uptake.lower_bound = 0
          uptake.upper_bound = old_lb*-1

          communitymodel.add_reactions([uptake])

      else:

          # rename to uptake
          rxn.id = 'EXCom_uptake_' + name_rxn

          # update bounds
          old_lb = rxn.lower_bound
          rxn.lower_bound = 0


          # create export reaction
          export = Reaction(id='EXCom_export_' + name_rxn)
          export.add_metabolites({met: -1*coeff})
          export.lower_bound = 0
          export.upper_bound = -1*old_lb

          communitymodel.add_reactions([export])
    else:
      rxn = communitymodel.reactions[Ex_index[i]]
      met, coeff = next(iter(rxn.metabolites.items()))
      name_rxn = rxn.id.split('EX_', 1)[1]

      if coeff < 0:
         rxn.id = 'EXCom_export_' + name_rxn
      else:
          rxn.id = 'EXCom_uptake_' + name_rxn
  return communitymodel


SAVEDIR = Path("D:\\Gut microbiome models")

models_name = ['Bacteroides_caccae_ATCC_43185.mat', 'Bifidobacterium_longum_infantis_ATCC_15697.mat']

model_path = os.path.join(SAVEDIR, f"{models_name[0]}")
file_path = Path(model_path)
suffix = file_path.suffix.lower()

if suffix == ".mat":
        print("Reading MAT file...")
        communitymodel = load_matlab_model(model_path)

        model_path = os.path.join(SAVEDIR, f"{models_name[1]}")
        model = load_matlab_model(model_path)

    # -------- XML / SBML file (COBRA) --------
elif suffix == ".xml":
        print("Reading SBML model using COBRA Toolbox...")

        communitymodel = read_sbml_model(model_path)
        model_path = os.path.join(SAVEDIR, f"{models_name[1]}")
        model = read_sbml_model(model_path)


else:
        raise ValueError("Only .mat and .xml (SBML) files are supported.")


for rxn in communitymodel.reactions:
    if abs(rxn.lower_bound) < 1e-12:
        rxn.lower_bound = 0.0
    if abs(rxn.upper_bound) < 1e-12:
        rxn.upper_bound = 0.0

for rxn in model.reactions:
    if abs(rxn.lower_bound) < 1e-12:
        rxn.lower_bound = 0.0
    if abs(rxn.upper_bound) < 1e-12:
        rxn.upper_bound = 0.0
# List of model names

for rxn in communitymodel.reactions:
    rxn.id += "_species1"

for met in communitymodel.metabolites:
    met.id += "_species1"


i=0
communitymodel = add_model_to_community(
        model,
        communitymodel,
        species_id=i+2
    )

biomass_id=['biomass536_species1','biomass525_species2'] # Biomass human 2
#biomass_id=['biomass223_species1','biomass525_species2'] # Biomass human7
#biomass_id=['biomass525_species1','biomass525_species2'] # Biomass human5
#biomass_id=['biomass223_species1','biomass525_species2'] # Biomass human6
#biomass_id=['biomass345_species1','biomass345_species2'] #Biomass human3

diets_name=["YCAG.xlsx","YCGD.xlsx","mMCB.xlsx","mMCB.xlsx","YCGMS.xlsx"]
diet_path = os.path.join(SAVEDIR, f"{diets_name[0]}")
diet = pd.read_excel(diet_path)

for i in range(len(diet)):
    met_name = f"{diet.iloc[i, 0]}[e]"
    ub_value = diet.iloc[i, 1]

    # Try shared extracellular metabolite
    met = communitymodel.metabolites.get_by_id(met_name) \
          if met_name in communitymodel.metabolites else None

    if met is not None:
        for rxn in met.reactions:
            if 'EX_uptake' in rxn.id:
                rxn.upper_bound = ub_value
            elif 'EX_export' in rxn.id or 'EXCom_export' in rxn.id:
                rxn.upper_bound = 0

    else:
        # Species-specific metabolites
        for species in ['species1', 'species2']:
            met_name_sp = f"{diet.iloc[i, 0]}[e]_{species}"
            if met_name_sp in communitymodel.metabolites:
                met_sp = communitymodel.metabolites.get_by_id(met_name_sp)

                for rxn in met_sp.reactions:
                    if 'EXCom_uptake' in rxn.id:
                        rxn.upper_bound = ub_value
                    elif 'EXCom_export' in rxn.id:
                        rxn.upper_bound = 0

m,n=create_stoichiometric_matrix(communitymodel).shape

rxns = [rxn.id for rxn in communitymodel.reactions]
mets = [met.id for met in communitymodel.metabolites]

S = create_stoichiometric_matrix(communitymodel, array_type='dense')
Sij = pd.DataFrame(
    S,
    index=mets,
    columns=rxns
)


lb = [rxn.lower_bound for rxn in communitymodel.reactions]
ub = [rxn.upper_bound for rxn in communitymodel.reactions]


substring='EXCom_uptake'
indices = [i for i, s in enumerate(rxns) if substring in s]
uptake_rxns_indices = [rxns[i] for i in indices]

uptake_species1 = [n for n in uptake_rxns_indices if 'species1' in n]

uptake_species2=[n for n in uptake_rxns_indices if 'species2' in n]

substring='EXCom_export'
indices = [i for i, s in enumerate(rxns) if substring in s]
export_rxns_indices = [rxns[i] for i in indices]

big_M=10000000000000

Exc_indice = [i for i, s in enumerate(rxns) if 'EXCom_' in s]
Exc = [rxns[i] for i in Exc_indice]

F1=1
F2=1


model = ConcreteModel()

model.N = Set(initialize=rxns)

model.M = Set(initialize=mets)


model.biomass = Set(initialize=biomass_id, within=model.N)

model.uptake_rxns_species1=Set(initialize=uptake_species1,within=model.N)

model.uptake_rxns_species2=Set(initialize=uptake_species2,within=model.N)

model.export_rxns=Set(initialize=export_rxns_indices,within=model.N)

    ##############################################################################################
    #              Create Parameters
    ##############################################################################################

model.lb = pyo.Param(model.N, initialize={rxn: lb[j]  for j, rxn in enumerate(rxns)})
model.ub = pyo.Param(model.N, initialize={rxn: ub[j]  for j, rxn in enumerate(rxns)})

    ##############################################################################################
    #              Create Variables
    ##############################################################################################

index_to_list = {k: i for i, k in enumerate(model.N)}

def variable_bounds_from_list(model, index):
  idx = index_to_list[index]
 # if 'EX_' not in index:
  #  lower_bound = lb[idx]*0.5
  #  upper_bound = ub[idx]*0.5
  #else:
  lower_bound = lb[idx]
  upper_bound = ub[idx]

  return (lower_bound, upper_bound)

model.v = pyo.Var(model.N, bounds=variable_bounds_from_list)

model.u_1=pyo.Var(model.N)

model.lamda_1=pyo.Var(model.M)

model.eta_UB_1=pyo.Var(model.N,domain=pyo.NonNegativeReals)

model.eta_LB_1=pyo.Var(model.N,domain=pyo.NonNegativeReals)


model.u_2=pyo.Var(model.N)

model.lamda_2=pyo.Var(model.M)

model.eta_UB_2=pyo.Var(model.N,domain=pyo.NonNegativeReals)

model.eta_LB_2=pyo.Var(model.N,domain=pyo.NonNegativeReals)

##############################################################################################
    #                   The objective function
##############################################################################################

model.obj=pyo.Objective(expr=sum(model.v[i] for i in model.biomass), sense=pyo.maximize)

#S = {(r, c): Sij.at[r, c] * ( F1 if 'species1' in c else F2 ) if c in Exc else Sij.at[r, c] for r in Sij.index for c in Sij.columns}
S = {(r, c):Sij.at[r,c] for r in  Sij.index for c in Sij.columns}

def massbalance_rule(model,m):
     return sum(S[m,n]*model.v[n] for n in model.N)==0
model.massbalance=pyo.Constraint(model.M,rule=massbalance_rule)

model.Biomass_species1=pyo.Constraint(expr=model.v[biomass_id[0]]>= model.u_1[biomass_id[0]] )

model.Biomass_species2=pyo.Constraint(expr=model.v[biomass_id[1]]>=model.u_2[biomass_id[1]] )
big_M=100000000

def Export_couple_Biomass_rule(model, export):
     if 'species1' in export:
        return model.v[export] <= big_M* model.v[biomass_id[0]]
     return model.v[export] <= big_M* model.v[biomass_id[1]]

model.Export_couple_biomass=pyo.Constraint(model.export_rxns,rule=Export_couple_Biomass_rule)


def Zero_const_species_rule(model, r):
    if r in uptake_rxns_indices:
        rr = r.replace('uptake', 'export')
        return model.v[r] * model.v[rr] == 0.000000000000000
    else:
        return pyo.Constraint.Skip

#model.zero_species2 = pyo.Constraint(model.N, rule=Zero_const_species_rule)
model.zero_species2 =pyo.Constraint(expr=sum(model.v[r] * model.v[r.replace('uptake', 'export')] for r in uptake_rxns_indices)==0)


##############################################################################################
    #                   (Inner-problem species1) constraint F.A for species 1
###############################################################################################
def massbalance_species1_rule(model,m):
     return sum(S[m,n]*model.u_1[n] for n in model.N)==0

model.massbalance_species1=pyo.Constraint(model.M,rule=massbalance_species1_rule)  #Su^1=0

def UB_const_species1_rule(model, r):
        if r not in model.uptake_rxns_species1:
            return model.u_1[r] <= model.ub[r]
        return  model.u_1[r] <= model.v[r]
model.ub_species1= pyo.Constraint(model.N, rule=UB_const_species1_rule)

def LB_const_species1_rule(model, r):
     return -1*model.u_1[r] <= -1*model.lb[r]
model.lb_species1= pyo.Constraint(model.N, rule=LB_const_species1_rule)


i_biomass=0
def dual_const_species1_rule(model, n):
    if n not in  biomass_id[i_biomass]:
        return sum(S[m,n]*model.lamda_1[m] for m in model.M)+ model.eta_UB_1[n]-model.eta_LB_1[n]==0

    return sum(S[m,biomass_id[i_biomass]]*model.lamda_1[m] for m in model.M)+model.eta_UB_1[biomass_id[i_biomass]]-model.eta_LB_1[biomass_id[i_biomass]]==1

model.dual_species1= pyo.Constraint(model.N, rule=dual_const_species1_rule)

model.dual_eq_primal_species1=pyo.Constraint(expr= (sum(model.v[n]*model.eta_UB_1[n] for n in model.uptake_rxns_species1)+ \
                  sum(-1*model.lb[n]*model.eta_LB_1[n]  for n in model.uptake_rxns_species1)+ \
                  sum(model.ub[n]*model.eta_UB_1[n] for n in model.N if n not in model.uptake_rxns_species1)+\
                  sum(-1*model.lb[n]*model.eta_LB_1[n] for n in model.N if n not in model.uptake_rxns_species1))==model.u_1[biomass_id[i_biomass]])

 ##############################################################################################
    #                   (Inner-problem species2) constraint F.A for species 2
###############################################################################################

def massbalance_species2_rule(model,m):
    return sum(S[m,n]*model.u_2[n] for n in model.N)==0
model.massbalance_species2=pyo.Constraint(model.M,rule=massbalance_species2_rule)

def UB_const_species2_rule(model, r):
    if r not in model.uptake_rxns_species2:
        return model.u_2[r] <= model.ub[r]
    return model.u_2[r] <= model.v[r]
model.ub_species2= pyo.Constraint(model.N, rule=UB_const_species2_rule)

def LB_const_species2_rule(model, r):
      return -1*model.u_2[r] <= -1*model.lb[r]
model.lb_species2= pyo.Constraint(model.N, rule=LB_const_species2_rule)

i_biomass=1
def dual_const_species2_rule(model, n):
    if n not in  biomass_id[i_biomass]:
       return sum(S[m,n]*model.lamda_2[m] for m in model.M)+ model.eta_UB_2[n]-model.eta_LB_2[n]==0
    return sum(S[m,biomass_id[i_biomass]]*model.lamda_2[m] for m in model.M)+ \
            model.eta_UB_2[biomass_id[i_biomass]] -model.eta_LB_2[biomass_id[i_biomass]]==1

model.dual_species2= pyo.Constraint(model.N, rule=dual_const_species2_rule)

model.dual_eq_primal_species2=pyo.Constraint(expr= (sum(model.v[n]*model.eta_UB_2[n] for n in model.uptake_rxns_species2)+ \
                  sum(-1*model.lb[n]*model.eta_LB_2[n] for n in model.uptake_rxns_species2)+ \
                  sum(model.ub[n]*model.eta_UB_2[n] for n in model.N if n not in model.uptake_rxns_species2)+\
                  sum(-1*model.lb[n]*model.eta_LB_2[n] for n in model.N if n not in model.uptake_rxns_species2))==model.u_2[biomass_id[i_biomass]])


start_time = time.time()
print('time is: ',start_time)
# solve using the nonlinear solver ipopt
result=SolverFactory('ipopt').solve(model)
end_time = time.time()
execution_time = end_time - start_time
print(value(model.v[biomass_id[0]]),'   ',value(model.v[biomass_id[1]]),'  ',  value(model.obj))
print(f"Execution time: {execution_time:.6f} seconds", '  ', result.solver.status)



