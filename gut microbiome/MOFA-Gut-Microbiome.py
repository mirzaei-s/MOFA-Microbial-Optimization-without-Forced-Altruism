import sys
import pandas as pd
import re
import numpy as np
import pyomo.environ as pyo
from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.util.infeasible import log_infeasible_constraints
import time

Sij = pd.read_excel('D:\CommunityGrowthRate\S_human3.xlsx', 'Sheet1',index_col=0) # stochiometric matrix

m= pd.read_excel('D:\CommunityGrowthRate\S_human3.xlsx', 'Sheet3')
r= pd.read_excel('D:CommunityGrowthRate\S_human3.xlsx', 'Sheet2')
l=pd.read_excel('D:\CommunityGrowthRate\S_human3.xlsx', 'Sheet4')
u=pd.read_excel('D:\CommunityGrowthRate\S_human3.xlsx', 'Sheet5')

rxns=r['rxns'].tolist()
mets=m['mets'].tolist()
lb=l['lb'].tolist()
ub=u['ub'].tolist()

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

#biomass_id=['biomass536_species1','biomass525_species2'] # Biomass human 2
#biomass_id=['biomass223_species1','biomass525_species2'] # Biomass human7
#biomass_id=['biomass525_species1','biomass525_species2'] # Biomass human5
#biomass_id=['biomass223_species1','biomass525_species2'] # Biomass human6
biomass_id=['biomass345_species1','biomass345_species2']

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



