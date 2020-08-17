import pandas as pd
import csv

# This script enumerates the pseudo-binary design-space for half-Heusler compounds (XYZ). 
# Only pseudobinaries based on mixing on one of the three sites are considered. 
# The possible semi-conducting end-member compositions are first enumerated by applying the Valence balanced rule.
 
  
# List of common transition metal based elements for Half Heusler structure along with their number of valence electrons

heusler_elements = {
    'X': {'Li':1, 'Mg':2, 'Ca':2, 'Sc':3, 'Y':3, 'Ti':4, 'Zr':4, 'Hf':4, 'V':5, 'Nb':5, 'Ta':5, 'La':3, 'Ce':3, 'Pr':3, 'Nd':3, 'Sm':3, 'Gd':3, 'Tb':3, 'Dy':3, 'Ho':3, 'Er':3, 'Tm':3, 'Yb':3, 'Lu': 3},
    'Y': {'Fe':-2, 'Ru':-2, 'Os':-2, 'Co':-1, 'Rh':-1, 'Ir':-1, 'Ni':0, 'Pd':0, 'Pt':0, 'Cu':1, 'Ag':1, 'Au':1, 'Zn':2, 'Cd':2},
    'Z': {'Al':-5, 'Ga':-5, 'In':-5, 'Si':-4, 'Ge':-4, 'Sn':-4, 'Pb':-4, 'As':-3, 'Sb':-3, 'Bi':-3}
}


df=pd.DataFrame(columns=["formula_A","formula_B","mixing_elements","common_elements"])
df1=pd.DataFrame(columns=["End-members"])

EM=[]					#list for end-member compounds
x=[]					#list for X-site element of end-member compound
y=[]					#list for Y-site element of end-member compound
z=[]					#list for Z-site element of end-member compound


for elementX, vecX in heusler_elements["X"].items():
	for elementY, vecY in heusler_elements["Y"].items():
		for elementZ, vecZ in heusler_elements["Z"].items():
			if vecX+vecY+vecZ == 0:				#Applying the 18-electron rule constraint to filter down the possible number of end members.
				EM.append(elementX+elementY+elementZ)
				df1=df1.append({"End-members":elementX+elementY+elementZ},ignore_index = True)
				x.append(elementX)
				y.append(elementY)
				z.append(elementZ)

				for i in range(len(EM)):				#Iterating over end-members to find instances with common Y and Z-site elements.  
					if (i!=len(EM)-1) and (y[i]==y[len(EM)-1]) and (z[i]==z[len(EM)-1]):
						df = df.append({"formula_A":EM[i],"formula_B":EM[len(EM)-1],"mixing_elements":[x[i],x[len(EM)-1]],"common_elements":[y[i],z[i]]},ignore_index = True)

				for i in range(len(EM)):				#Iterating over end-members to find instances with common X and Z-site elements.
					if (i!=len(EM)-1) and (x[i]==x[len(EM)-1]) and (z[i]==z[len(EM)-1]):
						df = df.append({"formula_A":EM[i],"formula_B":EM[len(EM)-1],"mixing_elements":[y[i],y[len(EM)-1]],"common_elements":[x[i],z[i]]},ignore_index = True)

				for i in range(len(EM)):				#Iterating over end-members to find instances with common X and Y-site elements.
					if (i!=len(EM)-1) and (x[i]==x[len(EM)-1]) and (y[i]==y[len(EM)-1]):
						df = df.append({"formula_A":EM[i],"formula_B":EM[len(EM)-1],"mixing_elements":[z[i],z[len(EM)-1]],"common_elements":[x[i],y[i]]},ignore_index = True)

df1.to_csv('end_member_Including_Lanthides.csv',sep = ',')
df.to_csv('pseudobinary_design_space_Including_Lanthnides.csv',sep = ',', index=False)
