import pandas as pd
import csv

# This script enumerates the pseudo-binary design-space for half-Heusler compounds (XYZ) based only on transition-metal elements. 
# Pseudobinaries based on substitution on two of the three sites are considered. 
# The possible semi-conducting end-member compositions are first enumerated by applying the 18-electron rule.
# CAVEAT: The 18-electron rule holds only for transition metal based elements. In the case of rare-earth based compounds, a more general 'valence-balanced rule should be applied.
# Psudobinary end-member pairs are subsequently enumerated by identifying compounds which have two elemnts in common
 
  
# List of common transition metal based elements for Half Heusler structure along with their number of valence electrons

heusler_elements = {
    'X': {'Li':1, 'Mg':2, 'Ca':2, 'Sc':3, 'Y':3, 'Ti':4, 'Zr':4, 'Hf':4, 'V':5, 'Nb':5, 'Ta':5},
    'Y': {'Fe':8, 'Ru':8, 'Os':8, 'Co':9, 'Rh':9, 'Ir':9, 'Ni':10, 'Pd':10, 'Pt':10, 'Cu':11, 'Ag':11, 'Au':11, 'Zn':12, 'Cd':12},
    'Z': {'Al':3, 'Ga':3, 'In':3, 'Si':4, 'Ge':4, 'Sn':4, 'Pb':4, 'As':5, 'Sb':5, 'Bi':5}
}


df=pd.DataFrame(columns=["Formula_A","Formula_B","Mixing_Elements_on_Mixing_Site1","Mixing_Elements_on_Mixing_Site2","Common_Elements"])


EM=[]					#list for end-member compounds
x=[]					#list for X-site element of end-member compound
y=[]					#list for Y-site element of end-member compound
z=[]					#list for Z-site element of end-member compound


for elementX, vecX in heusler_elements["X"].items():
	for elementY, vecY in heusler_elements["Y"].items():
		for elementZ, vecZ in heusler_elements["Z"].items():
			if vecX+vecY+vecZ == 18:				#Applying the 18-electron rule constraint to filter down the possible number of end members.
				EM.append(elementX+elementY+elementZ)
				x.append(elementX)
				y.append(elementY)
				z.append(elementZ)

				for i in range(len(EM)):				#Iterating over end-members to find instances with common Y and Z-site elements.  
					if (x[i]!=x[len(EM)-1]) and (y[i]!=y[len(EM)-1]) and (z[i]==z[len(EM)-1]):
						df = df.append({"Formula_A":EM[i],"Formula_B":EM[len(EM)-1],"Mixing_Elements_on_Mixing_Site1":[x[i],x[len(EM)-1]],"Mixing_Elements_on_Mixing_Site2":[y[i],y[len(EM)-1]],"Common_Elements":[z[i]]},ignore_index = True)

				for i in range(len(EM)):				#Iterating over end-members to find instances with common X and Z-site elements.
					if (x[i]!=x[len(EM)-1]) and (y[i]==y[len(EM)-1]) and (z[i]!=z[len(EM)-1]):
						df = df.append({"Formula_A":EM[i],"Formula_B":EM[len(EM)-1],"Mixing_Elements_on_Mixing_Site1":[x[i],x[len(EM)-1]],"Mixing_Elements_on_Mixing_Site2":[z[i],z[len(EM)-1]],"Common_Elements":[y[i]]},ignore_index = True)

				for i in range(len(EM)):				#Iterating over end-members to find instances with common X and Y-site elements.
					if (x[i]==x[len(EM)-1]) and (y[i]!=y[len(EM)-1]) and (z[i]!=z[len(EM)-1]):
						df = df.append({"Formula_A":EM[i],"Formula_B":EM[len(EM)-1],"Mixing_Elements_on_Mixing_Site1":[y[i],y[len(EM)-1]],"Mixing_Elements_on_Mixing_Site2":[z[i],z[len(EM)-1]],"Common_Elements":[x[i]]},ignore_index = True)


df.to_csv('pseudobinary_design_space_DOUBLE.csv',sep = ',', index=False)
