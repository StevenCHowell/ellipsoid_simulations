import sys,os,string
import numpy as np
import time

pwd=os.popen('pwd').readlines()[0].strip() #current directory
sys.path.append(pwd+'/periodictable-1.3.0')

import periodictable as pt
import periodictable.xsf as xxx

import isums
import igolden_vectors

# This program calculates the scattering profile I(q) of a molecule or system of molecules from a pdb file using the golden vector method. 
# Reading the pdb file and assembling the scattering lengths are done in this file, written in python. The calculation of I(q) itself is 
# handed off to a FORTRAN subroutine for efficiency. We use the module 'periodictable' to look up the q-dependent scattering lengths for X-rays. 
# More information on the periodic table cn be found at http://www.reflectometry.org/danse/elements.html.
#
# If you use this code, please cite:
#					MC Watson and JE Curtis. Rapid calculation of small-angle scattering intensities using the golden ratio.
#					J. Appl. Crystallogr. 46 (2013)
#
# The user's parameters are:
#				'option,' :specifies X-ray/neturon and which atoms will be used. 
#				'pdbfile' :the name of the input pdb file
#				'output_file' : the name of the output file
#				'qmax': the maximum magnitude of the scattering vector
#				'nqvalues': number of q-values at which the scattering profile will be calculated
#				'num_golden_vects' : the number of scattering directions to be used in the golden vector method. 

# Below, a longer explanation is provided with each parameter. Excluding the time spent reading the pdb and tabulating
# the scattering lengths, the processing time of the golden vector method scales as nqvals*num_golden_vects*(number of atoms).
#The current program only handles the atoms sulfur, phosphorous, oxygen, nitrogen, carbon, and hydrogen for X-Rays, and also deuterium for neutrons.
# For x-rays, more elements can be easily added by modifying the code around lines 170-190. For neutrons, you can update the dictionary bcoh_atom on line 98.

# option 1: calculate I(q) using X-ray scattering lengths for all atoms. 
# option 2: calculate I(q) using X-ray scattering lengths for all atoms except hydrogens. The scattering lengths of the non-H atoms will be the same as in option 1.
# option 3: calculate I(q) using neutron scattering lengths for all atoms.
# option 4: calculate I(q) using neutron scattering lengths for all atoms except hydrogens. The scattering lengths of the non-H atoms will be the same as in option 3.
# option 5: calculate I(q) using neutron scattering lengths for the amino acids, taking the positions of the C-alpha ('CA') atoms

option=3

#name of the .pdb file used for input

pdb_file='./lysozyme.pdb'

#name of the output file. The first column contains the values of q.
# The second column is the scattering profile on an absolute scale.
# The third column is I(q)/I(0)

output_file='output.dat'

#maximum value of the magnitude of the scattering vector, in units of inverse Angstroms,
# where q=4*pi*sin(theta)/lambda. The minimum value is set to qmin=0.

qmax=0.5

#number of q-values at which the scattering profile will be calculated:

nqvals=50

#The number of scattering directions used to calculate the scattering profile using the golden vector method
# This corresponds to 'n' in our paper. It must be an odd number. 29 is a good place to start.
# If you're interested in q-values beyond about q>0.5 A^(-1), num_golden_vects should be increased for greater accuracy.
# num_golden_vects must be an odd number. It's even, the program will add 1 to make it odd.

num_golden_vects=35

#############################################################################################################################
#### The actual code is below. Unless you want to customize the program, you don't need to change anything beyond this point.
#############################################################################################################################

if(option==1): NEUTRON=0;	CALC_ATOM=1;	CALC_AMINO=0;	NO_H=0
if(option==2): NEUTRON=0;	CALC_ATOM=1;	CALC_AMINO=0;	NO_H=1
if(option==3): NEUTRON=1;	CALC_ATOM=1;	CALC_AMINO=0;	NO_H=0
if(option==4): NEUTRON=1;	CALC_ATOM=1;	CALC_AMINO=0;	NO_H=1
if(option==5): NEUTRON=1;	CALC_ATOM=0;	CALC_AMINO=1;	NO_H=0

if(option>5):
	
	print 'the option must be 1,2,3,4,5'
	sys.exit()
	
if(int(num_golden_vects)%2==0):
	
	print 'num_golden_vects must be an odd number. It was changed to ',num_golden_vects+1
	num_golden_vects += 1
	
nqvals=int(nqvals) # ensure the number is an integer

print 'number of scattering magnitudes=',nqvals
print 'number of golden vectors = ',num_golden_vects

qmin=0 

qlist=np.arange(nqvals)*(qmax-qmin)/(nqvals-1) # evenly spaced on a linear scale

bcoh_amino={'GLY':1.72, 'ALA':1.64, 'VAL':1.47, 'LEU':1.39, 'ILU':1.39, 'ILE':1.39, 'PHE':4.13, 'TYR':4.71, 'TRP':6.02, 'ASP':3.84, 'GLU':3.76,
	'SER':2.22, 'THR':2.14, 'ASN':3.45, 'GLN':3.36, 'LYS':1.58, 'ARG':3.45, 'HIS': 4.96, 'HSE':4.96, 'MET':1.76, 'CYS':1.93, 'PRO':2.22}
# in units of 10^(-12) cm

bcoh_atom={'H': -.37390, 'D':.6671, 'C':.6646, 'S':.2847, 'P':.513, 'N':.936, 'O':.5803}
# in units of 10^(-12) cm

if(CALC_ATOM ==1):  bcoh=bcoh_atom
if(CALC_AMINO==1):  bcoh=bcoh_amino

start_time=time.time()

###### read pdb file:

f=open(pdb_file)
count=0
for line in f.readlines():
    
    words = line.strip().split()
    
    found_line=0
    
    if(line[0:4]=='ATOM'):
        
        atom_name=string.strip(line[12:16])
    
        if(CALC_ATOM==1):
		
		if(NO_H==1 and atom_name[0]!='H'): found_line=1
		
		if(NO_H==0): found_line=1
    
        if(CALC_AMINO==1 and atom_name=='CA'): # take C-alpha

            res_name=string.strip(line[17:21])
            found_line=1


	if (found_line==1):
		
		add_row=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
		#print('residue= ',words[3])
		
		if(CALC_AMINO):		add_scatt_length=bcoh[res_name]
		
		if(CALC_ATOM):		add_scatt_length=bcoh[atom_name[0]] # taking just the first letter
        
		
		if(count==0):
			
			positions=add_row
			batom=np.array([add_scatt_length])
			batom_names=[atom_name[0]]
			count +=1
        
		else:
			positions=np.vstack([ positions,add_row ])
			batom=np.append(batom,add_scatt_length)
			batom_names.append(atom_name[0])
			count+=1

f.close()

xcoor=positions[:,0];   ycoor=positions[:,1];   zcoor=positions[:,2]

natoms=len(batom) # number of atoms to be used in calculation

print 'number of atoms to be used in calculation=',natoms

Iofq=np.zeros(qlist.size) # the scattering intensity

# look up the scattering lengths

print 'tabulating scattering lengths...'

# for x-rays, make a 1D array of scattering lengths for each atom at each value of q
# for example, if the protein consists of the atoms C, N, P, the array reads
#  [C(q1) N(q1) P(q1) C(q2) N(q2) P(q2) ...]

if(NEUTRON==0):
	
	sulfur=pt.S;		phosphorus=pt.P;		oxygen=pt.O
	nitrogen=pt.N;		carbon=pt.C;			hydrogen=pt.H
	
	batom_big=np.zeros(natoms*qlist.size)
	
	for j in range(qlist.size):

		for i in range(natoms):

			if(batom_names[i]=='S'): batom_big[natoms*j+i]=xxx.Xray(sulfur).f0(qlist[j])
			if(batom_names[i]=='P'): batom_big[natoms*j+i]=xxx.Xray(phosphorus).f0(qlist[j])
			if(batom_names[i]=='O'): batom_big[natoms*j+i]=xxx.Xray(oxygen).f0(qlist[j])
			if(batom_names[i]=='N'): batom_big[natoms*j+i]=xxx.Xray(nitrogen).f0(qlist[j])
			if(batom_names[i]=='C'): batom_big[natoms*j+i]=xxx.Xray(carbon).f0(qlist[j])
			if(batom_names[i]=='H'): batom_big[natoms*j+i]=xxx.Xray(hydrogen).f0(qlist[j])

# calculate I(q) using GV method

print 'calculating scattering profile...'

start_time_GV=time.time()

gold_vects=igolden_vectors.golden_vectors2(num_golden_vects)
			
goldx=gold_vects[:,0];	goldy=gold_vects[:,1];	goldz=gold_vects[:,2]

if(NEUTRON==0):
    Iofq = isums.debye_gv_xray(xcoor,ycoor,zcoor,batom_big,qlist,goldx,goldy,goldz,1.0/float(num_golden_vects),num_golden_vects,(qlist.size),natoms,natoms*(qlist.size))

if(NEUTRON==1):
    Iofq = isums.debye_gv(xcoor,ycoor,zcoor,batom,qlist,goldx,goldy,goldz,1.0/float(num_golden_vects),num_golden_vects,(qlist.size),natoms)
		
end_time_GV=time.time()

# write data to file
	
data=np.zeros((qlist.size,3))
	
data[:,0]=qlist
data[:,1]=Iofq
data[:,2]=Iofq/Iofq[0]
	
np.savetxt(output_file,data,delimiter='\t')

print 'time spent reading pdb file and gathering scattering lengths=',start_time_GV-start_time,' seconds'
print 'time spent calculating the scattering profile=',end_time_GV-start_time_GV,' seconds'


