#
# Hailiang Zhang
# NIST & UTK
#

import sys,os,locale,glob
import random
import numpy as np

import sassie.sasmol.sasmol as sasmol
import sascalc

pwd=os.popen('pwd').readlines()[0].strip() #current directory
sys.path.append(pwd+'/periodictable-1.3.0')

import periodictable as pt
import periodictable.xsf as xxx

def get_atomic_radii(mol, scale=1.0):
    '''
    Functionality: get the radii for each individual atom
    Input: mol--sasmol object
           scale--the scaling factor to be applied on the VDW radius (default 1.0)
    Return: a list of radii for each atoms
    Note: VDW radius is used
    '''
    R_VDW={'H':1.2*scale, 'C':1.7*scale, 'N':1.55*scale, 'O':1.52*scale, 'S':1.8*scale, 'P':1.8*scale} # vdw radius
    radii = [R_VDW[atom_name[0]] for atom_name in mol.name()]
    return radii

def calculate_solvent_sld(fraction_D2O):
    '''
    Functionality: calculate the H2O/D2O solvent scattering length density based on the fraction of D2O
    Input: fraction_D2O--fraction of D2O
    Return: solvent scattering length density
    Note: room temperature is assumed
    '''
    fraction_H2O = 1.0-fraction_D2O
    b_H2O = -.37406*2+.5804
    b_D2O = .6671*2+.5804
    b_solvent = fraction_H2O*b_H2O + fraction_D2O*b_D2O
    density_solvent = fraction_H2O*1. + fraction_D2O*1.1
    mw_solvent = fraction_H2O*18. + fraction_D2O*20.
    sld_solvent = b_solvent*(density_solvent/mw_solvent*6.02e23)/1.e24
    return sld_solvent

def calculate_solvent_electron_density(fraction_D2O):
    '''
    Functionality: calculate the H2O/D2O solvent electron density based on the fraction of D2O
    Input: fraction_D2O--fraction of D2O
    Return: solvent electron density (#/A^3)
    Note: room temperature is assumed
    '''
    fraction_H2O = 1.0-fraction_D2O
    electrons_of_H2O = 10.0
    density_solvent = fraction_H2O*1. + fraction_D2O*1.1
    mw_solvent = fraction_H2O*18. + fraction_D2O*20.
    electron_density_of_solvent = electrons_of_H2O*(density_solvent/mw_solvent*6.02e23)/1.e24
    return electron_density_of_solvent

def find_non_aliphatic_Hs(mol):
    '''
    Functionality: find the non-aliphatic Hs for a molecule
    Input: mol--sasmol object
    Return: list of masks for the non-aliphatic Hs (0--is NOT a non-aliphatic H; 1--is a non-aliphatic H)
    Note: it is assumed that Hs always follows its attached heavy atoms in the pdb file
    '''
    natoms = mol.natoms()
    mask_non_aliphatic_Hs = [0]*natoms
    for i in range(natoms):
        atom_full_name = mol.name()[i]
        if atom_full_name[0] == 'H':
            if previous_nonH_atom_short_name!='C':
                 mask_non_aliphatic_Hs[i] = 1
            else:
                 mask_non_aliphatic_Hs[i] = 0
        else:
            mask_non_aliphatic_Hs[i] = 0
            previous_nonH_atom_short_name = atom_full_name[0]
    return mask_non_aliphatic_Hs

def get_atomic_neutron_bs(mol, HD_exchange_ratio, ofile, deuteration_option="all"):
    '''
    Functionality: deuterate the exchangeable H's based on the HD exhcange ratio in various ways
    Input: mol--sasmol object
           HD_exchange_ratio--H/D exchange ratio 
           deuteration_option="all"--deuterate all non-aliphatic H's including buried (default)
           deuteration_option="none"--no deuteration
           deuteration_option="surface"--deuterate all surface non-aliphatic H's (not implemented yet)
           deuteration_option="cryson"--deuterate by cryson default (backbone HN exchanged at 0.9*D2O%, sidechain HN exchanged at 1.*D2O%) (not implemented yet)
    Return: list of scattering length for each atom with H's deutrated based on the input deuteration_option
    Note: this method is named in consistent with the CPP code and may be misleading, because there is no atomic atom type changed from H to D in the sasmol object
          it just returns the scattering length of all the atoms with H's properly deuterated
    '''
    B={'H': -0.37390, 'D': 0.6671, 'C': 0.6646, 'S': 0.2847, 'P': 0.513, 'N': 0.936, 'O': 0.5803}
    mask_non_aliphatic_Hs = find_non_aliphatic_Hs(mol)
    natoms = mol.natoms()
    b = [B[atom_name[0]] for atom_name in mol.name()]
    count=0
    of = open(ofile,'w')
    if deuteration_option=="none": # no deuteration
        return b # b was already initialized as no-deuteration
    elif deuteration_option=="all": #deuterate all non-aliphatic H's including buried
        for i in range(natoms):
            if mask_non_aliphatic_Hs[i] and random.random()<HD_exchange_ratio:
            #if mask_non_aliphatic_Hs[i] and count<(mask_non_aliphatic_Hs.count(1))*HD_exchange_ratio: # sequentially deuterate Hs for debugging purpose only
                of.write("Deuterate atom #%5d"%(i+1)+" (atom name: %5s"%(mol.name()[i])+" in residue: %5s"%(mol.resname()[i])+")\n")
                b[i] = B['D']
                count += 1
            else:
                b[i] = B[mol.name()[i][0]]
        return b
    else:
        print("deuteration_option '"+deuteration_option+"' not implemented yet!")
        exit(0)

def get_atomic_Xray_fq(mol, qlist):
    '''
    Functionality: get the X-ray f(q) for each individual atom
    Input: mol--sasmol object
    Return: a list of fs at different q values for each atoms
    '''
    natoms = mol.natoms()
    
    sulfur=pt.S;        phosphorus=pt.P;        oxygen=pt.O
    nitrogen=pt.N;        carbon=pt.C;            hydrogen=pt.H
    
    B=np.zeros(natoms*qlist.size)
    
    elelment_names = [atom_name[0] for atom_name in mol.name()]

    for j in range(qlist.size):

        for i in range(natoms):

            if(elelment_names[i]=='S'): B[natoms*j+i]=xxx.Xray(sulfur).f0(qlist[j])
            if(elelment_names[i]=='P'): B[natoms*j+i]=xxx.Xray(phosphorus).f0(qlist[j])
            if(elelment_names[i]=='O'): B[natoms*j+i]=xxx.Xray(oxygen).f0(qlist[j])
            if(elelment_names[i]=='N'): B[natoms*j+i]=xxx.Xray(nitrogen).f0(qlist[j])
            if(elelment_names[i]=='C'): B[natoms*j+i]=xxx.Xray(carbon).f0(qlist[j])
            if(elelment_names[i]=='H'): B[natoms*j+i]=xxx.Xray(hydrogen).f0(qlist[j])

    return B



"""
The following is for devoloper testing/debugging purpose
"""
if __name__=='__main__':

    # setup the samol object
    mol = sasmol.SasMol(0)
    pdb_filename = 'new_lysozyme.pdb'
    mol.read_pdb(pdb_filename)
    natoms = mol.natoms()
    nframes = mol.number_of_frames()

    # get the Xray fs
    sascalc = sascalc.SasCalc(20, 0.5)
    qlist = sascalc.get_Q()
    B = get_atomic_Xray_fq(mol, qlist)
    import pprint
    pprint.pprint(B)

    # get solvent D2O fraction & HD exchange ratio
    fraction_D2O = 1.0
    HD_exchange_ratio = 0.1

    # calculate solvent sld
    sld_solvent = calculate_solvent_sld(fraction_D2O)

    # get the atomic radii
    radii = get_atomic_radii(mol)

    # deuterate the exchangeable H's
    deuteration_option = "all"
    b = get_atomic_neutron_bs(mol, HD_exchange_ratio, "HDexchange_Info_D2O%%=%d_HDexchange%%=%d.txt"%(int(fraction_D2O*100),int(HD_exchange_ratio*100)), deuteration_option)

    # print summary
    print("\n=====================================================================\nCondition summary")
    print("D2O content: %4.1f %%"%(fraction_D2O*100))
    print("H/D exchange ratio: %4.1f"%HD_exchange_ratio)
    print("sld_solvent: %10.6f"%sld_solvent)
