#
# Hailiang Zhang
# NIST & UTK
#

import sys,os,locale

import sassie.sasmol.sasmol as sasmol

import interface
import atomic_info
import sascalc

if __name__=="__main__":

    # parse the input parameter file
    xon, input_folder, pdb_filename, fraction_D2O, HD_exchange_ratio, deuteration_option, num_Q, Qmax, sascalc_option, gv_option, num_golden_vectors, starting_num_golden_vectors, max_num_golden_vectors, tolerance, output_folder = interface.parse_input_parameter_file(sys.argv[1])

    # setup the samol object
    mol = sasmol.SasMol(0)
    mol.read_pdb(os.path.join(input_folder,pdb_filename))
    natoms = mol.natoms()
    frame = 0

    # calculate solvent sld
    if (xon=='neutron'):
        sld_solvent = atomic_info.calculate_solvent_sld(fraction_D2O)
    elif (xon=='xray'):
        sld_solvent = atomic_info.calculate_solvent_electron_density(fraction_D2O)

    # get the atomic radii
    radii = atomic_info.get_atomic_radii(mol)

    # setup sascalc
    sascalc = sascalc.SasCalc(num_Q, Qmax)

    # get the atomic scattering information
    if xon=='xray':
        b = atomic_info.get_atomic_Xray_fq(mol, sascalc.get_Q())
    elif xon=='neutron':
        b = atomic_info.get_atomic_neutron_bs(mol, HD_exchange_ratio, os.path.join(output_folder,"HDexchange_Info_D2O%%=%d_HDexchange%%=%d.txt"%(int(fraction_D2O*100),int(HD_exchange_ratio*100))), deuteration_option)

    # print summary
    print("\n=====================================================================\nCondition summary")
    print
    print("D2O content: %4.1f %%"%(fraction_D2O*100))
    print("H/D exchange ratio: %4.1f"%HD_exchange_ratio)
    print("Deuteration option: "+deuteration_option)
    print("Solvent scattering length density: %10.6f"%sld_solvent)
    print
    print("Number of q values: %d"%num_Q)
    print("Maximum q values: %f"%Qmax)
    print
    print("Number of golden vectors: %d"%num_golden_vectors)
    print("Starting number of golden vectors: %d"%starting_num_golden_vectors)
    print("Maximum number of golden vectors: %d"%max_num_golden_vectors)
    print("Tolerance value of GVVV runtime average convergence: %f"%tolerance)
    print
    print("Results folder: "+output_folder)
    print("=====================================================================\n")

    # run sascalc
    print("Running sascalc of "+sascalc_option+" based on golden vector method with fixed number of golden vectors...\n")
    if gv_option=='fixed':
        sascalc.calc_GV_fix_n(mol, frame, num_golden_vectors, xon=xon, option=sascalc_option, B=b, R=radii, sld_solvent=sld_solvent)
    elif gv_option=='converge':
        sascalc.calc_GV_converge_n(mol, frame, starting_num_golden_vectors, max_num_golden_vectors, tolerance, xon=xon, option=sascalc_option, B=b, R=radii, sld_solvent=sld_solvent)
    sascalc.write_to_file(os.path.join(output_folder,"Iq_"+xon+"_"+sascalc_option+"_"+gv_option+"_D2O%%=%d_HDexchange%%=%d.txt"%(int(fraction_D2O*100),int(HD_exchange_ratio*100))))
