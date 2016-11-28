#
# Hailiang Zhang
# NIST & UTK
#

import sys,os,locale

import sassie.sasmol.sasmol as sasmol

import atomic_info
import sascalc

def parse_input_parameter_file(file):

    xon = None
    input_folder = None
    pdb_filename = None
    fraction_D2O = None
    HD_exchange_ratio = None
    deuteration_option = None
    num_Q = None
    Qmax = None
    sascalc_option = None
    gv_option = None
    num_golden_vectors = None
    starting_num_golden_vectors = None
    max_num_golden_vectors = None
    tolerance = None
    output_folder = None

    # parse the input parameter file
    for line in open(file):
        words = line.split("#")
        if len(words) != 2:
            continue
        if words[1].strip()=="x-ray or neutron":
            xon = words[0].strip()
        elif words[1].strip()=="input folder":
            input_folder = words[0].strip()
        elif words[1].strip()=="input pdb filename":
            pdb_filename = words[0].strip()
        elif words[1].strip()=="D2O fraction":
            fraction_D2O = locale.atof(words[0])
        elif words[1].strip()=="HD exchange ratio":
            HD_exchange_ratio = locale.atof(words[0])
        elif words[1].strip()=="deuteration option":
            deuteration_option = words[0].strip()
        elif words[1].strip()=="number of q values":
            num_Q = locale.atoi(words[0])
        elif words[1].strip()=="maximum q value":
            Qmax = locale.atof(words[0])
        elif words[1].strip()=="sascalc option":
            sascalc_option = words[0].strip()
        elif words[1].strip()=="golden vector method option":
            gv_option = words[0].strip()
        elif words[1].strip()=="number of golden vectors":
            num_golden_vectors = locale.atoi(words[0])
        elif words[1].strip()=="starting number of golden vectors":
            starting_num_golden_vectors = locale.atoi(words[0])
        elif words[1].strip()=="maximum number of golden vectors":
            max_num_golden_vectors = locale.atoi(words[0])
        elif words[1].strip()=="tolerance value of GVVV runtime average convergence":
            tolerance = locale.atof(words[0])
        if words[1].strip()=="results folder":
            output_folder = words[0].strip()

    # set default values
    if xon == None:
        xon = 'neutron'
    if input_folder == None:
        input_folder = './'
    if output_folder == None:
        output_folder = './'
    if num_Q == None:
        num_Q = 50
    if Qmax == None:
        Qmax = 0.5
    if sascalc_option == None:
        sascalc_option = 'vacuum'
    if gv_option == None:
        gv_option = 'fixed'

    if gv_option == 'fixed':
        if num_golden_vectors == None:
            num_golden_vectors = 35
    if gv_option == 'converge':
        if starting_num_golden_vectors == None:
            starting_num_golden_vectors = 11
        if  max_num_golden_vectors == None:
            max_num_golden_vectors = 149
        if tolerance == None:
            tolerance = 0.001

    # check for missing parameters
    if pdb_filename == None:
        print '****************************************************************************************'
        print '"input pdb filename" must be provided in the input parameter file. Program will quit!'
        print
        exit(0)
    if xon == 'neutron':
        if fraction_D2O == None:
            print '****************************************************************************************'
            print '"D2O fraction" must be provided in the input parameter file. Program will quit!'
            print
            exit(0)
        if HD_exchange_ratio == None:
            print '****************************************************************************************'
            print '"HD exchange ratio" must be provided in the input parameter file. Program will quit!'
            print
            exit(0)

    # check for redandencies
    if xon == 'xray':
        if fraction_D2O != None:
            print '****************************************************************************************'
            print 'Warning: "D2O fraction" will be ignored for X-ray scattering'
            print
        if HD_exchange_ratio != None:
            print '****************************************************************************************'
            print 'Warning: "HD exchange ratio" will be ignored for X-ray scattering'
            print
    if gv_option == 'fixed':
        if starting_num_golden_vectors!=None:
            print '****************************************************************************************'
            print 'Warning: "starting number of golden vectors" will be ignored for "golden vector method option" as fixed'
            print
        if max_num_golden_vectors!=None:
            print '****************************************************************************************'
            print 'Warning: "maximum number of golden vectors" will be ignored for "golden vector method option" as fixed'
            print
        if tolerance!=None:
            print '****************************************************************************************'
            print 'Warning: "tolerance value of GVVV runtime average convergence" will be ignored for "golden vector method option" as fixed'
            print
    if gv_option == 'converge':
        if num_golden_vectors != None:
            print '****************************************************************************************'
            print 'Warning: "number of golden vectors" will be ignored for "golden vector method option" as converge'
            print

    # check/setup folder
    if not os.path.isdir(input_folder):
        print '****************************************************************************************'
        print 'Folder "input_folder" does not exist as "input folder". Program will quit!'
        print
        exit(0)
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    return xon, input_folder, pdb_filename, fraction_D2O, HD_exchange_ratio, deuteration_option, num_Q, Qmax, sascalc_option, gv_option, num_golden_vectors, starting_num_golden_vectors, max_num_golden_vectors, tolerance, output_folder
