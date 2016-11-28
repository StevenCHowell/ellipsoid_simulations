#
# Hailiang Zhang
# NIST & UTK
#

import sys
import numpy
import GV

class SasCalc:
    '''
    Home of Q, I, and SAS calculators
    '''

    def __init__(self, num_Q=50, Qmax=0.5):
        '''
        Functionality: initialize Q, I arrays
        Input: num_Q -- number of 'q' values; Qmax -- qmax value
        Return: None
        Note: a number of num_Q 'q' values will be equally placed between [0,Qmax]
        '''
        self.__Qmax = Qmax
        self.__num_Q = num_Q
        self.__setup_Q()
        self.__setup_I()

    def get_I(self):
        '''
        Functionality: get the I array
        Input: None
        Return: I array
        '''
        return self.__I

    def get_Q(self):
        '''
        Functionality: get the Q array
        Input: None
        Return: Q array
        '''
        return self.__Q

    def get_num_Q(self):
        '''
        Functionality: get number of 'q' values
        Input: None
        Return: number of 'q' values
        '''
        return self.__num_Q

    def get_Qmax(self):
        '''
        Functionality: get max 'q' values
        Input: None
        Return: max 'q' values
        '''
        return self.__Qmax

    def write_to_file(self, file):
        '''
        Functionality: write Q and I arrays to the external file
        Input: file -- external file name
        Return: None
        Note: two columes ('q' and 'I(q)') will be written to the file
        '''
        of = open(file,'w')
        for q,i in zip(self.get_Q(), self.get_I()):
            of.write("%f %f\n"%(q,i))
        of.close()


    def calc_GV_fix_n(self, mol, frame=0, num_golden_vectors=11, **additional_parameters):
        '''
        Functionality: calculate I(q) based on golden vector method for a fixed number of golden vectors
        Input: mol -- sasmol object
               frame -- frame number (default 0)
               num_golden_vectors -- number of golden vectors (default 51) (an even number input will be reset to be the next odd number)
               additional_parameters:
                   xon="xray" -- x-ray scattering calculation
                   xon="neutron" -- neutron scattering calculation
                   option="vacuum" -- vacuum scattering calculation
                   option="solvent" -- solvent scattering calculation
                   option="complete" -- vacuum+solvent scattering calculation
                   R -- the list of atomic radii (if not provided for option="solvent" or "complete", program will quit)
                   B -- the list of atomic scattering length (if not provided for option="vacuum" or "complete", program will quit)
                   sld_solvent -- solvent scattering length density (if not provided for option="solvent" or "complete", program will quit)
        Return: None
        '''
        xon, option, R, B, sld_solvent = self.__parse_additional_parameters(additional_parameters)
        #
        if num_golden_vectors == num_golden_vectors/2*2:
            print('=============================\nWarning:')
            print('The number of golden vectors should be an odd number')
            print('so we will change it from %d to %d'%(num_golden_vectors,num_golden_vectors+1)) 
            num_golden_vectors += 1
        num_Q = self.__num_Q
        Q = self.get_Q()
        I = self.get_I()
        gv = GV.GV(num_golden_vectors)
        golden_vectors = gv.get_golden_vectors()
        natoms = mol.natoms()
        coor = mol.coor()
        for i in range(num_Q):
            I[i]=0.0
            for j in range(num_golden_vectors):
                qx = Q[i]*golden_vectors[j,0]
                qy = Q[i]*golden_vectors[j,1]
                qz = Q[i]*golden_vectors[j,2]
                A_real = 0.0
                A_imag = 0.0
                for k in range(natoms):
                    q_dot_r = qx*coor[frame,k,0]+qy*coor[frame,k,1]+qz*coor[frame,k,2]
                    if xon=='neutron':
                        bvac = B[k]
                    elif xon=='xray':
                        bvac = B[natoms*i+k]
                    if option=='vacuum':
                        b = bvac
                    elif option=='solvent':
                        b = -sld_solvent * (4./3.*numpy.pi*R[k]*R[k]*R[k]) 
                    elif option=='complete':
                        b = bvac - sld_solvent * (4./3.*numpy.pi*R[k]*R[k]*R[k]) 
                    else:
                        print("Unknown option: ",option)
                        sys.exit(0)
                    A_real += b*numpy.cos(q_dot_r)
                    A_imag += b*numpy.sin(q_dot_r)
                I[i] += A_real*A_real + A_imag*A_imag
        I /= num_golden_vectors 

    def calc_GV_converge_n(self, mol, frame=0, start_num_golden_vectors=11, max_num_golden_vector=151, tolerance=0.01, **additional_parameters):
        '''
        Functionality: calculate I(q) based on golden vector method by converging the runtime avarage
        Input: mol -- sasmol object
               frame -- frame number (default 0)
               start_num_golden_vectors -- starting number of golden vectors (default 11)
               max_num_golden_vectors -- maximum number of golden vectors (default 151)
               tolerance -- the smallest value of difference between the runtime avearge that will allow the convergence running to quit
               additional_parameters:
                   xon="xray" -- x-ray scattering calculation
                   xon="neutron" -- neutron scattering calculation
                   option="vacuum" -- vacuum scattering calculation
                   option="solvent" -- solvent scattering calculation
                   option="complete" -- vacuum+solvent scattering calculation
                   R -- the list of atomic radii (if not provided for option="solvent" or "complete", program will quit)
                   B -- the list of atomic scattering length (if not provided for option="vacuum" or "complete", program will quit)
                   sld_solvent -- solvent scattering length density (if not provided for option="solvent" or "complete", program will quit)
        Return: None
        '''
        xon, option, R, B, sld_solvent= self.__parse_additional_parameters(additional_parameters)
        #
        num_golden_vectors = start_num_golden_vectors
        Q = self.get_Q()
        I = self.get_I()
        converged = False
        count = 0
        while not converged and num_golden_vectors<max_num_golden_vector:
            print "Number of golden vectors: ",num_golden_vectors
            count += 1
            self.calc_GV_fix_n(mol, frame, num_golden_vectors, xon=xon, B=B, R=R, sld_solvent=sld_solvent)
            if num_golden_vectors == start_num_golden_vectors:
                I_run_average = numpy.copy(I)
            else:
                I_run_average_previous = numpy.copy(I_run_average)
                I_run_average = (I_run_average_previous*(count-1)+I)/count
                diff=numpy.mean(numpy.abs(I_run_average-I_run_average_previous)/I_run_average_previous)
                if(diff<tolerance): 
                    converged=True
            num_golden_vectors+=2
        I = I_run_average

    def __setup_Q(self):
        '''
        Functionality: set the Q array
        Input: None
        Return: Q array
        '''
        num_Q = self.get_num_Q()
        Qmax = self.get_Qmax()
        self.__Q = numpy.zeros(num_Q)
        for i in range(num_Q):
            self.__Q[i] = i*(Qmax/(num_Q-1));
        return

    def __setup_I(self):
        '''
        Functionality: set the I array
        Input: None
        Return: None
        '''
        num_Q = self.get_num_Q()
        self.__I = numpy.zeros(num_Q)
        return

    def __parse_additional_parameters(self, additional_parameters):
        '''
        Functionality: parse the additional input parameters for golden vector based sas calculation
        Input: additional_parameters:
                   xon="xray" -- x-ray scattering calculation
                   xon="neutron" -- neutron scattering calculation
                   option="vacuum" -- vacuum scattering calculation
                   option="solvent" -- solvent scattering calculation
                   option="complete" -- vacuum+solvent scattering calculation
                   R -- the list of atomic radii (if not provided for option="solvent" or "complete", program will quit)
                   B -- the list of atomic scattering length (if not provided for option="vacuum" or "complete", program will quit)
                   sld_solvent -- solvent scattering length density (if not provided for option="solvent" or "complete", program will quit)
        Return: option, R, B, sld_solvent
        '''
        # parse x-ray or neutron input
        if 'xon' not in additional_parameters.keys():
            xon = 'neutron'
        else:
            xon = additional_parameters['xon']
        # parse option input
        if 'option' not in additional_parameters.keys():
            option = 'complete'
        else:
            option = additional_parameters['option']
        # parse atomic radii input
        if 'R' not in additional_parameters.keys():
            if option in ['solvent','complete']:
                print('atomic radii was not provided when solvent/complete scattering calculation is requested')
                sys.exit(0)
            else:
                R=None
        else:
            R = additional_parameters['R']
        # pass atomic scattering length input
        if 'B' not in additional_parameters.keys():
            if option in ['vacuum','complete']:
                print('atomic scattering length was not provided when vacuum/complete scattering calculation is requested')
                sys.exit(0)
            else:
                B=None
        else:
            B = additional_parameters['B']
        # pass solvent scattering length density input
        if 'sld_solvent' not in additional_parameters.keys():
            if option in ['solvent','complete']:
                print('solvent scattering length density was not provided when solvent/complete scattering calculation is requested')
                sys.exit(0)
            else:
                sld_solvent = None
        else:
            sld_solvent = additional_parameters['sld_solvent']
        # return the parsed parameters
        return xon, option, R, B, sld_solvent


"""
The following is for devoloper testing/debugging purpose
"""
if __name__=='__main__':
    import atomic_info
    import sassie.sasmol.sasmol as sasmol

    sascalc = SasCalc(20,0.5)
    print "Q:\n",sascalc.get_Q()
    print "I:\n",sascalc.get_I()

    # setup the samol object
    mol = sasmol.SasMol(0)
    pdb_filename = 'new_lysozyme.pdb'
    mol.read_pdb(pdb_filename)
    natoms = mol.natoms()
    nframes = mol.number_of_frames()

    # get solvent D2O fraction & HD exchange ratio
    fraction_D2O = 0.0
    HD_exchange_ratio = 0.0

    # calculate solvent sld
    sld_solvent = atomic_info.calculate_solvent_sld(fraction_D2O)
    print 'sld_solvent: ',sld_solvent

    # get the atomic radii
    radii = atomic_info.get_atomic_radii(mol)
    print 'radii (first 10 atoms): ',radii[:10]

    # deuterate the molecule
    deuteration_option = "all"
    b = atomic_info.get_atomic_neutron_bs(mol, HD_exchange_ratio, "HDexchange_Info_D2O%%=%d_HDexchange%%=%d.txt"%(int(fraction_D2O*100),int(HD_exchange_ratio*100)), deuteration_option)
    print 'b (first 10 atoms): ',b[:10]

    # setup sascalc
    sascalc.calc_GV_fix_n(mol, 0, 11, option='vacuum', B=b, R=radii, sld_solvent=sld_solvent)
    print "I after calc_GV_fix_n in vacuum: ",sascalc.get_I()
    sascalc.calc_GV_fix_n(mol, 0, 11, option='solvent', B=b, R=radii, sld_solvent=sld_solvent)
    print "I after calc_GV_fix_n in solvent: ",sascalc.get_I()
    sascalc.calc_GV_fix_n(mol, 0, 11, option='complete', B=b, R=radii, sld_solvent=sld_solvent)
    print "I after calc_GV_fix_n in complete: ",sascalc.get_I()
    #sascalc.calc_GV_converge_n(mol, 0, 11, 51, 0.001, option='complete', B=b, R=radii, sld_solvent=sld_solvent)
    #print "I after calc_GV_converge_n in complete: ",sascalc.get_I()


    # test x-ray scattering
    b = atomic_info.get_atomic_Xray_fq(mol, sascalc.get_Q())
    sascalc.calc_GV_fix_n(mol, 0, 11, xon='xray',option='vacuum', B=b, R=radii, sld_solvent=sld_solvent)
    print "I after calc_GV_fix_n in vacuum for X-ray scattering: ",sascalc.get_I()
