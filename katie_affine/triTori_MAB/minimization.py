'''
This file contains modules for affine minimization
Reference: ...
Hailiang Zhang
Jan 2015
'''
import sys,os,string,locale,numpy
import numpy

import shapes
import affine 
import fileIO
#import cryson_driver
import GV

from scipy import optimize

import sassie.sasmol.sasmol as sasmol
#import matplotlib.pylab as plt

if len(sys.argv)==1: NGPU=0
else: NGPU = locale.atoi(sys.argv[1])


def calcX2(parameters, Nq, qmax, Q,Iq_exp, Iq_model, Err_exp, bs, coordinate_matrix):
    '''
    reduced chi-square as defined in Curtis CPC 2011.
    '''

    Nq = len(Iq_exp)
    Npoints = len(bs)

    affine_transformation_matrix = affine.buildAffineMatrix(parameters[0], parameters[1], parameters[2], 0.0, parameters[3], parameters[4], parameters[5], 0.0, parameters[6], parameters[7], parameters[8], 0.0)
    affine.fillShapeByReverseAffineTransformation(Npoints, affine_transformation_matrix, shapes.sphere_equation, coordinate_matrix)

    GV.GV(coordinate_matrix, bs, I_model, Npoints, Nq, qmax, NGPU)
    Iq_model /= Iq_model[0]
    #plt.semilogy(Q,Iq_model)
    #print 'Iq_model: \n',Iq_model

    X2 = 0
    diff=Iq_model-Iq_exp
    diff2 = diff*diff
    X2=numpy.sum(diff2/Err_exp)/(Nq-1) 
    print 'parameters: '
    numpy.savetxt(sys.stdout, numpy.array(parameters).reshape((3,3)), '%5.2f')
    print 'X2: ',X2

    return X2

 

if __name__ == '__main__':

    # GV calculation on GPU
    o = sasmol.SasMol(0)
    o.read_pdb(os.path.join('data',os.getcwd().split('/')[-1]+'.pdb'))
    Nq = 50
    qmax = 0.3
    Q = numpy.linspace(0,qmax,Nq)
    I_exp = numpy.zeros(Nq)
    B={'H': -0.37390, 'D': 0.6671, 'C': 0.6646, 'S': 0.2847, 'P': 0.513, 'N': 0.936, 'O': 0.5803}
    bs = map(lambda x: B[x], [item[0] for item in o.element()])
    GV.GV(o._coor[0].T, bs, I_exp, o.natoms(), Nq, qmax)
    I_exp /= I_exp[0]
    Err_exp =I_exp/100.
    #plt.plot(Q,I_exp,label='exp')
    #plt.legend(loc='best')
    #plt.show()
    print 'I exp:\n',I_exp
    fileIO.writeIq('output/Iq_exp.txt',Q,I_exp)
    exit(0)

    # run minimization 
    Npoints = 10000
    parameters = (20.0, 0.0, 0.0,  0.0,20.0,0.0, 0.0,0.0,20.0)
    bs = numpy.ones(Npoints)
    I_model = numpy.zeros(Nq)
    coordinate_matrix = numpy.matrix(numpy.zeros((3,Npoints)))
    calcX2(parameters, Nq, qmax, Q, I_exp, I_model, Err_exp, bs, coordinate_matrix)
    fileIO.make_and_write_pdb(coordinate_matrix, 'output', 'model_start.pdb')
    fileIO.writeIq('output/Iq_model_start.txt',Q,I_model)

    results = optimize.minimize(calcX2, parameters, args=(Nq, qmax, Q, I_exp, I_model, Err_exp, bs, coordinate_matrix),method='Powell')
    print(results)

    fileIO.make_and_write_pdb(coordinate_matrix, 'output', 'model_final.pdb')
    fileIO.writeIq('output/Iq_model_final.txt',Q,I_model)


