'''
This file contains modules for affine minimization
Reference: ...
Hailiang Zhang
Jan 2015
'''
import sys,os,glob,string,locale,numpy
import time

sys.path.append('/home/hailiang/work/LibrariesSources/python')
from affine import affine 
from affine import fileIO
from affine import GV

from affine.shapes import hollowCylinder_equation as shape_equation

from scipy import optimize

import sassie.sasmol.sasmol as sasmol
#import matplotlib.pylab as plt

if len(sys.argv)==1: NGPU=0
else: NGPU = locale.atoi(sys.argv[1])


def calcX2(parameters, Nq, qmax, Q,Iq_exp, Iq_model, Err_exp, bs, coordinate_matrix):

    Nq = len(Iq_exp)
    Npoints = len(bs)

    count = 0
    sx = parameters[count]; count+=1
    sy  = parameters[count]; count+=1
    sz  = parameters[count]; count+=1
    inner_radius = parameters[count]; count+=1
    extra_args = [inner_radius]

    hxy = parameters[count]; count+=1
    hxz = parameters[count]; count+=1
    dx  = 0.0

    hyx = parameters[count]; count+=1
    hyz = parameters[count]; count+=1
    dy  = 0.0

    hzx = parameters[count]; count+=1
    hzy = parameters[count]; count+=1
    dz  = 0.0

    if sx<=0.0 or sy<=0.0 or sz<=0.0 or inner_radius<=0.0 or inner_radius>=1.0:
        return numpy.inf

    affine_transformation_matrix = affine.buildAffineMatrix(sx,hxy,hxz,dx, hyx,sy,hyz,dy, hzx,hzy,sz,dz)
    affine.fillShapeByReverseAffineTransformation(Npoints, affine_transformation_matrix, shape_equation, coordinate_matrix, extra_args)

    start = time.time()
    GV.GV(coordinate_matrix, bs, I_model, Npoints, Nq, qmax, NGPU)
    end = time.time()
    print 'time used: ',end-start
    Iq_model /= Iq_model[0]
    #plt.semilogy(Q,Iq_model)
    #print 'Iq_model: \n',Iq_model

    X2 = 0
    diff=(Iq_model-Iq_exp)/Err_exp
    diff2 = diff*diff
    X2=numpy.sum(diff2)/(Nq-1) 
    print 'parameters: '
    print '%12.8f %12.8f %12.8f'%(sx,hxy,hxz)
    print '%12.8f %12.8f %12.8f'%(hyx,sy,hyz)
    print '%12.8f %12.8f %12.8f'%(hzx,hzy,sz)
    print 'extra_parameters: ',extra_args
    print 'reduced X2: ',X2

    return X2

 

if __name__ == '__main__':

    Nq = 50
    qmax = 0.3

    '''
    # Experimental GV calculation on GPU
    o = sasmol.SasMol(0)
    o.read_pdb(os.path.join('data',os.getcwd().split('/')[-1]+'.pdb'))
    Q = numpy.linspace(0,qmax,Nq)
    I_exp = numpy.zeros(Nq)
    B={'H': -0.37390, 'D': 0.6671, 'C': 0.6646, 'S': 0.2847, 'P': 0.513, 'N': 0.936, 'O': 0.5803}
    bs = map(lambda x: B[x], [item[0] for item in o.element()])
    GV.GV(o._coor[0].T, bs, I_exp, o.natoms(), Nq, qmax, NGPU)
    I_exp /= I_exp[0]
    Err_exp =I_exp/100.
    #plt.plot(Q,I_exp,label='exp')
    #plt.legend(loc='best')
    #plt.show()
    print 'I exp:\n',I_exp
    fileIO.writeIq('output/Iq_exp.txt',Q,I_exp)
    '''

    #get the exp data
    Q,I_exp,Err_exp = fileIO.readIqExp(glob.glob('output/Iq_exp_txt*.txt')[0])
    print 'Exp data:'
    for q,i,e in zip(Q,I_exp,Err_exp):
        print q,i,e

    # run minimization 
    Npoints = 10000
    parameters = (40.0, 40.0, 30.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    bs = numpy.ones(Npoints)
    I_model = numpy.zeros(Nq)
    coordinate_matrix = numpy.matrix(numpy.zeros((3,Npoints)))
    calcX2(parameters, Nq, qmax, Q, I_exp, I_model, Err_exp, bs, coordinate_matrix)
    fileIO.make_and_write_pdb(coordinate_matrix, 'output', 'model_start.pdb')
    fileIO.writeIq('output/Iq_model_start.txt',Q,I_model)

    results = optimize.minimize(calcX2, parameters, args=(Nq, qmax, Q, I_exp, I_model, Err_exp, bs, coordinate_matrix),method='Powell')
    #results = optimize.minimize(calcX2, parameters, args=(Nq, qmax, Q, I_exp, I_model, Err_exp, bs, coordinate_matrix),method='Nelder-Mead')
    #results = optimize.minimize(calcX2, parameters, args=(Nq, qmax, Q, I_exp, I_model, Err_exp, bs, coordinate_matrix),method='BFGS')
    print(results)

    fileIO.make_and_write_pdb(coordinate_matrix, 'output', 'model_final.pdb')
    fileIO.writeIq('output/Iq_model_final.txt',Q,I_model)


