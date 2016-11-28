'''
This file contains file IO handling
Reference: ...
Hailiang Zhang
Jan 2015
'''

import os
import numpy
import math
import locale 
import sassie.sasmol.sasmol as sasmol

def make_and_write_pdb(coordinate_matrix, output_path='./output', output_file_name='output.pdb'):
    '''
    method to write a pdb file via making a sasmol object
    '''
    m1 = sasmol.SasMol(0)
    natoms = coordinate_matrix.shape[1]
    m1.setNatoms = natoms

    atom=[] ; index=[] ; name=[] ; loc=[] ; resname=[] ; chain=[] ; resid=[] ; rescode=[]
    x=[] ; y=[] ; z=[]
    occupancy=[] ; beta=[] ; segname=[] ; element=[] ; charge=[] ; moltype=[]

    coor = numpy.zeros((1,natoms,3),numpy.float32)

    resid_count = -999
    for i in xrange(natoms):
        atom.append('ATOM')
        if(i<99998):
            index.append(i+1)
        else:
            index.append(99999)

        name.append('CA')
        loc.append(' ')
        resname.append('GLY')
        chain.append(' ')
        if(resid_count<9999):
            resid.append(resid_count)
            resid_count+=1
        else:
            resid.append(9999)
        rescode.append(' ')
        x.append(coordinate_matrix[0,i])    
        y.append(coordinate_matrix[1,i])
        z.append(coordinate_matrix[2,i])    
        occupancy.append("  0.00")    
        beta.append("  0.00")    
        element.append('C')
        segname.append('TOR')
        charge.append(' ')
        moltype.append('protein')

    m1.setAtom(atom)
    m1.setIndex(numpy.array(index,numpy.int))
    m1.setName(name)
    m1.setLoc(loc)
    m1.setResname(resname)
    m1.setChain(chain)
    m1.setResid(numpy.array(resid,numpy.int))
    m1.setRescode(rescode)
    m1.setOccupancy(occupancy)
    m1.setBeta(beta)
    m1.setElement(element)
    m1.setSegname(segname)
    m1.setCharge(charge)
    m1.setMoltype(moltype)

    x=numpy.array(x,numpy.float32)
    y=numpy.array(y,numpy.float32)
    z=numpy.array(z,numpy.float32)

    coor[0,:,0] = x ; coor[0,:,1] = y ; coor[0,:,2] = z

    m1.setCoor(coor)

    m1.write_pdb(os.path.join(output_path,output_file_name),0,'w')

    return

def getIqCryson(cryson_Iq_file):
    '''
    method to get the experimental Iq and error values
    NOTE: this code is hacked to read the cryson output from an all-atom calculation as the experimental Iq, and the error is 1% of the intensity
    return: Iq and Err
    '''
    I=[]
    for line in open(cryson_Iq_file).readlines()[1:]:
        words=line.split()
        I.append(locale.atof(words[2]))
    return numpy.array(I)

def readIqExp(ifile):
    '''
    method to write Q and I data
    '''
    Q,I,E=[],[],[]
    for line in open(ifile).readlines():
        words = line.split()
        Q.append(locale.atof(words[0]))
        I.append(locale.atof(words[1]))
        E.append(locale.atof(words[2]))
    return numpy.array(Q), numpy.array(I), numpy.array(E)

def writeIq(ofile, Q, I):
    '''
    method to write Q and I data
    '''
    fout = open(ofile,'w')
    for q,i in zip(Q,I):
        fout.write("%f %f\n"%(q,i))
