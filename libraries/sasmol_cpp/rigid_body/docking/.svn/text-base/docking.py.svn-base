'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import sys
import os
import string
import numpy
import time 	
import sassie.sasmol.sasmol as sasmol
import sassie.simulate.constraints.constraints as constraints
from sassie.simulate.rigid_body.docking.construct import construct_lmn
from sassie.simulate.rigid_body.docking.calc_threadholds import calc_threadholds

from collections import deque as Queue
import pprint,copy

#	DOCKING
#
#	ADAPTED FROM GRID(IASE) & TWO-BODY GRID
#
#	12/15/09	--	initial coding 				:	jc	
#	12/16/09	--	hard coded for (integrase)		:	jc
#	08/19/11	--	connected to gui				:	jc
#	11/21/14	--	removed two-body grid           :   jc
#
#LC	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''

		DOCKING is the program to move on molecule on a grid
	    and score the interactions by surface complementarity

		Molecule 1 is the reference molecule.

		Molecule 2 is the molecule to be moved on the grid 
			
		The molecules must have "residue" fields
	
		INPUT:
	
	 	filename1:	Is the PDB file for molecule 1 
	
	 	filename2:	Is the PDB file for molecule 2 (the molecule to be moved)

	 	outfile:	Name of file with the new coordinates for molecule 2
	
	 	txtOutput:	Object to send textual information back to the GUI


		OUTPUT:
	
	 	outfile:	A file with the aligned coordinates
		cartoutfile:	A file with x,y,z,thetax,thetay,thetaz coords of molecule 2

	 	txtOutput: Object text sent to GUI
	
'''

def print_failure(message,txtOutput):

	txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
	txtOutput.put(">>>> RUN FAILURE <<<<\n")
	txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
	txtOutput.put(message)
	
	return

def report_status(this_trial,total_number_of_trials,txtOutput):

                #if(((i+1)%(float(trials)/100.0)==0 or (trials<10))):
	fraction_done = (float(this_trial)/float(total_number_of_trials))
#	progress_string='\nCOMPLETED '+str(i+1)+' of '+str(trials)+' : '+str(fraction_done*100.0)+' % done'
#	print('%s\n' % progress_string)
#	print accepted,' configurations accepted out of ',nsteps,(float(accepted)/nsteps)*100.0,' %\n\n'
	report_string='STATUS\t'+str(fraction_done)
	txtOutput.put(report_string)

	return

def unpack_variables(variables):

	runname		=	variables['runname'][0]
	path 			=	variables['path'][0]
	pdbmol1		=	variables['pdbmol1'][0]	
	pdbmol2		=	variables['pdbmol2'][0]
	ofile			=	variables['ofile'][0]
	accpos		=	variables['accpos'][0]	
	pos			=	variables['pos'][0]		
	trans			=	variables['trans'][0]	
	dtrans		=	variables['dtrans'][0]	
	theta			=	variables['theta'][0]	
	dtheta		=	variables['dtheta'][0]	
	basis			=	variables['basis'][0]	
	cutoff		=	variables['cutoff'][0]	
	lowrg			=	variables['lowrg'][0]	
	highrg		=	variables['highrg'][0]	
	zflag			=	variables['zflag'][0]	
	zcutoff		=	variables['zcutoff'][0]	
	cflag			=	variables['cflag'][0]	
	confile		=	variables['confile'][0]	
	nexsegments1	=	variables['nexsegments1'][0]
	nsegments1		=	variables['nsegments1'][0]	
	reslow1		=	variables['reslow1'][0]	
	numcont1		=	variables['numcont1'][0]	
	nexsegments2	=	variables['nexsegments2'][0]
	nsegments2		=	variables['nsegments2'][0]	
	reslow2		=	variables['reslow2'][0]	
	numcont2		=	variables['numcont2'][0]	

	return runname,path,pdbmol1,pdbmol2,ofile,accpos,pos,trans,dtrans,theta,dtheta,basis,cutoff,lowrg,highrg,zflag,zcutoff,cflag,confile,nexsegments1,nsegments1,reslow1,numcont1,nexsegments2,nsegments2,reslow2,numcont2


def euler_rotation(m,phi,theta,psi):
   c11 = numpy.cos(theta)*numpy.cos(psi)
   c12 = numpy.cos(phi)*numpy.sin(psi) + numpy.sin(phi)*numpy.sin(theta)*numpy.cos(psi)
   c13 = numpy.sin(phi)*numpy.sin(psi) - numpy.cos(phi)*numpy.sin(theta)*numpy.cos(psi)
   c21 = -numpy.cos(theta)*numpy.sin(psi)
   c22 = numpy.cos(phi)*numpy.cos(psi) - numpy.sin(phi)*numpy.sin(theta)*numpy.sin(psi)
   c23 = numpy.sin(phi)*numpy.cos(psi) + numpy.cos(phi)*numpy.sin(theta)*numpy.sin(psi)
   c31 = numpy.sin(theta)
   c32 = -numpy.sin(phi)*numpy.cos(theta)
   c33 = numpy.cos(phi)*numpy.cos(theta)
   C = numpy.matrix([[c11,c12,c13],[c21,c22,c23],[c31,c32,c33]])
   newcoor = numpy.array(m.coor()[0]*C)
   m.setCoor(numpy.array([newcoor]))


def reo_slide(slide, N):
	for i in range(N):
		for j in range(N/2):
			tmp = slide[i][j+N/2]
			slide[i][j+N/2]=slide[i][j]
			slide[i][j]=tmp
	for i in range(N):
		for j in range(N/2):
			tmp = slide[j+N/2][i]
			slide[j+N/2][i]=slide[j][i]
			slide[j][i]=tmp



def merge_Q(Q, Qtmp, mMax):
	Q_new = numpy.append(Q,Qtmp)
	Q_new = numpy.sort(Q_new, order=['value'], kind='quicksort')
	if len(Q_new)<mMax:
		return Q_new
	else:
		return Q_new[-mMax:]

def fft_docking(m1,m2,m3,ofile,genpaths,accpos,pos,trans,dtrans,theta,dtheta,cutoff,basis,mask_array1,mask_array2,zflag,zcutoff,cflag,mask_a_array,mask_b_array,distance_array,type_array,txtOutput):
	nMax = 10
	Q_coarse=[]
	Q_fine=[]
	Q_coarse=numpy.array([(0, 0, 0, 0, 0.0)]*nMax,dtype=[('iphi',numpy.int8),('itheta',numpy.int8),('ipsi',numpy.int8),('index',int),('value',float)])
	Q_fine=numpy.array([(0, 0, 0, 0, 0.0)]*nMax,dtype=[('iphi',numpy.int8),('itheta',numpy.int8),('ipsi',numpy.int8),('index',int),('value',float)])

	print '\n\n=======================================\nFFT Docking\n\n'

	# Center and properly orient the acceptor
	m1.center(0)
	m2.center(0)

	# Set the box dimension (based on donor shape and size)
	
	# Assign parameters
	r = 1.8 # parameter to derive almn and blmn
	d = 2.5 # thickness of the surface layer
	#r = 2.4 # parameter to derive almn and blmn
	#d = 1.0 # thickness of the surface layer
	rou = -15 # molecular 1 interior parameter
	delta = 1 # molecular 2 interior parameter
	angle_step = 5*numpy.pi/180.0
	eta = 1.0 # grid step size, 1.0-1.2
	N = 90 # Number of grid points
	eta_fine = 0.7 # grid step size, 0.7-0.8
	N_fine = 128


	# Construct almn,AopqC
	almn = construct_lmn(m1.coor()[0], rou, N, r, d, eta, 0)
	AopqC = numpy.fft.fftn(almn,(N,N,N)).conj()
	almn_fine = construct_lmn(m1.coor()[0], rou, N_fine, r, d, eta_fine, 0)
	AopqC_fine = numpy.fft.fftn(almn_fine,(N_fine,N_fine,N_fine)).conj()

	# Scan Stage
	count=0
	threadholds = [0.005,0.01,0.02,0.05,0.1, 0.2, 0.5, 0.95]
	tot_nf_threadholds = [0.0]*len(threadholds)
	tot_nf = 0
	import datetime
	for phi in numpy.arange(0,numpy.pi*2.0,angle_step):
		print "PHI %6.2f"%phi
		start =datetime.datetime.now()
		for theta in numpy.arange(0,numpy.pi*1.0,angle_step):
			#print 'THETA %6.2f'%theta
			for psi in numpy.arange(0,numpy.pi*1.0,angle_step):
				#print 'PSI %6.2f'%psi
				m2tmp = copy.deepcopy(m2)
				euler_rotation(m2tmp,phi,theta,psi)
				blmn = construct_lmn(m2tmp.coor()[0], delta, N, r, d, eta,0)
				Bopq = numpy.fft.fftn(blmn,(N,N,N))
				Copq = AopqC*Bopq
				clmn = numpy.fft.ifftn(Copq,(N,N,N))
				
				clmn_tmp = clmn.real
				clmn_tmp = clmn_tmp[clmn_tmp>0.01]
				clmn_tmp = numpy.sort(clmn_tmp)[::-1]
				tot_nf += len(clmn_tmp)
				nf_threadholds = calc_threadholds(clmn_tmp,threadholds)
				tot_nf_threadholds += nf_threadholds
				
				count+=N**3
		end =datetime.datetime.now()
		print 'time used ',(end-start).seconds

	for i in range(len(threadholds)):
		print 'Number of frames with occupancy > ','%6.1f'%(threadholds[i]*100.0),'% is: ','%6d'%tot_nf_threadholds[i],' (','%6.2f'%(tot_nf_threadholds[i]/float(tot_nf)*100),'% out of ',tot_nf,' of frames with positive correlation)'
	print sum(tot_nf_threadholds), ' out of total number of positive frames ',tot_nf

def simple_docking(variables,txtOutput):

	runname,path,pdbmol1,pdbmol2,ofile,accpos,pos,trans,dtrans,theta,dtheta,basis,cutoff,lowrg,highrg,zflag,zcutoff,cflag,confile,nexsegments1,nsegments1,reslow1,numcont1,nexsegments2,nsegments2,reslow2,numcont2 = unpack_variables(variables)

	if(runname[-1]=='/'):
		lin=len(runname)
		runname=runname[:lin-1]

	direxist=os.path.exists(runname)
	if(direxist==0):
		os.system('mkdir -p '+runname)

	genpath=runname+'/two_body_grid'
	genpaths=genpath+'/'
	direxist=os.path.exists(genpath)
	if(direxist==0):
		os.system('mkdir -p '+genpath)

	m1=sasmol.SasMol(0)
	m2=sasmol.SasMol(1)
	m3=sasmol.SasMol(2)

	m1.read_pdb(path+'/'+pdbmol1)
	m2.read_pdb(path+'/'+pdbmol2)

	error = m3.merge_two_molecules(m1,m2)

	if(error!=[]):
		print 'ERROR:'+error[0]
		print 'ERROR:'+error[0]
		print 'ERROR:'+error[0]

	m3.write_pdb(genpaths+ofile+'.pdb',0,'w')

	cpst = 'cp '+path+'/'+pdbmol1+' '+genpaths
	os.system(cpst)
	cpst = 'cp '+path+'/'+pdbmol2+' '+genpaths
	os.system(cpst)

	frame=0

	mm1=m1.calcminmax()
	mm2=m2.calcminmax()

	'''	
	print 'mm1 = ',mm1
	print 'mm2 = ',mm2
	'''

	# set overlap basis for each molecule

	segment_names_1 = string.split(nsegments1,',')
	segment_names_2 = string.split(nsegments2,',')

	'''
	print 'segment_names_1 = ',segment_names_1
	print 'segment_names_2 = ',segment_names_2
	'''

	if(nexsegments1 > 0):

		for i in xrange(nexsegments1):
			if(i==0):	
				basis1st = "(name[i] == 'CA' and not (segname[i] == '"+segment_names_1[i]+"' and ( resid[i] >= "+str(reslow1[i])+" and resid[i] <= "+str(reslow1[i]+numcont1[i])+")))"
			else:
				basis1st = basis1st + " or (name[i] == 'CA' and not (segname[i] == '"+segment_names_1[i]+"' and ( resid[i] >= "+str(reslow1[i])+" and resid[i] <= "+str(reslow1[i]+numcont1[i])+" )))"

	else:
		basis1st = "name[i] == 'CA'"
	
	'''
	print 'basis1st = ',basis1st
	'''

	if(nexsegments2 > 0):

		for i in xrange(nexsegments2):
			if(i==0):	
				basis2st = "(name[i] == 'CA' and not (segname[i] == '"+segment_names_2[i]+"' and ( resid[i] >= "+str(reslow2[i])+" and resid[i] <= "+str(reslow2[i]+numcont2[i])+")))"
			else:
				basis2st = basis2st + " or (name[i] == 'CA' and not (segname[i] == '"+segment_names_2[i]+"' and ( resid[i] >= "+str(reslow2[i])+" and resid[i] <= "+str(reslow2[i]+numcont2[i])+" )))"

	else:
		basis2st = "name[i] == 'CA'"
	
	'''
	print 'basis2st = ',basis2st
	'''

	error,mask_array1 = m1.get_subset_mask(basis1st)
	error,mask_array2 = m2.get_subset_mask(basis2st)

	if(cflag == 1):
		filter_flag = 0
		error,constraint_basis1_array, constraint_basis2_array,distance_array,type_array = constraints.read_constraints(m3,confile,filter_flag)

		mask_a_array = [] ; mask_b_array = []

		for i in xrange(len(distance_array)):
			print constraint_basis1_array[i]
			print constraint_basis2_array[i]
			print distance_array[i]
			print type_array[i]

			error,local_mask_a_array = m3.get_subset_mask(constraint_basis1_array[i])
			error,local_mask_b_array = m3.get_subset_mask(constraint_basis2_array[i])
	
			mask_a_array.append(local_mask_a_array)	
			mask_b_array.append(local_mask_b_array)	

	else:
		mask_a_array = [] ; mask_b_array = []
		distance_array = [] ; type_array = []

	fft_docking(m1,m2,m3,ofile,genpaths,accpos,pos,trans,dtrans,theta,dtheta,cutoff,basis,mask_array1,mask_array2,zflag,zcutoff,cflag,mask_a_array,mask_b_array,distance_array,type_array,txtOutput)
			
	return

