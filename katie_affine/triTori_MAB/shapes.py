'''
This file contains modules for basic shapes
Reference: ...
Hailiang Zhang
Jan 2015
'''

import os
import numpy
import math
import random
import sassie.sasmol.sasmol as sasmol

def cube_equation(x,y,z,half_length=1.0):
    '''
    Check whether a point are within a unit cube or not
    Return: True if inside or on surface, False if outside
    '''
    hl = half_length
    if x>=-hl and x<=hl and y>=-hl and y<=hl and z>=-hl and z<=hl:
        return True
    else:
        return False

def cylinder_equation(x,y,z,radius=1.0,half_length=1.0):
    '''
    Check whether a point are within a unit cylinder or not
    Return: True if inside or on surface, False if outside
    '''
    if math.sqrt(x*x+y*y) <= radius and z>=-half_length and z<=half_length:
        return True
    else:
        return False

def sphere_equation(x,y,z,radius=1.0):
    '''
    Check whether a point are within a unit sphere or not
    Return: True if inside or on surface, False if outside
    '''
    if math.sqrt(x*x+y*y+z*z) <= radius:
        return True
    else:
        return False

def halfshell_equation(x,y,z,ro=1.0,ri=0.4,h=1.0): #ZHL note: affine transformation doesn't alter ro/ri; ri is hardwired for GROEL
    '''
    Check whether a point are within a unit half_shell or not
    Return: True if inside or on surface, False if outside
    '''
    if z<0 or z>h or x<-ro or x>ro or y<-ro or y>ro:
        return False
    ro2_slice = ro*ro - z*z
    ri2_slice = ri*ri - z*z
    x2y2 = x*x + y*y
    if x2y2 >= ri2_slice and x2y2 <= ro2_slice:
        return True
    else:
        return False

def fillSphere(Npoints,radius=1.0):
    '''
    Fill a sphere with points
    Return: a numpy array of 3*Npoints where Npoints is the number of points filled
    '''
    coordinates = numpy.zeros((3,Npoints))
    count = 0
    N1 = int(math.pow(Npoints,1./3.))
    for x in numpy.linspace(-radius,radius,N1):
        for y in numpy.linspace(-radius,radius,N1):
            for z in numpy.linspace(-radius,radius,N1):
                if math.sqrt(x*x+y*y+z*z) <= radius:
                    coordinates[:,count]=[x,y,z]
                    count += 1
    return coordinates
"""

def fillSphere(Npoints,radius=1.0):
    '''
    Fill a sphere with points
    Return: a numpy array of 3*Npoints where Npoints is the number of points filled
    '''
    coordinates = numpy.zeros((3,Npoints))
    count = 0
    while count < Npoints:
        x = random.uniform(-radius,radius)
        y = random.uniform(-radius,radius)
        z = random.uniform(-radius,radius)
        if math.sqrt(x*x+y*y+z*z) <= radius:
            coordinates[:,count]=[x,y,z]
            count += 1
            #print 'Filled %d out of %d points'%(count,Npoints)
    return coordinates
"""


if __name__ == '__main__':

    pass

