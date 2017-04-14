'''
This file contains modules for unit basic shapes
The unit basic shapes are the largest basic shape that can fit into a unit cube that is centered in the origin and has a length of 2
Hailiang Zhang
Jan 2015
'''

import os
import numpy
import math
import random


def cube_equation(x,y,z,**kwargs):
    '''
    Check whether a point are within a unit cube or not
    Return: True if inside or on surface, False if outside
    '''
    if -1.<=x<=1. and -1.<=y<=1. and -1.<=z<=1.:
        return True
    else:
        return False


def sphere_equation(x,y,z,**kwargs):
    '''
    Check whether a point are within a unit sphere or not
    Return: True if inside or on surface, False if outside
    '''
    if not cube_equation(x,y,z):
        return False

    if math.sqrt(x*x+y*y+z*z) <= 1.0:
        return True
    else:
        return False


def prism_equation(x,y,z,**kwargs):
    '''
    Check whether a point are within a unit prism or not
    Return: True if inside or on surface, False if outside
    '''
    if not cube_equation(x,y,z):
        return False

    if y>0.7320508075688772:
        return False

    if abs((y+1.)/x)<math.sqrt(3.):
        return False
    else:
        return True


def cappedPrism_equation(x,y,z,**kwargs):
    '''
    Check whether a point are within a unit capped prism or not
    Return: True if inside or on surface, False if outside
    '''
    if not prism_equation(x,y,z):
        return False

    portion_cap_length = kwargs['extra'][0]
    portion_handle_width = kwargs['extra'][1]

    if (1.-z)/2.0 <= portion_cap_length:
        return True
    else:
        y0 = 2./math.sqrt(3.)-1.0
        if (y-y0) > portion_handle_width/math.sqrt(3.) or (y-y0) < -portion_handle_width*2./math.sqrt(3.) or abs(x) > portion_handle_width:
            return False
        elif abs((y+1.-(1.-portion_handle_width)*2./math.sqrt(3.))/x)<math.sqrt(3.):
            return False
        else:
            return True


def pyramid_equation(x,y,z,**kwargs):
    '''
    Check whether a point are within a unit pyramid or not
    Return: True if inside or on surface, False if outside
    '''
    if not cube_equation(x,y,z):
        return False

    ratio = (1.-z)/2.
    if -ratio<=x<=ratio and -ratio<=y<=ratio and -ratio<=z<=ratio:
        return True
    else:
        return False


def cylinder_equation(x,y,z,**kwargs):
    '''
    Check whether a point are within a unit cylinder or not
    Return: True if inside or on surface, False if outside
    '''
    if not cube_equation(x,y,z):
        return False

    if math.sqrt(x*x+y*y) <= 1.0:
        return True
    else:
        return False


def hollowCylinder_equation(x,y,z,**kwargs):
    '''
    Check whether a point are within a unit hollowClinder or not
    Return: True if inside or on surface, False if outside
    '''
    if not cylinder_equation(x,y,z):
        return False

    inner_radius = kwargs['extra'][0]

    if math.sqrt(x*x+y*y) >= inner_radius:
        return True
    else:
        return False


def shell_equation(x,y,z,**kwargs):
    '''
    Check whether a point are within a unit half_shell or not
    Return: True if inside or on surface, False if outside
    '''
    if not sphere_equation(x,y,z):
        return False

    ri = kwargs['extra'][0]

    if (x*x + y*y) >= ri*ri*(1.-z*z) :
        return True
    else:
        return False


def halfShell_equation(x,y,z,**kwargs):
    '''
    Check whether a point are within a unit half_shell or not
    Return: True if inside or on surface, False if outside
    '''

    ri = kwargs['extra'][0]

    if z<0.0 or not shell_equation(x,y,z,extra=[ri]):
        return False
    else:
        return True


def torus_equation(x,y,z,**kwargs):
    '''
    Check whether a point are within a unit half_shell or not
    Return: True if inside or on surface, False if outside
    '''
    r_small = kwargs['extra'][0]

    if abs(z)>r_small or not sphere_equation(x,y,z):
        return False

    ratio = (1.-r_small) / math.sqrt(x*x + y*y)
    ring_orig_x = x*ratio
    ring_orig_y = y*ratio
    if math.sqrt((x-ring_orig_x)*(x-ring_orig_x) + (y-ring_orig_y)*(y-ring_orig_y) + z*z) <= r_small:
        return True
    else:
        return False


"""
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

