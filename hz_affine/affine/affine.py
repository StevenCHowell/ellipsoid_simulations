'''
This file contains modules for affine-transformation
Reference: ...
Hailiang Zhang
Jan 2015
'''

import numpy

def buildAffineMatrix(sx,hxy,hxz,dx, hyx,sy,hyz,dy, hzx,hzy,sz,dz):
    '''
    Build the affine matrix
    Return: affine matrix
    '''

    affine_matrix = numpy.matrix([[sx,hxy,hxz,dx], [hyx,sy,hyz, dy], [hzx,hzy,sz,dz], [0.0,0.0,0.0,1.0]])
    return affine_matrix

def padCoordinateMatrix(coordinate_matrix):
    '''
    Pad a coorindate matrix with a row of 1's so that it can be consumed by the affine transformation
    Return: the coordinate matrix which has been padded with a row of 1's
    '''
    dummy = numpy.ones((1,coordinate_matrix.shape[1]))
    return numpy.concatenate((coordinate_matrix,dummy),axis=0)

def unPadCoordinateMatrix(affine_coordinate_matrix):
    '''
    Un-pad a coorindate matrix with a row of 1's so that it contains only the coordinates information
    Return: the coordinate matrix which has been un-padded
    '''
    return affine_coordinate_matrix[0:3,:]

def affineTransformation(affine_transformation_matrix, affine_coordinate_matrix):
    '''
    Perform the affine transformation
    Return: the affine coordinate matrix after transformation
    '''
    result = affine_transformation_matrix*affine_coordinate_matrix
    return result

def fillShapeByReverseAffineTransformation(Npoints, full_affine_transformation_matrix, shape_equation, coordinate_matrix, extra_args):
    '''
    Fill a shape in the affine space
        and do the bondary check by reverse affine transformation
    Return: the affine coordinate matrix after transformation
    '''
    affine_transformation_matrix = full_affine_transformation_matrix[:3,:3]
    # get the boundaries
    # ZHL note: only check for 4 cube vertex because the other 4 is just the reflections; this needs change if translation is included in the affine transformation matrix'''
    cube_vertex_coordinate_matrix = numpy.matrix([[-1.0,1.0,-1.0,1.0],[-1.0,-1.0,1.0,1.0],[-1.0,-1.0,-1.0,-1.0]])
    affine_vertex_coordinate_matrix = affine_transformation_matrix * cube_vertex_coordinate_matrix
    affine_vertex_coordinate_matrix = numpy.concatenate((affine_vertex_coordinate_matrix,-affine_vertex_coordinate_matrix),axis=1)
    affine_vertex_coordinate_min = affine_vertex_coordinate_matrix.min(axis=1)
    affine_edge_max = (affine_vertex_coordinate_matrix.max(axis=1)-affine_vertex_coordinate_min).max()
    # check for points
    reverse_affine_transformation_matrix = affine_transformation_matrix.I
    points_found = 0
    while points_found < Npoints:
        points_remaining = Npoints - points_found
        coordinate_matrix_remaining = numpy.random.rand(3,points_remaining)
        coordinate_matrix_remaining[0,:] = affine_vertex_coordinate_min[0] + coordinate_matrix_remaining[0,:]*affine_edge_max
        coordinate_matrix_remaining[1,:] = affine_vertex_coordinate_min[1] + coordinate_matrix_remaining[1,:]*affine_edge_max
        coordinate_matrix_remaining[2,:] = affine_vertex_coordinate_min[2] + coordinate_matrix_remaining[2,:]*affine_edge_max
        coordinate_orth_matrix_remaining = reverse_affine_transformation_matrix * coordinate_matrix_remaining
        good_points = numpy.array([shape_equation(*tuple(numpy.array(coordinate).ravel()), extra=extra_args) for coordinate in coordinate_orth_matrix_remaining.T])
        new_points_found = points_found + len(filter(lambda x: x, good_points))
        coordinate_matrix[:, points_found:new_points_found] = coordinate_matrix_remaining.T[good_points].T
        points_found = new_points_found
        if points_found == 0:
            import pprint
            pprint.pprint(locals())
            exit(0)

if __name__ == '__main__':
    coordinate_matrix = numpy.matrix([[1.0],[2.0],[3.0]])
    print 'coordinate matrix:\n',coordinate_matrix
    affine_coordinate_matrix = padCoordinateMatrix(coordinate_matrix)
    print 'affine coordinate matrix:\n',affine_coordinate_matrix
    affine_transformation_matrix = buildAffineMatrix(1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0)
    print 'affine transformation matrix:\n',affine_transformation_matrix
    result = affineTransformation(affine_transformation_matrix, affine_coordinate_matrix)
    print 'result:\n',result

