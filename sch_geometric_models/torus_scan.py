#!/usr/bin/env python
# coding:utf-8
"""
    Author:  Steven C. Howell --<steven.howell@nist.gov>
    Purpose: generating an array of geometric models
    Created: 12/08/2016

00000000011111111112222222222333333333344444444445555555555666666666677777777778
12345678901234567890123456789012345678901234567890123456789012345678901234567890

"""

import errno
import logging
import os
import sys
import time

import numpy as np
import sasmol.sasmol as sasmol

sys.path.append('/home/schowell/data/code/ellipsoid_simulations/hz_affine/')
from affine.shapes import torus_equation as shape_eq
from affine import fileIO

sys.path.append('/home/schowell/data/scratch/sascalc_demo/')
import sascalc

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)


class AffineParameters:

    def __init__(self, sx, sy, sz, hxy=0.0, hxz=0.0, hyx=0.0, hyz=0.0, hzx=0.0,
                 hzy=0.0, dx=0.0, dy=0.0, dz=0.0):
        self.sx = sx
        self.sy = sy
        self.sz = sz
        self.hxy = hxy
        self.hxz = hxz
        self.hyx = hyx
        self.hyz = hyz
        self.hzx = hzx
        self.hzy = hzy
        self.dx = dx
        self.dy = dy
        self.dz = dz


def fill_shape(full_affine_matrix, shape_eq, coor, extra_args):
    '''
    Purpose:
        Fill a shape in the affine space and do the bondary check
        by reverse affine transformation.

    Return:
        the affine coordinate matrix after transformation
    '''
    n_coor = len(coor[0])
    affine_matrix = full_affine_matrix[:3, :3]

    # get the boundaries
    # ZHL note: only check for 4 cube vertex because the other 4 is just the
    # reflections; this needs change if translation is included in the affine
    # transformation matrix'''
    cube_vertex_coor_matrix = np.matrix([[-1.0, 1.0, -1.0, 1.0],
                                         [-1.0, -1.0, 1.0, 1.0],
                                         [-1.0, -1.0, -1.0, -1.0]])

    affine_vertex_coor_matrix = affine_matrix * cube_vertex_coor_matrix
    affine_vertex_coor_matrix = np.concatenate((affine_vertex_coor_matrix,
                                                -affine_vertex_coor_matrix),
                                               axis=1)
    affine_vertex_coor_min = affine_vertex_coor_matrix.min(axis=1)
    affine_edge_max = (affine_vertex_coor_matrix.max(axis=1) -
                       affine_vertex_coor_min).max()
    # check for points
    reverse_trans_matrix = affine_matrix.I
    n_coor_found = 0
    while n_coor_found < n_coor:
        points_remaining = n_coor - n_coor_found

        # create random points
        trial_coor = np.random.rand(3, points_remaining)
        trial_coor *= affine_edge_max  # rescale the points
        trial_coor += affine_vertex_coor_min  # offset the points

        trial_coor_orth = reverse_trans_matrix * trial_coor

        good_coor = [shape_eq(*tuple(np.array(coor_i).ravel()), extra=extra_args)
                     for coor_i in trial_coor_orth.T]   # this may need to be changed
        good_coor = np.array(good_coor)

        n_new_coor_found = n_coor_found + sum(good_coor)
        coor[0, n_coor_found:n_new_coor_found, :] = trial_coor.T[good_coor]
        n_coor_found = n_new_coor_found

        assert n_coor_found, 'ERROR - check inputs:\ {}'.format(locals())


def scan_values(sx_range, sy_range, sz_range, r_small_range, mol, output_dir,
                n_scan=50):

    hxy = hxz = hyx = hyz = hzx = hzy = dx = dy = dz = 0
    sx = sy = sz = 1
    affine_matrix = np.matrix([[sx, hxy, hxz, dx],
                               [hyx, sy, hyz, dy],
                               [hzx, hzy, sz, dz],
                               [0.0, 0.0, 0.0, 1.0]])
    sx0 = min(sx_range)
    sx1 = max(sx_range)
    sy0 = min(sy_range)
    sy1 = max(sy_range)
    sz0 = min(sz_range)
    sz1 = max(sz_range)
    r0 = min(r_small_range)
    r1 = max(r_small_range)

    dcd_fname = 'sx{}to{}_sy{}to{}_sz{}to{}_r{}to{}_by{}.dcd'.format(sx0, sx1,
                                                                     sy0, sy1,
                                                                     sz0, sz1,
                                                                     r0, r1,
                                                                     n_scan)
    dcd_outfile = mol.open_dcd_write(os.path.join(output_dir, dcd_fname))

    tic0 = time.time()
    n = 0
    for sx in np.linspace(sx0, sx1, n_scan):
        affine_matrix[0, 0] = sx
        tic1 = time.time()

        for sy in np.linspace(sy0, sy1, n_scan):
            affine_matrix[1, 1] = sy

            for sz in np.linspace(sz0, sz1, n_scan):
                affine_matrix[2, 2] = sz

                for r_small in np.linspace(r0, r1, n_scan):
                    fill_shape(affine_matrix, shape_eq, mol.coor(),
                               [r_small])
                    n += 1
                    mol.write_dcd_step(dcd_outfile, 0, n)
        toc1 = time.time() - tic1
        logging.debug('time to scan at sx={}: {} s or {} h'.format(sx, toc1,
                                                                   toc1/3600))


    toc0 = time.time() - tic0
    logging.debug('calculation time: {} s or {} h'.format(toc0, toc0/3600))
    mol.close_dcd_write()


def mkdir_p(path):
    '''
    make directory recursively
    adapted from http://stackoverflow.com/questions/600268/
    '''
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


if __name__ == '__main__':

    output_dir = 'run5'
    mkdir_p(output_dir)

    mol = sasmol.SasMol(0)
    try:
        mol.read_pdb(os.path.join(output_dir, 'dummy.pdb'))
    except IOError:
        n_atoms = 10000
        coor = np.zeros((n_atoms, 3), order='F')
        fileIO.make_and_write_pdb(coor.T, output_dir, 'dummy.pdb')
        mol.read_pdb(os.path.join(output_dir, 'dummy.pdb'))

    sx_range = [10, 40]
    sy_range = [10, 40]
    sz_range = [10, 40]
    r_small_range = [0.05, 0.95]
    n_scan = 5
    scan_values(sx_range, sy_range, sz_range, r_small_range, mol, output_dir,
                n_scan)

