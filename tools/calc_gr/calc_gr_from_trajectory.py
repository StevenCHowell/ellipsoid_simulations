import os
import sys
import string
import locale
import numpy
import math
import random
import matplotlib.pyplot as plt
import sasmol.sasmol as sasmol

sys.path.append('./')
import gr as fortran_gr

def get_box_length_list(xst_file_name, stride, sigma=1):

    boxlength_list = []

    boxlength_file = open(xst_file_name, 'r').readlines()
    number_of_lines = len(boxlength_file)
    print 'number_of_lines = ', number_of_lines

    for line in boxlength_file:
        this_line = string.split(line)
        boxlength_list.append(locale.atof(this_line[1])*sigma)

    return boxlength_list

def update_gr(xcoor, ycoor, zcoor, box_length, nbins, natoms2, deltag):

    tgr = numpy.zeros(nbins, numpy.float)

    for i in xrange(natoms2 - 1):
        for j in xrange(i + 1, natoms2):

            xr = xcoor[i] - xcoor[j]
            yr = ycoor[i] - ycoor[j]
            zr = zcoor[i] - zcoor[j]

            xr = xr - box_length * ((xr / box_length).round())
            yr = yr - box_length * ((yr / box_length).round())
            zr = zr - box_length * ((zr / box_length).round())

            r = math.sqrt((xr * xr) + (yr * yr) + (zr * zr))

            if (r < box_length / 2.0):
                ig = int(r / deltag)
                tgr[ig] = tgr[ig] + 2

    return tgr


def main(pdb_file_name, dcd_file_name, xst_file_name, stride, sigma, nskip=0,
         show=False):

    m1 = sasmol.SasMol(0)
    m1.read_pdb(pdb_file_name)

    if(dcd_file_name != False):
        dcdfile = m1.open_dcd_read(dcd_file_name)
        number_of_frames = dcdfile[2]
    else:
        number_of_frames = 1

    print '> found ', number_of_frames, ' frames in dcd file'

    if(xst_file_name != False):
        # box_length_list = get_box_length_list(xst_file_name, stride, sigma)
        box_length = numpy.loadtxt(xst_file_name)

    #box_length_list = [box_length]

    assert len(box_length) == number_of_frames + 1, (
        '# dcd should be one less than the lengeth of box length: ',
        'len(box_length_list), number_of_frames = {}, {}'.format(
            len(box_length), number_of_frames))

    print 'box_length = ', box_length[0]
    print 'box_length = ', box_length[-1]

    # # not sure what this was for
    # min_box_length = min(box_length)
    # max_box_length = max(box_length)

    # print 'mininum box length = ', min_box_length
    # print 'maximum box length = ', max_box_length

    # print '2 PI / boxlength = ', 2.0 * numpy.pi / min_box_length

    box_length_sum = 0.0

    natoms = m1.natoms()

    print 'number of atoms = ', natoms

    nbins = 1000
    deltag = max_box_length / (2.0 * nbins)

    print '> number of bins = ', nbins
    print '> delta g(r) = ', deltag

    sum_gr = numpy.zeros(nbins, numpy.float)
    fsum_gr = numpy.zeros(nbins, numpy.float)
    gr = numpy.zeros(nbins, numpy.float)
    fgr = numpy.zeros(nbins, numpy.float)
    r = numpy.zeros(nbins, numpy.float)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_ylabel('g(r)')
    ax1.set_xlabel('r')

    count = 0


    for i in xrange(nskip, number_of_frames):
        print i + 1,
        sys.stdout.flush()
        if(dcd_file_name != False):
            m1.read_dcd_step(dcdfile, i)
        xcoor = m1.coor()[0, :, 0]*sigma
        ycoor = m1.coor()[0, :, 1]*sigma
        zcoor = m1.coor()[0, :, 2]*sigma

        ftgr = numpy.zeros(nbins, numpy.float)
        ftgr = fortran_gr.update_gr(
            xcoor, ycoor, zcoor, box_length[i], nbins, deltag)
        fsum_gr = fsum_gr + ftgr

        box_length_sum += box_length[i]
        count += 1

    average_box_length = box_length_sum / count

    rho = natoms / (average_box_length**3.0)

    for i in xrange(nbins):
        r[i] = deltag * (i + 0.5)
        vb = (((i + 1)**3) - (i**3)) * (deltag**3)
        nid = (4.0 / 3.0) * numpy.pi * vb * rho
        fgr[i] = fsum_gr[i] / (count * natoms * nid)

    line, = ax1.plot(r, fgr, color='red', lw=3)

    plt.savefig('test_500to1000_by' + str(stride) + '.png')
    if show:
        plt.show()
    else:
        plt.close('all')

    with open('test' + str(stride) + '.dat', 'w') as outfile:
        for i in range(len(r)):
            outfile.write('%f\t%f\n' % (r[i], fgr[i]))


if __name__ == '__main__':
    import os.path as op
    sigma = 3.405

    run_path = '../../simulations/lj_sphere_monomer/runs/p_0p14/output'
    pdb_file_name = 'run1.pdb'
    dcd_file_name = 'run1_mod.dcd'
    xst_file_name = 'box_length.txt'
    pdb_file_name = op.join(run_path, pdb_file_name)
    dcd_file_name = op.join(run_path, dcd_file_name)
    xst_file_name = op.join(run_path, xst_file_name)

    stride = 1

    main(pdb_file_name, dcd_file_name, xst_file_name, stride, sigma,
         nskip=1000)
