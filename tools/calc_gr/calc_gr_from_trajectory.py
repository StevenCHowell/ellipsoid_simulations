import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import sasmol.sasmol as sasmol

sys.path.append('./')
import gr as fortran_gr


def get_box_length_list(xst_file_name, stride, sigma=1):
    import string
    import locale

    boxlength_list = []

    boxlength_file = open(xst_file_name, 'r').readlines()
    number_of_lines = len(boxlength_file)
    print 'number_of_lines = ', number_of_lines

    for line in boxlength_file:
        this_line = string.split(line)
        boxlength_list.append(locale.atof(this_line[1])*sigma)

    return boxlength_list


def update_gr(xcoor, ycoor, zcoor, box_length, n_bins, natoms2, delta_g):

    tgr = np.zeros(n_bins, np.float)

    for i in xrange(natoms2 - 1):
        for j in xrange(i + 1, natoms2):

            xr = xcoor[i] - xcoor[j]
            yr = ycoor[i] - ycoor[j]
            zr = zcoor[i] - zcoor[j]

            xr = xr - box_length * ((xr / box_length).round())
            yr = yr - box_length * ((yr / box_length).round())
            zr = zr - box_length * ((zr / box_length).round())

            r = np.sqrt((xr * xr) + (yr * yr) + (zr * zr))

            if (r < box_length / 2.0):
                ig = int(r / delta_g)
                tgr[ig] = tgr[ig] + 2

    return tgr


def main(pdb_fname, run_log_fname, stride=1, sigma=1, dcd_fname=None,
         n_skip=0, show=False):

    box_mol = sasmol.SasMol(0)
    box_mol.read_pdb(pdb_fname)

    if dcd_fname:
        dcd_file = box_mol.open_dcd_read(dcd_fname)
        n_frames = dcd_file[2]
    else:
        n_frames = 1

    run_log = np.loadtxt(run_log_fname)

    assert len(run_log) == n_frames + 1, (
        '# dcd should be one less than the lengeth of box length: ',
        'len(box_length_list), number_of_frames = {}, {}'.format(
            len(run_log), n_frames))

    box_length = run_log[:, 1]
    print 'box_length: [ fist last] = {}'.format(box_length[[0, -1]])

    n_atoms = box_mol.natoms()
    n_bins = 1000
    delta_g = box_length / (2.0 * n_bins)
    rho = n_atoms / box_length ** 3.0
    n_frames_used = n_frames - n_skip
    gr_all = np.empty((n_frames_used, n_bins))  # one g(r) for ecah dcd frame
    gr = np.empty((n_bins, 2))
    r = np.empty(n_bins)

    fsum_gr = np.zeros(n_bins, np.float)
    fgr = np.zeros(n_bins, np.float)

    count = 0


    box_length_sum = 0.0

    bin_index = np.arange(n_bins) + 1
    bin_volume = np.outer(delta_g ** 3, (bin_index + 1) ** 3 - bin_index ** 3)

    # for i in xrange(n_skip, n_frames):
    for i in xrange(n_skip, n_skip + 10):
        sys.stdout.flush()

        if dcd_fname:
            box_mol.read_dcd_step(dcd_file, i)

        x_coor = box_mol.coor()[0, :, 0] * sigma
        y_coor = box_mol.coor()[0, :, 1] * sigma
        z_coor = box_mol.coor()[0, :, 2] * sigma

        gr_all[i] = fortran_gr.update_gr(x_coor, y_coor, z_coor, box_length[i],
                                    n_bins, delta_g[i])

        for i in xrange(n_bins):
            nid = (4.0 / 3.0) * np.pi * bin_volume * rho

    gr[:, 0] = r
    gr[:, 1] = np.mean(gr_all, axis=0)

    for i in xrange(n_bins):
        r[i] = delta_g * (i + 0.5)
        fgr[i] = fsum_gr[i] / (count * n_atoms * nid)


    np.savetxt(gr, 'test{}.dat'.format(stride))



def plot_gr(gr):

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_ylabel('g(r)')
    ax1.set_xlabel('r')
    line, = ax1.plot(r, fgr, color='red', lw=3)

    plt.savefig('test_500to1000_by' + str(stride) + '.png')

    if show:
        plt.show()
    else:
        plt.close('all')


if __name__ == '__main__':
    import os.path as op
    sigma = 3.405

    run_path = '../../simulations/lj_sphere_monomer/runs/p_0p14/output'
    pdb_fname = 'run1.pdb'
    dcd_file_name = 'run1_mod.dcd'
    xst_file_name = 'box_length.txt'
    pdb_fname = op.join(run_path, pdb_fname)
    dcd_file_name = op.join(run_path, dcd_file_name)
    xst_file_name = op.join(run_path, xst_file_name)

    gr = main(pdb_fname, xst_file_name, sigma=sigma, dcd_fname=dcd_file_name,
              n_skip=1000)

    plot_gr(gr)