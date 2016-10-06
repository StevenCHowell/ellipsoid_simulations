import logging
import os
import sys
import numpy as np
import sasmol.sasmol as sasmol

sys.path.append('./')
import gr as fortran_gr

FORMAT = "%(asctime)-15s: %(message)s"
logging.basicConfig(format=FORMAT, level=logging.DEBUG)


def update_gr(xcoor, ycoor, zcoor, box_length, n_bins, n_atoms, delta_r):

    gr = np.zeros(n_bins, np.float)

    for i in xrange(n_atoms - 1):
        for j in xrange(i + 1, n_atoms):

            dx = xcoor[i] - xcoor[j]
            dy = ycoor[i] - ycoor[j]
            dz = zcoor[i] - zcoor[j]

            dx -= box_length * ((dx / box_length).round())
            dy -= box_length * ((dy / box_length).round())
            dz -= box_length * ((dz / box_length).round())

            r = np.sqrt((dx * dx) + (dy * dy) + (dz * dz))

            if (r < box_length / 2.0):
                ig = int(r / delta_r)  # round down
                gr[ig] += 2

    return gr


def main(pdb_fname, run_log_fname, stride=1, sigma=1, dcd_fname=None,
         n_skip=0, n_bins=333):

    box_mol = sasmol.SasMol(0)
    box_mol.read_pdb(pdb_fname)

    if dcd_fname:
        dcd_file = box_mol.open_dcd_read(dcd_fname)
        n_frames = dcd_file[2]
    else:
        n_frames = 1

    run_log = np.loadtxt(run_log_fname)

    if len(run_log) != n_frames:
        if len(run_log) == n_frames + 1:
            run_log = run_log[:-1]
            logging.warning('dcd file  had one more frame than the log file, '
                         'discarding last line\ndcd_fname:\t{}\n'
                         'run_log_fname:\t{}'.format(dcd_fname, run_log_fname))
        else:
            logging.error('mismatch between dcd and log file \n'
                          'dcd_fname:\t{}\nrun_log_fname:\t{}'.format(
                              dcd_fname, run_log_fname))

    n_atoms = box_mol.natoms()
    n_gr = n_frames - n_skip  # number of frames, or g(r) curves, to averages

    box_length = run_log[n_skip:, 1] * sigma
    print('box_length: (min, max) = ({}, {})'.format(box_length.min(),
                                                     box_length.max()))

    gr_all = np.zeros((n_gr, n_bins))  # one g(r) for ecah dcd frame
    gr = np.zeros((n_bins, 2))

    # # using a different r_grid for each frame
    # delta_r = box_length / (2.0 * n_bins)
    # bin_index = np.arange(n_bins) + 1
    # rho = (n_atoms / box_length ** 3.0).reshape(-1, 1)  # particle density

    # bin_volume = np.outer(delta_r ** 3, (bin_index + 1) ** 3 - bin_index ** 3)
    # r = np.outer(delta_r, (bin_index + 0.5))
    # n_ideal = (4.0 / 3.0) * np.pi * bin_volume * rho  # expected n for ideal gas

    # using the same r_grid for each frame
    delta_r = box_length.max() / (2.0 * n_bins)  # delg in F&S
    bin_index = np.arange(n_bins) + 1
    rho = (n_atoms / box_length ** 3.0).reshape(-1, 1)  # frame specific density

    # bin_volume = ((bin_index + 1) ** 3 - bin_index ** 3) * delta_r ** 3
    r = (bin_index - 0.5) * delta_r
    bin_volume = 4.0 / 3.0 * np.pi * ((r + delta_r / 2) ** 3 -
                                      (r - delta_r / 2) ** 3)
    n_ideal = (4.0 / 3.0) * np.pi * bin_volume * rho  # expected n for ideal gas

    # for i in xrange(n_skip, n_frames):
    for i in xrange(n_skip, n_skip + 10):
        sys.stdout.flush()

        if dcd_fname:
            box_mol.read_dcd_step(dcd_file, i)

        x_coor = box_mol.coor()[0, :, 0] * sigma
        y_coor = box_mol.coor()[0, :, 1] * sigma
        z_coor = box_mol.coor()[0, :, 2] * sigma

        gr_all[i] = fortran_gr.update_gr(x_coor, y_coor, z_coor, box_length[i],
                                         n_bins, delta_r)
        # gr_all[i] = update_gr(x_coor, y_coor, z_coor, box_length[i], n_bins,
                              # n_atoms, delta_r)
        gr_all[i] /= n_ideal[i]  # normalize expected n for ideal gas

    gr[:, 0] = r
    gr[:, 1] = np.mean(gr_all, axis=0) / n_atoms  # normalize by frames and atoms

    np.savetxt('gr_all_stride{}.dat'.format(stride), gr, fmt='%.14f')
    gr_cutoff = gr[gr[:, 0] < box_length.min()/2]
    np.savetxt('gr_cutoff_stride{}.dat'.format(stride), gr_cutoff, fmt='%.14f')

    return gr_cutoff


def plot_gr(gr, stride=1, show=False):
    import matplotlib.pyplot as plt

    gr_dat = np.loadtxt('argon_85K_gr.dat')

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.set_ylabel('g(r)')
    ax1.set_xlabel('r')
    scale_factor = 1.0 / gr[-1, 1]
    ax1.plot(gr[:, 0], scale_factor * gr[:, 1], color='red', lw=2)
    ax1.plot(gr_dat[:, 0], gr_dat[:, 1], color='blue', lw=2)

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
