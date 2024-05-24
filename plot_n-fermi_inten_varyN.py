#!/usr/bin/env python3

import sys
import h5py
import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

col_rainbow = [ 'red', 'tab:orange', 'green', 'blue', 'purple' ]
# col_oranges = [matplotlib.cm.Oranges(127), matplotlib.cm.Oranges(159), matplotlib.cm.Oranges(191), matplotlib.cm.Oranges(223), matplotlib.cm.Oranges(255)]
col_oranges = [matplotlib.cm.Oranges(155), matplotlib.cm.Oranges(180), matplotlib.cm.Oranges(205), matplotlib.cm.Oranges(230), matplotlib.cm.Oranges(255)]

col = col_oranges

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, nargs=1, help="where to find the .h5 file")
    parser.add_argument("-L", "--length", type=int, default=20, help="lattice sidelength")
    parser.add_argument("-p", "--particles", type=int, default=101, help="number of particles")
    # parser.add_argument("-n", "--n_flavor", type=int, default=100, help="degeneracy of fermions")
    # parser.add_argument("-U", "--U_hub", type=float, default=0.5, help="Hubbard interaction")
    parser.add_argument("-N", "--num_tsteps", type=int, default=201, help="number of time steps")
    parser.add_argument("--T_start", type=float, default=0.0, help="start time")
    parser.add_argument("--T_end", type=float, default=10.0, help="end time")
    parser.add_argument("-s", "--save", nargs="?", const="dump", default="", help="save to pdf file")
    parser.add_argument("--noview", action="store_true", default=False, help="no visual output")
    parser.add_argument("--notitle", action="store_true", default=False, help="no title on plot")
    args = parser.parse_args()
    
    path = args.path[0]
    L = args.length
    N = args.particles
    # U = args.U_hub
    # nc = args.n_flavor
    
    ncs = (2, 10, 100, 1000)
    U = "5e-1"

    dt = (args.T_end - args.T_start) / (args.num_tsteps-1)
    s_i_t = round(args.T_start / dt)
    e_i_t = round(args.T_end / dt) + 1
    
    data_fermi = np.loadtxt(path + f"/n_fermi_U{U}")
    data_inten = np.loadtxt(path + f"/n_inten_U{U}")

    # how the data is presented
    # colors = [matplotlib.cm.Oranges(155), matplotlib.cm.Oranges(180), matplotlib.cm.Oranges(205), matplotlib.cm.Oranges(230), matplotlib.cm.Oranges(255)]
    # colors = [ 'blue', 'orange', 'green', 'pink', 'brown', 'purple', 'gray', 'red', 'yellow' ]
    pfmt_pth = {'color': 'black', 'linewidth': 2.0, 'linestyle': 'dashed', 'label': 'p.th.'}
    # pfmt_num = [ {'color': colors[i_U], 'linewidth': 2.0} for i_U in range(len(Us)) ]
    pfmt_num = [ { 'linewidth': 2.0 } for i_n in range(len(ncs)) ]

    font = {'size': 20}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('axes', labelsize=27)
    plt.rcParams.update({'text.latex.preamble' : r'\usepackage{amsmath}'})
    
    fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace': 0, 'wspace': 0}, figsize=(10, 8), dpi=150)

    # if not args.notitle:
    #     ax.set_title(r'$\Delta n_{k_F}(t) = n_{k_F}(t) - n_{k_F}(0)$ [$N_{rep} = '+str(n_samples)+'$]')
    # ax.set_xlabel(r'time $t$')

    for i_n, nc in enumerate(ncs):
        label = f'$N = {nc}$'
        axs[0].plot(data_fermi[s_i_t:e_i_t,0], data_fermi[s_i_t:e_i_t,2+i_n], **pfmt_num[i_n], label=label)
        axs[1].plot(data_inten[s_i_t:e_i_t,0], data_inten[s_i_t:e_i_t,2+i_n], **pfmt_num[i_n], label=label)
        
        axs[1].set(xlabel = r'time $t/t_\text{h}^{-1}$', xlim = (args.T_start, args.T_end))

    # perturbation theory
    axs[0].plot(data_fermi[:,0], data_fermi[:,1], **pfmt_pth)
    axs[1].plot(data_inten[:,0], data_inten[:,1], **pfmt_pth)
 

        # axs[0].plot(times, nc/float(U)**2 * (fermi_data[nc][s_i_t:e_i_t,2] - fermi_data[nc][s_i_t,2]), color=col[i_n], label='$N = '+str(nc)+'$', linewidth=2.0)
        # axs[0].set(xlabel = r'time $t ~ [1/t_\text{hop}]$', ylabel = r'jump $N U^{-2} \langle \Delta n_{k_F} \rangle$', xlim = (args.T_start, args.T_end))

        # axs[1].plot(times, nc/(L*L*float(U)*float(U)) * (energy_data[nc][s_i_t:e_i_t,3] - energy_data[nc][s_i_t,3]), color=col[i_n], label='$N = '+str(nc)+'$', linewidth=2.0)
        # # axs[1].set(xlabel = r'time $t ~ [1/t_\text{hop}]$', ylabel = r'int. en. $N U^{-1} \langle \rho_{ii}^2 \rangle(t)$', xlim = (args.T_start, args.T_end))
        # axs[1].set(xlabel = r'time $t ~ [1/t_\text{hop}]$', ylabel = r'int. en. $N U^{-1} \Delta\langle n_i^2 \rangle$', xlim = (args.T_start, args.T_end))

        # axs[0].label_outer()
        # axs[1].label_outer()
    
    # pt_fermi = np.loadtxt(path + '/pt2nd_L20_U5e-1_T20/n_fermi.txt')
    # axs[0].plot(times, 1.0/0.5 * (pt_fermi[s_i_t:e_i_t,2] - pt_fermi[s_i_t,2]), color='black', linestyle='dashed', label='p.th.')

    # pt_energy = np.loadtxt(path + '/pt2nd_L20_U5e-1_T20/n_energy.txt')
    # for i_n, nc in enumerate(ncs):
    #     axs[1].plot(times, -1.0/(L*L * 0.5) * (pt_energy[s_i_t:e_i_t,2] - pt_energy[s_i_t,2]), color='black', linestyle='dashed', label='p.th.')

    axs[0].set(ylabel = 'Fermi surface')
    axs[1].set(ylabel = 'interaction energy')

    plt.text(0.2, 0.935, r'$N U^{-2} \big( n_{k_F}(t) - n_{k_F}(0) \big)$', horizontalalignment='center', verticalalignment='center', transform=axs[0].transAxes)
    plt.text(0.2, 0.935, r'$N U^{-1} \big( \langle n_i^2 \rangle(t) - \langle n_i^2 \rangle(0) \big)$', horizontalalignment='center', verticalalignment='center', transform=axs[1].transAxes)

    handles, labels = axs[0].get_legend_handles_labels()
    # fig.legend(handles, labels, loc='upper center', ncol=5, fontsize=18)
    # fig.legend(handles, labels, loc=(0.2, 0.6), ncol=2, fontsize=22)
    axs[0].legend(handles, labels, loc='lower left', ncol=2, frameon=False, fontsize=22)

    plt.subplots_adjust(left=0.12, right=0.985, top=0.98, bottom=0.115)

    if not args.save == "":
        plt.savefig(args.save + ".pdf", format='pdf')
        plt.savefig(args.save + ".eps", format='eps')

    if not args.noview:
        plt.show()

