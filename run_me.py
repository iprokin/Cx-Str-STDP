#!/usr/bin/env python2

# (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
# INRIA Rhone-Alpes
# STDP model : script that solves equations and produces plots

import numpy as np
import matplotlib
matplotlib.use('Agg') # to be able to generate plots without X server
import matplotlib.pylab as plt

from copy import deepcopy
import solve_py
from functions import nesteddict_to_lists_r, find_for_cnd, find_STDP, interpolate_and_gsmooth
from ST_vars_and_ic import ST_vars_and_ic, CaMKII_ic, si
from paramets import paramets as pars
import os

y0 = np.hstack([CaMKII_ic, ST_vars_and_ic.values()])

lv, lk = nesteddict_to_lists_r(pars)
t = np.arange(pars['integration']['t_start'], pars['integration']['t_end'], pars['integration']['t_step'])

fig_out_dir='./pngs/'
if not os.path.isdir(fig_out_dir):
    os.mkdir(fig_out_dir)

######################
##### time series ####
######################


def compute_and_plot_ca_dt(dt):
    lv1 = deepcopy(lv)
    lv1[lk.index('stimulation, dt_stim')] = dt
    y, ISTATE = solve_py.lsoda_calc(y0, lv1, t)

    fig, axs = plt.subplots(1, 2, figsize=(15, 5), sharey=True)
    for ax in axs:
        ax.plot(t, y[si('Ca_cyt'), :], 'k')
        ax.set_xlabel('Time, s')
        ax.set_title('$\Delta t_{STDP}$ = %g ms' % (lv1[lk.index('stimulation, dt_stim')]*1000))
    axs[0].set_ylabel('Calcium, $\mu M$')
    axs[0].set_xlim([0, 120])
    axs[1].set_xlim([0, 1])
    fig.savefig(fig_out_dir+'Fig_ca_dt{}ms.png'.format(dt*1000), bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

dts = [0.015, 0.03, 0.2]
dts = dts + list(-np.array(dts))
for dt in dts:
    compute_and_plot_ca_dt(dt)



cnds = map(
        dict,
        zip(
            zip(['stimulation, num_stim']*4, [10, 50, 75, 100][::-1])*2,
            zip(['stimulation, dt_stim']*4, [-0.015]*4) + zip(['stimulation, dt_stim']*4, [0.015]*4)
            )
        )

def funcs(y, lv):
    return {
            'ctrl1': pars['ECb']['kCB1R']*y.T[:, si('o_CB1R')] + pars['DA']['gamma1DA']*pars['DA']['DA'],
            'CaMKIIpho': solve_py.extras.camkiipho(y.T[:, :13]),
#            'CaM': solve_py.extras.cam_conc_func(y.T[:, si('Ca_cyt')], lv)
            }

extras = [None]*len(cnds)
for i, cnd in enumerate(cnds):
    y, lv1 = find_for_cnd(lv, lk, y0, t, cnd)
    extras[i] = funcs(y, lv1)

fig, axs0 = plt.subplots(1, 2, figsize=(12,5))
axs = dict(zip([-0.015, 0.015], axs0))
cc = map(plt.get_cmap('copper'), np.linspace(0, 1, 4))
for i, cnd in enumerate(cnds):
    ax = axs[cnd['stimulation, dt_stim']]
    ax.set_color_cycle(cc)
    ax.plot(t, extras[i]['CaMKIIpho'], label="$N_{pairings}$=%i" % cnd['stimulation, num_stim'], alpha=0.5)
    ax.set_title('$\Delta t_{STDP}$ = %g ms' % (cnd['stimulation, dt_stim']*1000))
    ax.set_xlabel('Time, s')
axs0[0].set_ylabel('CaMKII act., $\mu M$')
plt.legend()
fig.savefig(fig_out_dir+'Fig_CaMKII.png', bbox_inches='tight', pad_inches=0.1)
plt.close(fig)

fig, axs0 = plt.subplots(1, 2, figsize=(12,5))
axs = dict(zip([-0.015, 0.015], axs0))
for i, cnd in enumerate(cnds):
    ax = axs[cnd['stimulation, dt_stim']]
    ax.plot(t, extras[i]['ctrl1'], label="$N_{pairings}$=%i" % cnd['stimulation, num_stim'], alpha=0.5)
    ax.set_title('$\Delta t_{STDP}$ = %g ms' % (cnd['stimulation, dt_stim']*1000))
    ax.set_xlabel('Time, s')
axs0[0].set_ylabel('CB1R activation ($y_{CB1R}$) (a.u.)')
plt.legend()
fig.savefig(fig_out_dir+'Fig_yCB1R.png', bbox_inches='tight', pad_inches=0.1)
plt.close(fig)


######################
##### STDP curve #####
######################


def stdp_n(n, lv, dts, dts_new):
    t = np.arange(0.0, 150.0+n/pars['stimulation']['Freq'], pars['integration']['t_step'])
    # raw STDP
    dstdp = find_STDP(t, y0, si('fpre'), lv, lk.index('stimulation, dt_stim'), dts, parallel=True)
    # smoothing
    ftot = interpolate_and_gsmooth(dts, dstdp['ftot'], dts_new, 0.003)
    # clip to 3
    ftot = np.clip(ftot, 0, 3)
    return ftot

fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True, sharex=True)

dts = np.linspace(-0.04, 0.04, 80)
dts_new = np.linspace(-0.04, 0.04, 300)
lv1 = deepcopy(lv)
for i, n in enumerate([10, 100]):
    lv1[lk.index('stimulation, num_stim')] = n
    ftot = stdp_n(n, lv1, dts, dts_new)
    axs[i].plot(dts_new*1000, ftot*100)
    axs[i].set_title("$N_{pairings}$=%i" % n)

axs[0].set_xlabel('$\Delta t_{STDP}$ (ms)')
axs[0].set_ylabel('$W_{total}$ (%)')
fig.savefig(fig_out_dir+'Fig_STDP.png', bbox_inches='tight', pad_inches=0.1)
plt.close(fig)
