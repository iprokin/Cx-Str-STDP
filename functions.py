# (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
# INRIA Rhone-Alpes
# STDP model : helper functions

from copy import deepcopy
import numpy as np
import solve_py

def pmap(f, xs):
    from multiprocessing import cpu_count, Pool
    p = Pool(cpu_count())
    Y = p.map(f, xs)
    p.close()
    p.join()
    return Y


def interpolate_and_gsmooth(x, y, xnew, sdx):
    # x,y,xnew, sdx - in units of x/xnew
    from scipy.ndimage.filters import gaussian_filter
    ynew = np.interp(xnew, x, y)
    ys = gaussian_filter(ynew, sdx/(xnew[1]-xnew[0]))
    return ys


def nesteddict_to_lists_r(d, keys=None):
    from collections import OrderedDict
    """converts nested dict to the list of keys and the list of values
    keys - is the selecton of top level keys to be considered
    (all top keys taken if not set)
    """
    lv, lk = [], []
    if keys is None:
        keys = d.keys()
    for k in keys:
        if type(d[k]) in [dict, OrderedDict]:
            lvx, lkx = nesteddict_to_lists_r(d[k])
            lv += lvx
            lk += [k+', '+lkx1 for lkx1 in lkx]
        else:
            lv.append(d[k])
            lk.append(k)
    return lv, lk


def find_y0(y0, lv, lk, t_end=200, t_step=1e-3):
    lv0 = deepcopy(lv)
    for ind in [lk.index('stimulation, pre_on'), lk.index('stimulation, post_on'), lk.index('stimulation, num_stim')]:
        lv0[ind] = 0
    t = np.arange(0, t_end, t_step)
    y, ISTATE = solve_py.lsoda_calc(y0, lv0, t)
    return y[:,-1] # that's it!


def find_for_cnd(lv, lk, y0, t, cnd):
    lv1 = deepcopy(lv)
    lv1 = np.array(lv1)
    lv1[map(lambda k: lk.index(k), cnd.keys())] = cnd.values()
    y, ISTATE = solve_py.lsoda_calc(find_y0(y0, lv1, lk), lv1, t)
    return y, lv1


def fpre_fpost_y_4F(y, sifpre,
                    fpost_func=lambda CaMKact_t: CaMKact_t[-1]*3.5/164.6+1):
    y = y[:, [0, -1]]
    CaMKact_t = solve_py.extras.camkiipho_4f(y[:13, :])
    fpost = fpost_func(CaMKact_t)
    fpre = y[sifpre, -1]
    return fpre, fpost


def find1point(t, y0, sifpre, parvals0, idx_dt, dt,
               fpre_fpost_y_4F=fpre_fpost_y_4F):
    parvals1 = deepcopy(parvals0)
    parvals1[idx_dt] = dt
    y, ISTATE = solve_py.lsoda_calc(y0, parvals1, t)
    if ISTATE > 0:  # no errors
        return fpre_fpost_y_4F(y, sifpre)
    else:
        return np.nan, np.nan


def find_STDP(t, y0, sifpre, parvals0, idx_dt, dts,
              parallel=True,
              fpre_fpost_y_4F=fpre_fpost_y_4F):
    from functools import partial
    result = dict()
    if parallel:
        umap = pmap
        parvals1 = parvals0
    else:
        umap = map
        parvals1 = deepcopy(parvals0)
    X = umap(partial(find1point,
                     t, y0, sifpre, parvals1, idx_dt,
                     fpre_fpost_y_4F=fpre_fpost_y_4F),
             dts)
    X = map(lambda *a: list(a), *X)
    result['fpre'], result['fpost'] = np.array(X[0]), np.array(X[1])
    result['ftot'] = result['fpost'] * result['fpre']
    return result

