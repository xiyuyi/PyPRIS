import os
import pandas as pd
import numpy as np
import copy
try:
    from matplotlib import pyplot as plt
    plt.switch_backend('agg')
except RuntimeError:
    pass

def fetch_saved_objects(path, key):
    objs = next(os.walk(path))[2]
    # first, find out how many pris iterations are there based on pris objects.
    prisn = -1
    for f in objs:
        i = f.split('pris')[-1].split('.')[0]
        if len(i.split('_')) == 1:
            try:
                prisn = np.max([prisn, int(i)])
            except:
                pass

    d = pd.DataFrame([x for x in np.arange(prisn + 1)], index=['pris' + str(x) for x in np.arange(prisn + 1)],
                     columns=[key])
    d[key]['pris0'] = 'a'  # dont know why, have to initiate this, otherwise the first element would be nan.
    plist = [str(x) for x in np.arange(7)]
    for pris_N in plist:
        # find the linbreg for current thres and current pris_N
        linbreg_n = 0
        for f in objs:
            if 'pris' + str(pris_N) in f:
                try:
                    lbn_str = f.split('.')[-2].split('_')[-1]  # linbreg iteration number string.
                    try:
                        linbreg_n = np.max([int(lbn_str), linbreg_n])
                    except:
                        pass
                except:
                    pass

        for f in objs:
            # find the linbreg for current thres and current pris_N
            if 'pris' + str(pris_N) + '_' + str(linbreg_n) + '.file' in f:
                linbreg = path + '/' + f

            # find the sensing_mx for current thres and current pris_N
            if 'pris' + str(pris_N) + '_SensingMx.file' in f:
                sensingmx = path + '/' + f

            # find the pris obj for current thres and current pris_N
            if '_pris' + str(pris_N) + '.file' in f:
                prisobj = path + '/' + f

        # construct the dict:
        e = {'linbreg': linbreg,
             'sensing matrix': sensingmx,
             'pris': prisobj}
        # print(e)
        d[key]['pris' + str(pris_N)] = e

    return d


def linbreg_report(linbreg, msg):
    v = linbreg.candidate_vis()
    vis = v[:, :, :]
    prj_ax0 = copy.deepcopy(np.mean(vis, axis=0))
    prj_ax1 = copy.deepcopy(np.mean(vis, axis=1))
    prj_ax2 = copy.deepcopy(np.mean(vis, axis=2).T)
    patch = np.zeros((vis.shape[0], vis.shape[0]))
    cat1 = np.concatenate([prj_ax0, prj_ax2], axis=1)
    cat2 = np.concatenate([prj_ax1, patch], axis=1)
    cat = np.concatenate([cat1, cat2], axis=0)

    b = linbreg.b.reshape(linbreg.obs_dim0, linbreg.obs_dim1)
    recb = linbreg.recb.reshape(linbreg.obs_dim0, linbreg.obs_dim1)
    print(msg)
    plt.figure(figsize=(15, 5))
    plt.subplot(1, 3, 1)
    plt.imshow(cat)
    t = plt.title('projections')

    plt.subplot(1, 3, 2)
    plt.imshow(b)
    t = plt.title('input blur')

    plt.subplot(1, 3, 3)
    plt.imshow(recb)
    t = plt.title('recovered blur '+str(linbreg.PyPRIS_iter)+', it'+str(linbreg.it_count))
    return [b, recb, cat, prj_ax0, prj_ax1, prj_ax2]


def candidates_init_pop_and_roll(p):
    """
    this function should go over candidates pools in a tiling fashion.
    For each tile, the corresponding sensing matrix is going to be calculated, then multiply the observation with the
    inverse of the sensing matrix, and then pop out empty candidates.
    Refine, and pop. until it reaches a desired candidate pool size (user defined.).
    Then update the pypris object with the current ReF, A, candidates, candidates intervals, etc.
    And return the updated pypris object.

    :param p: a PyPRIS object.
    :return: a modified PyPRIS object with the following updated fields:
        current_relReF (should be an appropriate number), because we want to smart initiate the pool with original pixel size.
        current_A: should reduce to the lean candidate pool.
        current_candidates
        current_candidates_intervals

    """
    print('OK')
    pass
