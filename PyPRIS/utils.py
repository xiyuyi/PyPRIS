import os
import pandas as pd
import numpy as np
import copy
from csv import writer
from sklearn.cluster import DBSCAN
import pickle

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


def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)


def get_clusters(locs, brightnesses, epsv, ticketIndex, csvpath, saveoption):
    db = DBSCAN(eps=epsv, min_samples=1).fit(locs)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    clusters = []
    # now check the kth cluster
    for k in set(labels):
        class_member_mask = (labels == k)
        xy = locs[class_member_mask & core_samples_mask]
        br = brightnesses[class_member_mask & core_samples_mask]
        xyc = np.mean(xy, axis=0)
        brc = np.sum(br)
        clusters.append([xyc[0], xyc[1], brc])
        therow = [ticketIndex, xyc[0], xyc[1], brc]
        if saveoption is True:
            append_list_as_row(file_name=csvpath, list_of_elem=therow)
    return clusters


def get_clusters_fromlinbreg(lbpath, ticketIndex, csvpath, D, save=False):
    """
    get clusters from the full path of a single linbreg object.
    :param lbpath: full path to the targeted linbreg object
    :param ticketIndex: this will be the index of the ticket that goes into the csv vile.
    :param csvpath: full path of the csv file to append the fitting results
    :param D: dimension [int]
    :return: fitted clusters as a list
    """
    if D == 2:
        # perform clustering using 2D only
        with open(lbpath, 'rb') as f:
            lb = pickle.load(f)
        locs=np.asarray(lb.candidate_coords)[:,1:3]
        brightnesses=lb.x[0:-1]
        epsv=np.linalg.norm(lb.candidate_intervals[1:])
        clusters = get_clusters(locs,brightnesses, epsv=epsv, ticketIndex=ticketIndex, csvpath=csvpath, saveoption=save)
    else:
        print('Sorry, only 2D clustering is available as of yet')
    return clusters


def find_psf_matrix_offset(pypris):
    """
    find the offset of the PSF matrix from the exact center of the matrix.
    :param pypris:
    :return:
    """
    center = np.round((np.max(pypris.current_candidates, axis=0) - np.min(pypris.current_candidates, axis=0)) / 2).astype('float32')+0.5
    range_ind0 = np.arange(center[0], center[0] + 1 - pypris.current_candidates_intervals[0]/2, pypris.current_candidates_intervals[0])
    range_ind1 = np.arange(center[1], center[1] + 1 - pypris.current_candidates_intervals[1]/2, pypris.current_candidates_intervals[1])
    range_ind2 = np.arange(center[2], center[2] + 1 - pypris.current_candidates_intervals[2]/2, pypris.current_candidates_intervals[2])
    pypris.current_candidates = list()  # erase the current pool of candidates, and make the new set
    if pypris.species_n == 1:
        for i0 in range_ind0:
            for i1 in range_ind1:
                for i2 in range_ind2:
                    pypris.current_candidates.append([i0, i1, i2, 1])

    # Now generate the observation of all those candidates locations by generating the sensing matrix
    pypris.generate_sensing_mx()

    # find the candidate that is located at a center pixel in observation space.
    # sort observation pixels
    l = (np.sort(pypris.current_A[:, :-1], axis=0))
    # calculate the ratio of the brightese pixel and the second brightest pixel
    p1 = l[-1, :]
    p2 = l[-2, :]
    ratio = list(p1 / p2)

    # find the index with the highest ratio
    ind = ratio.index(np.max(ratio))

    # print out the offset
    offset = pypris.current_candidates[ind][0:3] - center
    return offset