import os
import pandas as pd

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
                linbreg = count_fp + '/saved_objects/' + f

            # find the sensing_mx for current thres and current pris_N
            if 'pris' + str(pris_N) + '_SensingMx.file' in f:
                sensingmx = count_fp + '/saved_objects/' + f

            # find the pris obj for current thres and current pris_N
            if '_pris' + str(pris_N) + '.file' in f:
                prisobj = count_fp + '/saved_objects/' + f

        # construct the dict:
        e = {'linbreg': linbreg,
             'sensing matrix': sensingmx,
             'pris': prisobj}
        # print(e)
        d[key]['pris' + str(pris_N)] = e

    return d