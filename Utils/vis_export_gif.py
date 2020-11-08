import imageio
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
thedir = 'G:\\DH_localization\\PyPRIS_tickets_F1_F2_results5.5\\PyPRIS_tickets_set6'
f = [name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name))]

for fitem in f:
    path_d= thedir + "/"+fitem+"/saved_objects"
    r1 = imageio.get_writer(path_d + '_3projections_pris.gif', mode='I')
    for prisIter in np.arange(0,6):
        r2 = imageio.get_writer(path_d + '_plots_pris'+str(prisIter)+'.gif', mode='I')
        fnames1 = list()
        fnames2 = list()
        for itN in np.arange(1,908002,2000):
            fnames1.append("PyPRIS__Proj_vies_"+fitem+"_pris"+str(prisIter)+"_plots_it"+str(itN)+".png")
            fnames2.append("PyPRIS_"+fitem+"_pris"+str(prisIter)+"_plots_it"+str(itN)+"visualize.png")

        for fname in fnames1:
            try:
                im = imageio.imread(path_d+"/"+fname)
                f = plt.figure()
                ax = f.add_subplot(111)
                plt.imshow(im)
                plt.text(-0.3, 1.1,fname, ha='left', va='top', transform=ax.transAxes)
                # Make a random plot...
                # If we haven't already shown or saved the plot, then we need to
                # draw the figure first...
                f.canvas.draw()
                # Now we can save it to a numpy array.
                data = np.frombuffer(f.canvas.tostring_rgb(), dtype=np.uint8)
                data = data.reshape(f.canvas.get_width_height()[::-1] + (3,))
                plt.close()
                r1.append_data(data)
            except:
                pass

        for fname in fnames2:
            try:
                im = imageio.imread(path_d+"/"+fname)
                f = plt.figure()
                ax = f.add_subplot(111)
                plt.imshow(im)
                plt.text(-0.3, 1.1,fname, ha='left', va='top', transform=ax.transAxes)
                # Make a random plot...
                # If we haven't already shown or saved the plot, then we need to
                # draw the figure first...
                f.canvas.draw()
                # Now we can save it to a numpy array.
                data = np.frombuffer(f.canvas.tostring_rgb(), dtype=np.uint8)
                data = data.reshape(f.canvas.get_width_height()[::-1] + (3,))
                plt.close()
                r2.append_data(data)
            except:
                pass
