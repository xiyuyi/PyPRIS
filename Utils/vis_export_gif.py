import imageio
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np

fitem="bgSCF1.5_mu7.5e+09_alpha1.0e-09"
path_d="G:/DH_localization/PyPRIS_tickets_set4/"+fitem+"/saved_objects"
r = imageio.get_writer(path_d + '_3projections.gif', mode='I')
#r = imageio.get_writer(path_d + '_plots.gif', mode='I')


fnames=list()
for itN in np.arange(1,908002,2000):
    fnames.append("PyPRIS__Proj_vies_"+fitem+"_pris0_plots_it"+str(itN)+".png")
    #fnames.append("PyPRIS_"+fitem+"_pris0_plots_it"+str(itN)+"visualize.png")

for fname in fnames:
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
        r.append_data(data)
    except:
        pass
# produce a series of images into
# with imageio.get_writer(path_d+'_3projections_movie.gif', mode='I') as writer:
#     for filename in filenames:
#         image = imageio.imread(filename)
#         writer.append_data(image)d_data(image)