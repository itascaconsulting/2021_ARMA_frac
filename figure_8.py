import matplotlib; matplotlib.rcParams["savefig.directory"] = "."
import numpy as np
from matplotlib import pyplot as plt
import os

for i in range(5):
    for j in range(5):
        k = i*5 + j
        geom_file = f"evis{k}.geom"
        if os.path.isfile(geom_file):
            x,y = np.loadtxt(geom_file, skiprows=2, usecols=(1,2)).T
            plt.subplot(5,5,k+1)
            plt.scatter(x,y,marker="o", color="black", s=0.25, edgecolor="black")
            plt.gca().set_aspect(1)
            plt.ylim(-0.35,0.35)
            plt.xlim(-0.35,0.35)
            plt.xticks([])
            plt.yticks([])
            #plt.gca().set_frame_on(False)
plt.show()
