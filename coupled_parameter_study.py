from IPython import get_ipython; get_ipython().magic('reset -sf')
import itasca as it
it.command("python-reset-state false")
import numpy as np
import pylab as plt
import math
from scipy.constants import inch
from importlib import reload
import coupled_load
coupled_load = reload(coupled_load)

blast, apply_gas_force, pressure = coupled_load.blast, coupled_load.apply_gas_force, coupled_load.pressure

ct = 3e-4 # characteristic time
cp = 1e6  # characteristic pressure

material = "disc-mat.p3sav"

hole_radius = 1.5*inch/2.0

prefix = "dvis"
c=0
for pmult in [2.0]:
    for rmult in [0.1]:
        rise_time = ct * rmult
        peak_pressure = cp * pmult
        pressure_decay_time = 5*ct
        decay = -math.log(0.5)/pressure_decay_time

        blast.setup(material, rise_time, peak_pressure,
                    decay, hole_radius, c, prefix)
        run_time = rise_time + 1.5e-3
        print ("running case {} for {}".format(c, run_time))
        blast.run(run_time)
        blast.save()
        c += 1
        coupled_load.show_cracks()