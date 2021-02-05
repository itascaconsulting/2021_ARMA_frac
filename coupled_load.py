import itasca as it
from itasca import ballarray as ba
from vec import vec
import numpy as np
import pylab as plt
import math

class cube_blast(object):
    def __init__(self):
        pass
    def setup(self, material, rise_time, peak_pressure,
                 pressure_decay, hole_radius,case,prefix):
        it.command("""
        model new
        MODEL LARGE-STRAIN on
        model res '{savefile}'
        ball attr damp 0
        ball ini disp 0
        ball delete range cyl end-1 0 0 -100 end-2 0 0 100 rad {hole_radius}
        contact prop dp_nratio 0.1
        contact prop dp_sratio 0.1
        ball prop "dp_nratio" 0.1
        ball prop "dp_sratio" 0.1

        """.format(savefile=material,hole_radius=hole_radius))


        inner_radius = it.fish.get("mv_W")/2.0
        assert it.fish.get("mv_W") == it.fish.get("mv_D")
        outer_radius = inner_radius*5
        thickness = it.fish.get("mv_H")

        it.command(f"""
        model config dynamic
        model domain extent {-1.1*outer_radius} {1.1*outer_radius} {-1.1*outer_radius} {1.1*outer_radius} {-1.1*thickness} {1.1*thickness} condition destroy
        zone create cylindrical-shell ...
          point 0 0                 {0}              {-thickness/2.0} ...
          point 1 {outer_radius}    {0}              {-thickness/2.0} ...
          point 2 0                 {0}              {thickness/2.0} ...
          point 3 0                 {-outer_radius}  {-thickness/2.0} ...
          point 8 {inner_radius}    {0}              {-thickness/2.0} ...
          point 9 0                 {-inner_radius}  {-thickness/2.0} ...
          point 10 {inner_radius}    {0}              {thickness/2.0} ...
          point 11 0                 {-inner_radius}  {thickness/2.0} ...
          size 25 4 25 ratio 1 1 1


        zone reflect origin 0 0 0 norm -1 0 0
        zone reflect origin 0 0 0 norm 0 -1 0

        zone cmodel assign elastic
        zone property young [pbm_emod+lnm_emod] poisson 0.25 
        ;zone cmodel assign mohr-coulomb-tension
        ;zone property young [pbm_emod+lnm_emod] poisson 0.25 cohesion {it.fish.get("pbm_coh_m")} tension {it.fish.get("pbm_ten_m")} friction 50 dilation 0 number-cracks 1
        
        ;zone cmodel assign mohr-coulomb
        ;zone property young [pbm_emod+lnm_emod] poisson 0.25 cohesion 1e100 tension 1e5 friction 50 dilation 0         
        ;zone cmodel assign strain-softening        
        ;table 'ConTen_C35' add 0.0000000 .75e6 ...
        ;0.0000005 1.0000 ...
        ;0.0000001 1.0000 ...
        ;0.0002000 1.0000
        ;zone property table-tension 'ConTen_C35' young [pbm_emod+lnm_emod] poisson 0.25 cohesion 1e100 friction 40 tension 0.75e6

        zone property density [cm_densityVal]
        wall-zone create name 'dem_boundary' range cylinder end-1 0 0 {-thickness/2.0} end-2 0 0 {thickness/2.0} rad {inner_radius}

        zone gridpoint fix velocity-z range position-z {thickness/2.0}
        zone gridpoint fix velocity-z range position-z {-thickness/2.0}
        zone face apply quiet-normal range cylinder end-1 0 0 {-thickness/2.0} end-2 0 0 {thickness/2.0} rad {.99*outer_radius} not
        zone face apply velocity-strike 0 range cylinder end-1 0 0 {-thickness/2.0} end-2 0 0 {thickness/2.0} rad {.99*outer_radius} not
        zone face apply velocity-dip 0 range cylinder end-1 0 0 {-thickness/2.0} end-2 0 0 {thickness/2.0} rad {.99*outer_radius} not

        contact cmat default model linearpbond property pb_ten 1e100 pb_coh 1e100 method deformability emod {it.fish.get('mv_emod')} kratio 1.0

        model clean all
        contact method bond gap 0 {0.1*it.ball.find(1).radius()} range contact type 'ball-facet'
        contact property lin_mode 1 pb_ten 1e100 pb_coh 1e100 range contact type 'ball-facet'
        """)

        self.case = case
        self.prefix = prefix
        self.rise_time = rise_time
        self.peak_pressure = peak_pressure
        self.pressure_decay = pressure_decay
        self.rad = float(ba.radius().mean())
        self.barea = self.rad**2 * math.pi

        x,y,z = ba.pos().T
        annulus = np.sqrt(x**2 + y**2) < (hole_radius + 2*self.rad)
        ba.set_extra(1, 1.0*annulus)
        # list of ball objects which will carry load
        self.load_balls = [it.ball.find(i) for i in ba.ids()[np.nonzero(annulus)[0]]]
        # initiall bonds
        self.initial_bond_list = set(((c.end1().id(), c.end2().id()) for c in \
                                 it.contact.list("mechanical") if
                                      (c.model()=='linearpbond' and
                                       c.prop("pb_state")==3)))

        self.starting_time = it.mech_age()
        self.data = []
        it.set_callback("apply_gas_force",1)

    def run(self, deltat):
        it.command("solve age {}".format(it.mech_age() + deltat))

    def save(self):
        it.command("model save '{}blast{}.p3sav'".format(self.prefix, self.case))
        np.savetxt("{}blast{}.txt".format(self.prefix, self.case), self.data)

    def show(self):
        t,p,c = np.array(self.data).T
        plt.subplot(211)
        plt.plot(t,p)
        plt.subplot(212)
        plt.plot(t,c)
        plt.show()

def pressure(time):
    "hole gas pressure as a function of time since explosion start."
    if time <= blast.rise_time:
        return blast.peak_pressure * time/blast.rise_time
    else:
        tstar = time - blast.rise_time
        return blast.peak_pressure * math.exp(-blast.pressure_decay * tstar)


def apply_gas_force(*args):
    dtime = it.mech_age() - blast.starting_time
    press = pressure(dtime)
    magnitude = press * blast.barea
    for ball in blast.load_balls:
        direction = (ball.pos()*vec((1,1,0))).norm()
        ball.set_force_app(direction * magnitude)

def broken_bonds():
    current_bonds = set(((c.end1().id(), c.end2().id()) for c in \
                         it.contact.list("mechanical") if
                            (c.model()=='linearpbond' and
                             c.prop("pb_state")==3)))
    return blast.initial_bond_list - current_bonds


def show_cracks(geom_name):
    cracks = broken_bonds()
    output = open(geom_name, "w")
    print("ITASCA GEOMETRY3D", file=output)
    print("NODES ; id x y z EXTRA 1 value", file=output)
    for i, c in enumerate(cracks):
        b1,b2 = c
        pos = 0.5*(it.ball.find(b1).pos() + it.ball.find(b2).pos())
        print("{} {} {} {} EXTRA 1 0".format(i+1,*pos), file=output)
    output.close()
    it.command("geom delete")
    it.command(f"geom import '{geom_name}' format geometry")
    return len(cracks)

# singleton
blast = cube_blast()
