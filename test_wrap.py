from IPython import get_ipython; get_ipython().magic('reset -sf')
import itasca as it
it.command("python-reset-state false")
import numpy as np
it.command('model restore "PFCmat-mat.p3sav"')
#it.command('model new')


inner_radius = 0.125
outer_radius = inner_radius*3
thickness = 0.05

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
  size 20 5 20 ratio 1.1 1.1 1


zone reflect origin 0 0 0 norm -1 0 0
zone reflect origin 0 0 0 norm 0 -1 0
""")


it.command(f"""
zone cmodel assign elastic
zone property young [pbm_emod+lnm_emod] poisson 0.25
zone property density [cm_densityVal]
zone face group "boundary" range cylinder end-1 0 0 {-thickness/2.0} end-2 0 0 {thickness/2.0} rad {inner_radius}
zone group "zboundary" range cylinder end-1 0 0 {-thickness/2.0} end-2 0 0 {thickness/2.0} rad {1.1*inner_radius}
wall-zone create name 'dem_boundary' group 'zboundary'

""")

