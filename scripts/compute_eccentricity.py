import numpy as np

import numpy as np


def compute_eccentricity(r, v, mu):
    rn = np.linalg.norm(r)
    vn = np.linalg.norm(v)
    print(f"rn = {rn}")
    print(f"vn = {vn}")
    v_r = np.dot(r / rn, v)
    print(f"v_r = {v_r}")
    v_p = np.sqrt(vn**2 - v_r**2)

    h = np.cross(r, v)
    print(h)
    hn = np.linalg.norm(h)
    print(hn)

    vxh = np.cross(v, h)
    vxh_mu = vxh / mu
    norm_r = r / rn
    print(f"vxh = {vxh}")
    print(f"vxh_mu = {vxh_mu}")
    print(f"norm_r = {norm_r}")
    e_vec = vxh_mu - norm_r
    e = np.linalg.norm(e_vec)
    print(f"e = {e_vec}")
    return e

r_km = np.array((1_000, 5_000, 7_000))  # km
v_km = np.array((3.0, 4.0, 5.0))  # km/s

r_m = r_km * 1000.0
v_m = v_km * 1000.0
mu_km = 3.986004418e5  # km^3/s^2
mu_m = 3.986004418e14  # m^3/s^2

e_km = compute_eccentricity(r_km, v_km, mu_km)
e_m = compute_eccentricity(r_m, v_m, mu_m)

print(f"e_km = {e_km}\ne_m  = {e_m}")

"""
rn = 8660.254037844386
vn = 7.0710678118654755
v_r = 6.697263122599659

h = 
  -3.0000e+03
   1.6000e+04
  -1.1000e+04
hn = 19646.8827043885

rn = 8660254.037844386
vn = 7071.067811865475
v_r = 6697.263122599659

h = 
  -3.0000e+09
   1.6000e+10
  -1.1000e+10
hn = 19646882704.3885

vxh = 
  -1.2400e+14
   1.8000e+13
   6.0000e+13
vxh_mu = 
  -0.0311
   0.0045
   0.0151
norm_r = 
   0.1155
   0.5774
   0.8083

e = 
  -0.1466
  -0.5728
  -0.7932

e_km = 0.9893688725705063
e_m  = 0.9893688725705063
"""