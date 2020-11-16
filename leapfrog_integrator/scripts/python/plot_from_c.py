# !/usr/bin/python3
#  -*- coding: utf-8 -*-
from functions import mag_vec, state_vector_to_semimajor
from canon_units import units
from util import (MSUN, RSUN, AU, PLANETS, DAY, GYEAR, G,
                  MEARTH, YEAR, HOUR, MYEAR, KYEAR, MIN, c_light)

import numpy as np
import matplotlib.pyplot as plt

import os

uM, uL, uT = units(uL=RSUN, uM=MSUN)

time_values = np.loadtxt("./../c/times.csv",
                         delimiter=";", dtype=float, usecols=(0),
                         unpack=True)

(x_vals_star,
 y_vals_star,
 x_vals_inner,
 y_vals_inner,
 x_vals_outer,
 y_vals_outer) = np.loadtxt("./../c/positions.csv",
                            delimiter=";", dtype=float, usecols=(0, 1, 2, 3, 4, 5),
                            unpack=True)

(vx_vals_star,
 vy_vals_star,
 vx_vals_inner,
 vy_vals_inner,
 vx_vals_outer,
 vy_vals_outer) = np.loadtxt("./../c/velocities.csv",
                             delimiter=";", dtype=float, usecols=(0, 1, 2, 3, 4, 5),
                             unpack=True)


print(max(time_values) * uT / YEAR)
#################################################
# #### SET INITIAL CONDITIONS OF THE SYSTEM #####
#################################################
m0 = MSUN / uM  # mass of star
m1 = PLANETS.Saturn.M / uM  # mass of inner
m2 = PLANETS.Jupiter.M / uM  # mass of outer

######################################
# #### COMPUTE ORBITAL PROPERTIES ####
######################################
r0 = mag_vec(x_vals_star, y_vals_star)
v0 = mag_vec(vx_vals_star, vy_vals_star)

r1 = mag_vec(x_vals_inner, y_vals_inner)
v1 = mag_vec(vx_vals_inner, vy_vals_inner)
a1 = state_vector_to_semimajor(m0, m1, r1, v1)

r2 = mag_vec(x_vals_outer, y_vals_outer)
v2 = mag_vec(vx_vals_outer, vy_vals_outer)
a2 = state_vector_to_semimajor(m0, m2, r2, v2)

##################################################
# #### COMPUTE THE TOTAL ENERGY OF THE SYSTEM ####
##################################################
kinetic_energy = 0.5 * (m0 * v0**2. + m1 * v1**2. + m2 * v2**2.)

# potential_energy = - G * m1 * m0 / SI_separations - G * m2 * m0 / SO_separations - G * m1 * m2 / IO_separations
potential_energy = - G * m1 * m0 / mag_vec(x_vals_star - x_vals_inner, y_vals_star - y_vals_inner) - G * m2 * m0 / mag_vec(x_vals_star - x_vals_outer, y_vals_star - y_vals_outer) - G * m1 * m2 / mag_vec(x_vals_inner - x_vals_outer, y_vals_inner - y_vals_outer)

total_energy = kinetic_energy + potential_energy
print(kinetic_energy[0], potential_energy[0])
########################################################
# #### COMPUTE THE RELATIVE ERROR OF THE SIMULATION ####
########################################################
rel_error = np.abs((total_energy - total_energy[0]) / total_energy[0])

####################
# #### PLOTTING ####
####################
# DIRECTORY FOR SAVING IMAGES
images_dir = f"../../figures/"
os.makedirs(images_dir, exist_ok=True)

index = 1

plt.figure()
plt.plot(x_vals_star[::index] * uL / AU, y_vals_star[::index] * uL / AU, "y.", ms=0.4, label="Star")
plt.plot(x_vals_inner[::index] * uL / AU, y_vals_inner[::index] * uL / AU, "b.", ms=0.4, label="Inner planet")
plt.plot(x_vals_outer[::index] * uL / AU, y_vals_outer[::index] * uL / AU, "r.", ms=0.4, label="Outer planet")
plt.legend(loc="upper right")
plt.title("Orbit Plots")
plt.xlabel("x position [au]")
plt.ylabel("y position [au]")
plt.legend(loc="upper left", ncol=3)
plt.grid(ls="--")
fig_name = os.path.join(images_dir, "positions.png")
plt.savefig(fig_name, dpi=300)
plt.show()

plt.figure()
plt.plot(time_values[::index] * uT / YEAR, r1[::index] * uL / AU, "b-", ms=0.5, label="Inner planet")
plt.ylabel("position vector [au]", fontsize=11)
plt.xlabel("time [yr]", fontsize=11)
plt.legend(loc="upper left")
fig_name = os.path.join(images_dir, "position_vector_inner.png")
# plt.xlim((0.0, total_t * uT / YEAR))
# plt.ylim((0.00995, 0.01015))
plt.savefig(fig_name, dpi=300)
plt.show()

plt.figure()
plt.plot(time_values[::index] * uT / YEAR, a1[::index] * uL / AU, "b-", ms=0.4, label="Inner planet")
plt.ylabel("Semi-major axis [au]", fontsize=11)
plt.xlabel("time [yr]", fontsize=11)
# plt.ylim((0.00975, 0.01175))
plt.legend(loc="upper left")
fig_name = os.path.join(images_dir, "semimajor_axis_inner.png")
plt.savefig(fig_name, dpi=300)
plt.show()

plt.figure()
plt.plot(time_values[::index] * uT / YEAR, a2[::index] * uL / AU, "r-", ms=0.4, label="Outer planet")
plt.ylabel("Semi-major axis [au]", fontsize=11)
plt.xlabel("time [yr]", fontsize=11)
plt.legend(loc="upper left")
fig_name = os.path.join(images_dir, "semimajor_axis_outer.png")
plt.savefig(fig_name, dpi=300)
plt.show()

plt.figure()
plt.semilogy(time_values[::index] * uT / YEAR, rel_error[::index], "g-", ms=0.4, label="Relative error")
plt.ylabel("Relative error", fontsize=11)
plt.xlabel("time [yr]", fontsize=11)
plt.legend(loc="upper left")
fig_name = os.path.join(images_dir, "rel_error.png")
plt.savefig(fig_name, dpi=300)
plt.show()
