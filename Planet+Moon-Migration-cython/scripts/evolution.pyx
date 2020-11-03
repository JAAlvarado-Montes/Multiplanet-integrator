#!/usr/bin/env python
# -*- coding:utf-8 -*-

from bodysolver import Simulation, Body
from canon_units import units

from modulecore import (global_differential_equation,
                        omegaCritic,
                        k2Q_planet_core,
                        k2Q_planet_envelope,
                        state_vector_to_semimajor,
                        initial_velocity,
                        mag_vec)

from util import (MSUN, RSUN, AU, PLANETS, DAY, GYEAR, GCONST,
                  MEARTH, YEAR, HOUR, MYEAR, KYEAR, MIN, c_light)

import astropy.units as u
import matplotlib.pyplot as plt

import numpy as np
import time
import os
import sys


# ############################################################
# PHYSICAL PARAMETERS
# ############################################################
uM, uL, uT = units(uL=1e-3 * AU, uM=1e-12 * MSUN)

sys.stdout.write(f"Canonical units:\nuM={uM} uL={uL} uT={uT}\n")
sys.stdout.flush()


Ms = 1.0 * MSUN / uM                # Mass of the star
Rs = 1.0 * RSUN / uL                # Radius of the star

c_light = c_light * uT / uL

inner_planet = PLANETS.Saturn
e1 = 0.001  # Initial eccentricity inner planet
a1 = 0.01 * AU / uL  # Initial semimajor axis inner planet
M1 = 0.3 * inner_planet.M / uM   # Mass of the inner planet
R1 = 0.3 * inner_planet.R / uL   # Radius of the inner planet
Protini = inner_planet.Prot / uT  # Initial planetary rotational period
alpha_1 = inner_planet.alpha  # This is Jupiter's alpha
beta_1 = inner_planet.beta  # This is Jupiter's beta
rigidity_1 = 0.0 * (uL * uT**2 / uM)  # This is Jupiter's rigidity
om_1_ini = 2 * np.pi / Protini  # Initial planetary rotational rate
epsilon_planet = om_1_ini / omegaCritic(M1, R1)
k2q1_e = k2Q_planet_envelope(alpha_1,
                             beta_1,
                             epsilon_planet)
k2q1_c = k2Q_planet_core(rigidity_1, alpha_1,
                         beta_1, M1, R1)
k2q1 = k2q1_e + k2q1_c

outer_planet = PLANETS.Jupiter
e2 = 0.1  # Initial eccentricity outer planet
a2 = 0.08 * AU / uL  # Initial semimajor axis outer planet
M2 = 1.0 * outer_planet.M / uM  # Mass of the outer planet
R2 = 1.0 * outer_planet.R / uL  # Radius of the outer planet
Protini = outer_planet.Prot / uT  # Initial planetary rotational period
alpha_2 = outer_planet.alpha  # This is Jupiter's alpha
beta_2 = outer_planet.beta  # This is Jupiter's beta
rigidity_2 = 4.46E10 * (uL * uT**2 / uM)  # This is Jupiter's rigidity
om_2_ini = 2 * np.pi / Protini  # Initial planetary rotational rate
epsilon_planet = om_2_ini / omegaCritic(M2, R2)
k2q2_e = k2Q_planet_envelope(alpha_2,
                             beta_2,
                             epsilon_planet)
k2q2_c = k2Q_planet_core(rigidity_2, alpha_2,
                         beta_2, M2, R2)
k2q2 = k2q2_e + k2q2_c


# ############################################################
# FLAG TO INCLUDE EVOLUTION
# ############################################################
key = 0

# ############################################################
# Initial Conditions
# ############################################################
# Initial positions, velocities, and vectors body 1
r1_ini = a1 * (1 - e1**2.) / (1 + e1)  # Initial position inner planet
b1_x0 = r1_ini
b1_y0 = 0.0 / uL
b1_z0 = 0.0 / uL
b1_ve = initial_velocity(Ms, M1, r1_ini, a1)  # Initial tangential velocity inner planet

b1_vx0 = 0.0
b1_vy0 = b1_ve
b1_vz0 = 0.0
v1_ini = mag_vec(b1_vx0, b1_vy0, b1_vz0)

# Initial positions, velocities, and vectors body 2
r2_ini = a2 * (1 - e2**2.) / (1 + e2)  # Initial position outer planet
b2_x0 = r2_ini
b2_y0 = 0.0 / uL
b2_z0 = 0.0 / uL
b2_ve = initial_velocity(Ms, M2, r2_ini, a2)  # Initial tangential velocity outer planet

b2_vx0 = 0.0
b2_vy0 = b2_ve
b2_vz0 = 0.0
v2_ini = mag_vec(b2_vx0, b2_vy0, b2_vz0)

star = Body(name='Star',
            mass=Ms * u.kg,
            x_vec=np.array([0, 0, 0]) * u.m,
            v_vec=np.array([0, 0, 0]) * u.m / u.s,
            om_vec=np.array([0., 0., 0.]) * u.s**-1.)

inner_planet = Body(name='Inner planet',
                    mass=M1 * u.kg,
                    x_vec=np.array([b1_x0, b1_y0, b1_z0]) * u.m,
                    v_vec=np.array([b1_vx0, b1_vy0, b1_vz0]) * u.m / u.s,
                    om_vec=np.array([0., om_1_ini, 0.]) * u.s**-1.)

outer_planet = Body(name='Outer planet',
                    mass=M2 * u.kg,
                    x_vec=np.array([b2_x0, b2_y0, b2_z0]) * u.m,
                    v_vec=np.array([b2_vx0, b2_vy0, b2_vz0]) * u.m / u.s,
                    om_vec=np.array([0., om_2_ini, 0.]) * u.s**-1.)

parameters = dict(Ms=star.mass, Rs=Rs, M1=inner_planet.mass,
                  M2=outer_planet.mass, R1=R1, R2=R2,
                  key=key, k2q1=k2q1, k2q2=k2q2,
                  alpha_1=alpha_1,
                  beta_1=beta_1,
                  rigidity_1=rigidity_1,
                  alpha_2=alpha_2,
                  beta_2=beta_2,
                  rigidity_2=rigidity_2,
                  c_light=c_light)

bodies = [inner_planet, outer_planet]

a1_ini = state_vector_to_semimajor(Ms, M1, r1_ini, v1_ini)
a2_ini = state_vector_to_semimajor(Ms, M2, r2_ini, v2_ini)

print("{}\n mass={}kg, a_ini={}au\n".format(bodies[0].name,
                                            bodies[0].mass * uM,
                                            a1_ini * uL / AU))
print("{}\n mass={}kg, a_ini={}au\n".format(bodies[1].name,
                                            bodies[1].mass * uM,
                                            a2_ini * uL / AU))

# ************************************************************
# INTEGRATION
# ************************************************************
t_scale = YEAR / uT * u.s
t = 10. * t_scale

dt_scale = HOUR / uT * u.s
dt = 1. * dt_scale

N_steps = t.value / dt.value

start_time = time.time()
sys.stdout.write(f"Running n-body simulation for {N_steps} steps...")
sys.stdout.flush()

simulation = Simulation(bodies)

simulation.set_integration_method("lsoda")

simulation.set_diff_eq(global_differential_equation, **parameters)

simulation.run(t, dt)

end_time = time.time()
exec_time = (end_time - start_time) / MIN
print(f"finished in {exec_time:.3f}m")
# ************************************************************
# INTEGRATION
# ************************************************************

# ############################################################
# PLOTS
# ############################################################
# DIRECTORY FOR SAVING IMAGES
images_dir = f"../figures/{t.value * uT / YEAR}yr_dt{dt.value * uT / DAY}d"

os.makedirs(images_dir, exist_ok=True)

# Get the times and solutions of the simulation
times, solutions = simulation.history

ind = 20

r1 = mag_vec(solutions[:, 0], solutions[:, 1], solutions[:, 2])
v1 = mag_vec(solutions[:, 3], solutions[:, 4], solutions[:, 5])

r2 = mag_vec(solutions[:, 9], solutions[:, 10], solutions[:, 11])
v2 = mag_vec(solutions[:, 12], solutions[:, 13], solutions[:, 14])

a1 = state_vector_to_semimajor(Ms, M1, r1, v1)
a2 = state_vector_to_semimajor(Ms, M2, r2, v2)

fig = plt.figure(figsize=(7.5, 5.0))
ax = fig.add_subplot(1, 1, 1)

ax.plot(solutions[::ind, 0] * uL / AU, solutions[::ind, 1] * uL / AU,
        "b.", ms=0.4, label="Inner planet")
ax.plot(solutions[::ind, 9] * uL / AU, solutions[::ind, 10] * uL / AU,
        "k.", ms=0.4, label="Outer planet")

ax.set_ylabel("y [au]", fontsize=11)
ax.set_xlabel("x [au]", fontsize=11)

ax.legend(loc="upper left")
fig_name = os.path.join(images_dir, "xy_positions.png")
fig.savefig(fig_name)
# exit(0)

fig = plt.figure(figsize=(7.5, 5.0))
ax = fig.add_subplot(1, 1, 1)
ax.plot(times[::ind] * uT / YEAR, a1[::ind] * uL / AU,
        "b.", ms=0.5, label="Inner planet")

ax.set_ylabel("Semi-major axis [au]", fontsize=11)
ax.set_xlabel("time [yr]", fontsize=11)

ax.legend(loc="upper left")
fig_name = os.path.join(images_dir, "semimajor_axis_inner_planet_.png")
fig.savefig(fig_name)


fig = plt.figure(figsize=(7.5, 5.0))
ax = fig.add_subplot(1, 1, 1)
ax.plot(times[::ind] * uT / YEAR, a2[::ind] * uL / AU,
        "k.", ms=0.5, label="Outer planet")

ax.set_ylabel("Semi-major axis [au]", fontsize=11)
ax.set_xlabel("time [yr]", fontsize=11)

ax.legend(loc="upper left")
fig_name = os.path.join(images_dir, "semimajor_axis_outer_planet_.png")
fig.savefig(fig_name)

dot_prod_1 = solutions[:, 0] * solutions[:, 3] + solutions[:, 1] * solutions[:, 4] + solutions[:, 2] * solutions[:, 5]
e1_x = (v1**2. / (GCONST * Ms) - 1. / r1) * solutions[:, 0] - dot_prod_1 / (GCONST * Ms) * solutions[:, 3]
e1_y = (v1**2. / (GCONST * Ms) - 1. / r1) * solutions[:, 1] - dot_prod_1 / (GCONST * Ms) * solutions[:, 4]
e1_z = (v1**2. / (GCONST * Ms) - 1. / r1) * solutions[:, 2] - dot_prod_1 / (GCONST * Ms) * solutions[:, 5]

dot_prod_2 = solutions[:, 9] * solutions[:, 12] + solutions[:, 10] * solutions[:, 13] + solutions[:, 11] * solutions[:, 14]
e2_x = (v2**2. / (GCONST * Ms) - 1. / r2) * solutions[:, 9] - dot_prod_2 / (GCONST * Ms) * solutions[:, 12]
e2_y = (v2**2. / (GCONST * Ms) - 1. / r2) * solutions[:, 10] - dot_prod_2 / (GCONST * Ms) * solutions[:, 13]
e2_z = (v2**2. / (GCONST * Ms) - 1. / r2) * solutions[:, 11] - dot_prod_2 / (GCONST * Ms) * solutions[:, 14]

ecc_1 = mag_vec(e1_x, e1_y, e1_z)
ecc_2 = mag_vec(e2_x, e2_y, e2_z)

fig = plt.figure(figsize=(7.5, 5.0))
ax = fig.add_subplot(1, 1, 1)
plt.plot(times[::ind] * uT / YEAR, ecc_1[::ind],
         "b.", ms=0.5, label="Inner planet")

ax.set_ylabel("Eccentricity", fontsize=11)
ax.set_xlabel("time [yr]", fontsize=11)

ax.legend(loc="upper left")
fig_name = os.path.join(images_dir, "eccentricity_inner_planet_.png")
fig.savefig(fig_name)

fig = plt.figure(figsize=(7.5, 5.0))
ax = fig.add_subplot(1, 1, 1)
plt.plot(times[::ind] * uT / YEAR, ecc_2[::ind],
         "k.", ms=0.5, label="Outer planet")

ax.set_ylabel("Eccentricity", fontsize=11)
ax.set_xlabel("time [yr]", fontsize=11)

ax.legend(loc="upper left")
fig_name = os.path.join(images_dir, "eccentricity_outer_planet_.png")
fig.savefig(fig_name)
