from util import (MSUN, RSUN, AU, PLANETS, DAY, GYEAR, G,
                  MEARTH, YEAR, HOUR, MYEAR, KYEAR, MIN, c_light)
from canon_units import units
from functions import (mag_vec, mag_vec_soft,
                       grav_force, initial_velocity,
                       state_vector_to_semimajor, grav_relat_force)

from tqdm import tqdm, trange

import numpy as np
import matplotlib.pyplot as plt
import time
import os
import sys

method_flag = sys.argv[1]

###################################################
# #### COMPUTE CANONICAL UNITS FOR SIMULATION #####
###################################################
uM, uL, uT = units(uL=RSUN, uM=MSUN)

#################################################
# #### SET INITIAL CONDITIONS OF THE SYSTEM #####
#################################################
m0 = MSUN / uM  # mass of star
m1 = PLANETS.Saturn.M / uM  # mass of inner
m2 = PLANETS.Jupiter.M / uM  # mass of outer

c_light = c_light * uT / uL

e_inner = 0.0  # Initial eccentricity inner planet
a_inner = 0.01 * AU / uL  # Initial semimajor axis inner planet
initial_position_inner = a_inner * (1 - e_inner**2.) / (1 + e_inner)  # Initial position inner planet
v_inner_initial = initial_velocity(m0, m1, initial_position_inner, a_inner)  # Initial tangential velocity

e_outer = 0.1  # Initial eccentricity outer planet
a_outer = 0.08 * AU / uL  # Initial semimajor axis outer planet
initial_position_outer = a_outer * (1 - e_outer**2.) / (1 + e_outer)  # Initial position outer planet
v_outer_initial = initial_velocity(m0, m2, initial_position_outer, a_outer)  # Initial tangential velocity

initial_position_inner = np.array((initial_position_inner, 0))  # Initial position vector of inner, with star at origin
initial_velocity_inner = np.array((0, v_inner_initial))  # Initial velocity vector of inner (travelling anticlockwise)

initial_position_outer = np.array((initial_position_outer, 0))  # Initial position vector of outer, with star at origin
initial_velocity_outer = np.array((0, v_outer_initial))  # Initial velocity vector of outer (travelling anticlockwise)

initial_position_star = np.array((0, 0)) / uL  # Still in star's rest frame
initial_velocity_star = np.array((0, 0)) * uT / uL


#############################################################
# #### COMPUTE INITIAL COORDINATES OF THE CENTER OF MASS ####
#############################################################
COM_x_initial = (m0 * initial_position_star[0] + m1 * initial_position_inner[0] + m2 * initial_position_outer[0]) / (m0 + m1 + m2)        # Finding the initial center of mass
COM_y_initial = (m0 * initial_position_star[1] + m1 * initial_position_inner[1] + m2 * initial_position_outer[1]) / (m0 + m1 + m2)
COM_initial = np.array((COM_x_initial, COM_y_initial))  # initial COM vector

initial_position_inner = initial_position_inner - COM_initial  # Moving to COM frame
initial_position_outer = initial_position_outer - COM_initial  # Moving to COM frame
initial_position_star = initial_position_star - COM_initial   # Moving to COM frame

######################################
# #### INTEGRATION TIME AND STEP #####
######################################
total_t = 10 * DAY / uT
dt = 0.1 * HOUR / uT
dt_store = 10 * HOUR / uT

time_vector = np.arange(0, total_t, dt)

if method_flag == "while":
    num_steps = 2

if method_flag == "for":
    num_steps = len(time_vector)
    time_values = time_vector

print("Integration time =", total_t * uT / YEAR, "y", end=" ")
print("timestep =", dt * uT / DAY, "d")

start_time = time.time()
print("Running n-body simulation...\n")

############################################
# #### INITIAL POSITIONS FOR THE SYSTEM ####
############################################
inner_positions = np.zeros((num_steps, 2))
inner_positions[0] = initial_position_inner

outer_positions = np.zeros((num_steps, 2))
outer_positions[0] = initial_position_outer

star_positions = np.zeros((num_steps, 2))
star_positions[0] = initial_position_star

############################################
# #### INITIAL VELOCITIES FOR THE SYSTEM ####
############################################
inner_velocities = np.zeros((num_steps, 2))
inner_velocities[0] = initial_velocity_inner

outer_velocities = np.zeros((num_steps, 2))
outer_velocities[0] = initial_velocity_outer

star_velocities = np.zeros((num_steps, 2))
star_velocities[0] = initial_velocity_star

#########################################
# #### SEPARATION VECTORS FOR INNER ####
#########################################
SI_separation_vectors = np.zeros((num_steps, 2))   # FROM star TO inner
SI_separation_vectors[0] = inner_positions[0] - star_positions[0]
IS_separation_vectors = - SI_separation_vectors

SI_separations = np.zeros((num_steps))
SI_separations[0] = mag_vec_soft(SI_separation_vectors[0][0], SI_separation_vectors[0][1])

########################################
# #### SEPARATION VECTORS FOR OUTER ####
########################################
SO_separation_vectors = np.zeros((num_steps, 2))   # FROM star TO inner
SO_separation_vectors[0] = outer_positions[0] - star_positions[0]
OS_separation_vectors = - SO_separation_vectors

SO_separations = np.zeros((num_steps))
SO_separations[0] = mag_vec_soft(SO_separation_vectors[0][0], SO_separation_vectors[0][1])

###############################################
# #### SEPARATIONS VECTORS FOR INNER-OUTER ####
###############################################
IO_separation_vectors = np.zeros((num_steps, 2))   # FROM inner TO outer
IO_separation_vectors[0] = outer_positions[0] - inner_positions[0]
OI_separation_vectors = - IO_separation_vectors

IO_separations = np.zeros((num_steps))
IO_separations[0] = mag_vec_soft(IO_separation_vectors[0][0], IO_separation_vectors[0][1])

########################################
# #### INITIAL FORCES OF THE SYSTEM ####
########################################
initial_force_inner = grav_force(m1, m0, SI_separation_vectors[0], SI_separations[0]) + grav_force(m1, m2, OI_separation_vectors[0], IO_separations[0]) + grav_relat_force(m1, m0, c_light, SI_separation_vectors[0], SI_separations[0], inner_velocities[0]) * (m0 + m1) / m0
initial_force_outer = grav_force(m2, m0, SO_separation_vectors[0], SO_separations[0]) + grav_force(m2, m1, IO_separation_vectors[0], IO_separations[0]) + grav_relat_force(m1, m0, c_light, SI_separation_vectors[0], SI_separations[0], inner_velocities[0]) / m0
initial_force_star = grav_force(m0, m1, IS_separation_vectors[0], SI_separations[0]) + grav_force(m0, m2, OS_separation_vectors[0], SO_separations[0])

#######################################
# #### FORCE VECTORS FOR THE INNER ####
#######################################
inner_forces = np.zeros((num_steps, 2))
inner_forces[0] = initial_force_inner

###################################
# #### FORCE VECTORS FOR OUTER ####
###################################
outer_forces = np.zeros((num_steps, 2))
outer_forces[0] = initial_force_outer

######################################
# #### FORCE VECTORS FOR THE STAR ####
######################################
star_forces = np.zeros((num_steps, 2))
star_forces[0] = initial_force_star

COM = np.zeros((num_steps, 2))
COM[0] = COM_initial

if method_flag == "while":

    t = np.float(0.0)
    t_save = np.float(0.0)

    SI_sep = []
    SO_sep = []
    IO_sep = []

    s_positions = [[], []]
    i_positions = [[], []]
    o_positions = [[], []]

    s_velocities = [[], []]
    i_velocities = [[], []]
    o_velocities = [[], []]

    time_values = []

    runs = int(total_t / dt)
    # pbar = tqdm(total=runs)

    for i in trange(len(time_vector) - 1):
    # while t <= total_t:
        # pbar.update(1)
        t = i * dt
        if t >= t_save:
            time_values.append(t_save)

            # Positions of the star
            s_positions[0].append(star_positions[0][0])
            s_positions[1].append(star_positions[0][1])

            # Positions of the inner planet to store
            i_positions[0].append(inner_positions[0][0])
            i_positions[1].append(inner_positions[0][1])

            # Positions of the outer planet to store
            o_positions[0].append(outer_positions[0][0])
            o_positions[1].append(outer_positions[0][1])

            # Velocities of the star to store
            s_velocities[0].append(star_velocities[0][0])
            s_velocities[1].append(star_velocities[0][1])

            # Velocities of the inner planet
            i_velocities[0].append(inner_velocities[0][0])
            i_velocities[1].append(inner_velocities[0][1])

            # Velocities of the outer planet to store
            o_velocities[0].append(outer_velocities[0][0])
            o_velocities[1].append(outer_velocities[0][1])

            # Separation vectores to store
            SI_sep.append(SI_separations[0])
            SO_sep.append(SO_separations[0])
            IO_sep.append(IO_separations[0])

            t_save = t_save + dt_store

        # Update positions
        star_positions[0] = star_positions[0] + star_velocities[0] * dt + 0.5 * star_forces[0] / m0 * dt**2
        inner_positions[0] = inner_positions[0] + inner_velocities[0] * dt + 0.5 * inner_forces[0] / m1 * dt**2
        outer_positions[0] = outer_positions[0] + outer_velocities[0] * dt + 0.5 * outer_forces[0] / m2 * dt**2

        COM[0] = (m0 * star_positions[0] + m1 * inner_positions[0] + m2 * outer_positions[0]) / (m0 + m1 + m2)

        # Reset positions to common frame of reference
        star_positions[0] = star_positions[0] - COM[0]
        inner_positions[0] = inner_positions[0] - COM[0]
        outer_positions[0] = outer_positions[0] - COM[0]

        # Calculate separation vectors
        SI_separation_vectors[0] = inner_positions[0] - star_positions[0]     # FROM star TO inner
        IS_separation_vectors[0] = - SI_separation_vectors[0]

        SO_separation_vectors[0] = outer_positions[0] - star_positions[0]     # FROM star TO outer
        OS_separation_vectors[0] = - SO_separation_vectors[0]

        IO_separation_vectors[0] = outer_positions[0] - inner_positions[0]    # FROM inner TO outer
        OI_separation_vectors[0] = - IO_separation_vectors[0]

        # Compute magnitude of the separation vectors
        SI_separations[0] = mag_vec_soft(SI_separation_vectors[0][0], SI_separation_vectors[0][1])
        SO_separations[0] = mag_vec_soft(SO_separation_vectors[0][0], SO_separation_vectors[0][1])
        IO_separations[0] = mag_vec_soft(IO_separation_vectors[0][0], IO_separation_vectors[0][1])

        # Compute new forces
        star_forces[1] = grav_force(m0, m1, IS_separation_vectors[0], SI_separations[0]) + grav_force(m0, m2, OS_separation_vectors[0], SO_separations[0])
        inner_forces[1] = grav_force(m1, m0, SI_separation_vectors[0], SI_separations[0]) + grav_force(m1, m2, OI_separation_vectors[0], IO_separations[0]) + grav_relat_force(m1, m0, c_light, SI_separation_vectors[0], SI_separations[0], inner_velocities[0]) * (m0 + m1) / m0
        outer_forces[1] = grav_force(m2, m0, SO_separation_vectors[0], SO_separations[0]) + grav_force(m2, m1, IO_separation_vectors[0], IO_separations[0]) + grav_relat_force(m1, m0, c_light, SI_separation_vectors[0], SI_separations[0], inner_velocities[0]) / m0

        # Update velocities
        star_velocities[0] = star_velocities[0] + 0.5 * (star_forces[0] + star_forces[1]) / m0 * dt
        inner_velocities[0] = inner_velocities[0] + 0.5 * (inner_forces[0] + inner_forces[1]) / m1 * dt
        outer_velocities[0] = outer_velocities[0] + 0.5 * (outer_forces[0] + outer_forces[1]) / m2 * dt

        star_forces[0] = star_forces[1]
        inner_forces[0] = inner_forces[1]
        outer_forces[0] = outer_forces[1]

        # t += dt

    ##########################################################
    # ### GET POSITIONS AND VELOCITIES FOR ALL THE BODIES ####
    ##########################################################
    x_vals_star = np.array(s_positions[0])
    y_vals_star = np.array(s_positions[1])
    x_vals_inner = np.array(i_positions[0])
    y_vals_inner = np.array(i_positions[1])
    x_vals_outer = np.array(o_positions[0])
    y_vals_outer = np.array(o_positions[1])

    vx_vals_star = np.array(s_velocities[0])
    vy_vals_star = np.array(s_velocities[1])
    vx_vals_inner = np.array(i_velocities[0])
    vy_vals_inner = np.array(i_velocities[1])
    vx_vals_outer = np.array(o_velocities[0])
    vy_vals_outer = np.array(o_velocities[1])

    time_values = np.array(time_values)

elif method_flag == "for":
    for i in trange(num_steps - 1):

        # Calculate positional change
        star_positions[i + 1] = star_positions[i] + star_velocities[i] * dt + 0.5 * star_forces[i] / m0 * dt**2
        inner_positions[i + 1] = inner_positions[i] + inner_velocities[i] * dt + 0.5 * inner_forces[i] / m1 * dt**2
        outer_positions[i + 1] = outer_positions[i] + outer_velocities[i] * dt + 0.5 * outer_forces[i] / m2 * dt**2

        COM[i + 1] = (m0 * star_positions[i + 1] + m1 * inner_positions[i + 1] + m2 * outer_positions[i + 1]) / (m0 + m1 + m2)

        # Reset positions to common frame of reference
        star_positions[i + 1] = star_positions[i + 1] - COM[i + 1]
        inner_positions[i + 1] = inner_positions[i + 1] - COM[i + 1]
        outer_positions[i + 1] = outer_positions[i + 1] - COM[i + 1]

        # Calculate separation vectors
        SI_separation_vectors[i + 1] = inner_positions[i + 1] - star_positions[i + 1]     # FROM star TO inner
        IS_separation_vectors[i + 1] = - SI_separation_vectors[i + 1]

        SO_separation_vectors[i + 1] = outer_positions[i + 1] - star_positions[i + 1]     # FROM star TO outer
        OS_separation_vectors[i + 1] = - SO_separation_vectors[i + 1]

        IO_separation_vectors[i + 1] = outer_positions[i + 1] - inner_positions[i + 1]
        OI_separation_vectors[i + 1] = - IO_separation_vectors[i + 1]

        SI_separations[i + 1] = mag_vec_soft(SI_separation_vectors[i + 1][0], SI_separation_vectors[i + 1][1])
        SO_separations[i + 1] = mag_vec_soft(SO_separation_vectors[i + 1][0], SO_separation_vectors[i + 1][1])
        IO_separations[i + 1] = mag_vec_soft(IO_separation_vectors[i + 1][0], IO_separation_vectors[i + 1][1])

        # Calculate forces
        inner_forces[i + 1] = grav_force(m1, m0, SI_separation_vectors[i + 1], SI_separations[i + 1]) + grav_force(m1, m2, OI_separation_vectors[i + 1], IO_separations[i + 1]) + (m0 + m1) / m0 * (grav_relat_force(m1, m0, c_light, SI_separation_vectors[i + 1], SI_separations[i + 1], inner_velocities[i + 1]))
        outer_forces[i + 1] = grav_force(m2, m0, SO_separation_vectors[i + 1], SO_separations[i + 1]) + grav_force(m2, m1, IO_separation_vectors[i + 1], IO_separations[i + 1]) + grav_relat_force(m1, m0, c_light, SI_separation_vectors[i + 1], SI_separations[i + 1], inner_velocities[i + 1]) / m0
        star_forces[i + 1] = grav_force(m0, m1, IS_separation_vectors[i + 1], SI_separations[i + 1]) + grav_force(m0, m2, OS_separation_vectors[i + 1], SO_separations[i + 1])

        # Calculate velocity
        star_velocities[i + 1] = star_velocities[i] + 0.5 * (star_forces[i] + star_forces[i + 1]) * dt / m0
        inner_velocities[i + 1] = inner_velocities[i] + 0.5 * (inner_forces[i] + inner_forces[i + 1]) * dt / m1
        outer_velocities[i + 1] = outer_velocities[i] + 0.5 * (outer_forces[i] + outer_forces[i + 1]) * dt / m2

    x_vals_star = star_positions[:, 0]
    y_vals_star = star_positions[:, 1]
    x_vals_inner = inner_positions[:, 0]
    y_vals_inner = inner_positions[:, 1]
    x_vals_outer = outer_positions[:, 0]
    y_vals_outer = outer_positions[:, 1]

    vx_vals_star = star_velocities[:, 0]
    vy_vals_star = star_velocities[:, 1]
    vx_vals_inner = inner_velocities[:, 0]
    vy_vals_inner = inner_velocities[:, 1]
    vx_vals_outer = outer_velocities[:, 0]
    vy_vals_outer = outer_velocities[:, 1]

else:
    exit(0)

end_time = time.time()
exec_time = (end_time - start_time) / MIN
print(f"\nSimulation finished in {exec_time:.3f} m")

#####################################################
# #### COMPUTE STATE VECTORS AND SEMIMAJOR AXES #####
#####################################################
r0 = mag_vec(x_vals_star, y_vals_star)
v0 = mag_vec(vx_vals_star, vy_vals_star)

r1 = mag_vec(x_vals_inner, y_vals_inner)
v1 = mag_vec(vx_vals_inner, vy_vals_inner)

r2 = mag_vec(x_vals_outer, y_vals_outer)
v2 = mag_vec(vx_vals_outer, vy_vals_outer)

a1 = state_vector_to_semimajor(m0, m1, r1, v1)
a2 = state_vector_to_semimajor(m0, m2, r2, v2)

##################################################
# #### COMPUTE THE TOTAL ENERGY OF THE SYSTEM ####
##################################################
kinetic_energy = 0.5 * (m0 * v0**2. + m1 * v1**2. + m2 * v2**2.)

# potential_energy = - G * m1 * m0 / mag_vec(x_vals_star - x_vals_inner, y_vals_star - y_vals_inner) - G * m2 * m0 / mag_vec(x_vals_star - x_vals_outer, y_vals_star - y_vals_outer) - G * m1 * m2 / mag_vec(x_vals_outer - x_vals_inner, y_vals_outer - y_vals_inner)
if method_flag == "while":
    potential_energy = - G * m1 * m0 / SI_sep - G * m2 * m0 / SO_sep - G * m1 * m2 / IO_sep

if method_flag == "for":
    potential_energy = - G * m1 * m0 / SI_separations - G * m2 * m0 / SO_separations - G * m1 * m2 / IO_separations

total_energy = potential_energy  # kinetic_energy #+ potential_energy
print(kinetic_energy[0], potential_energy[0])
########################################################
# #### COMPUTE THE RELATIVE ERROR OF THE SIMULATION ####
########################################################
rel_error = total_energy  # np.abs((total_energy - total_energy[0]) / total_energy[0])

# tables_dir = f"../tables/{total_t * uT / YEAR}yr_dt{dt * uT / HOUR}h/{method_flag}"
# os.makedirs(tables_dir, exist_ok=True)

# pos_tab_name = os.path.join(tables_dir, "positions.csv")
# with open(pos_tab_name, "w+") as file:
#     np.savetxt(file, np.column_stack((x_vals_star, y_vals_star, x_vals_inner, y_vals_inner, x_vals_outer, y_vals_outer)),
#                delimiter=",", fmt="%.20f",
#                header="positions")

# vel_tab_name = os.path.join(tables_dir, "velocities.csv")
# with open(vel_tab_name, "w+") as file:
#     np.savetxt(file, np.column_stack((vx_vals_star, vy_vals_star, vx_vals_inner, vy_vals_inner, vx_vals_outer, vy_vals_outer)),
#                delimiter=",", fmt="%.20f",
#                header="velocities")

# semimajor_tab_name = os.path.join(tables_dir, "semimajor_axes.csv")
# with open(semimajor_tab_name, "w+") as file:
#     np.savetxt(file, np.column_stack((a1, a2)),
#                delimiter=",", fmt="%.20f",
#                header="semimajor axis")

####################
# #### PLOTTING ####
####################
# DIRECTORY FOR SAVING IMAGES
images_dir = f"../figures/{total_t * uT / YEAR}yr_dt{dt * uT / HOUR}h/{method_flag}"
os.makedirs(images_dir, exist_ok=True)

if method_flag == "while":
    index = 1

if method_flag == "for":
    index = int(dt_store / dt)

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
plt.plot(time_values[::index] * uT / YEAR, a1[::index] * uL / AU, "b-", ms=0.4, label="Inner planet")
plt.ylabel("Semi-major axis [au]", fontsize=11)
plt.xlabel("time [yr]", fontsize=11)
plt.ylim((0.00975, 0.01175))
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
plt.plot(time_values[::index] * uT / YEAR, v1[::index], "b-", ms=0.4, label="Inner planet")
plt.ylabel("velocity [au]", fontsize=11)
plt.xlabel("time [yr]", fontsize=11)
plt.legend(loc="upper left")
fig_name = os.path.join(images_dir, "velocity_inner.png")
plt.savefig(fig_name, dpi=300)
plt.show()

plt.figure()
plt.plot(time_values[::index] * uT / YEAR, r1[::index] * uL / AU, "bo", ms=0.5, label="Inner planet")
plt.ylabel("position vector [au]", fontsize=11)
plt.xlabel("time [yr]", fontsize=11)
plt.legend(loc="upper left")
fig_name = os.path.join(images_dir, "position_vector_inner.png")
plt.xlim((0.0, total_t * uT / YEAR))
plt.ylim((0.00995, 0.01015))
plt.savefig(fig_name, dpi=300)
plt.show()

plt.figure()
plt.plot(time_values[::index] * uT / YEAR, rel_error[::index], "g-", ms=0.4, label="Relative error")
plt.ylabel("Relative error", fontsize=11)
plt.xlabel("time [yr]", fontsize=11)
plt.legend(loc="upper left")
fig_name = os.path.join(images_dir, "rel_error.png")
plt.savefig(fig_name, dpi=300)
plt.show()

# plt.figure()
# plt.plot(time_values * uT / YEAR, total_energy, "k-", label="Total energy")
# plt.legend()
# fig_name = os.path.join(images_dir, "total_energy.png")
# plt.savefig(fig_name, dpi=300)
# plt.show()
