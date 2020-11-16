# !/usr/bin/python3
#  -*- coding: utf-8 -*-
from util import (PI, G, PLANETS, sigmasb, DP, BP, KP,
                  YEAR, GYEAR, LSUN, MSUN, RSUN)
import numpy as np


def mag_vec(a, b):
    mag = (a**2. + b**2.)**0.5

    return mag


def mag_vec_soft(a, b, soft=0.5):
    mag = (a**2. + b**2. + soft**2.)**0.5

    return mag


def state_vector_to_semimajor(M, m, r, v):
    mu = G * (M + m)
    eps = v**2. / 2 - mu / r
    a = -mu / (2. * eps)
    # a = 1 / (2 / r - v**2. / mu)
    return a


def escape_velocity(M, r):
    ev = (2 * G * M / r)**0.5

    return ev


def orbital_velocity(M, r):
    ov = (G * M / r)**0.5

    return ov


def initial_velocity(M, m, r0, a):
    v0 = ((M + m) * (2. / r0 - 1. / a))**0.5

    return v0


def polar_mom_iner(M, R):
    return 2 * M * R**2. / 5.


def grav_force(m, M, sep_vec, pos):
    g_force = - G * m * M * sep_vec / pos**3

    return g_force


def grav_relat_force(m, M, c, sep_vec, pos, vel_vec):
    coeff = G * m * M / (c**2. * pos**3.)
    first_term = (4. * G * M / pos - np.dot(vel_vec, vel_vec)) * sep_vec
    second_term = 4. * np.dot(sep_vec, vel_vec) * vel_vec

    return coeff * (first_term + second_term)


def new_positions(x, v, f, m, t):
    new_p = x + v * t + 0.5 * f / m * t**2

    return new_p
