#!/usr/bin/env python
# -*- coding:utf-8 -*-
#############################################################
# EXTERNAL PACKAGES
#############################################################
from scipy import constants as const
from collections import namedtuple

import astropy.constants as c
import numpy as np

#############################################################
# CORE ROUTINES
#############################################################


class dict2obj(object):
    def __init__(self, dic={}):
        self.__dict__.update(dic)

    def __add__(self, other):
        for attr in other.__dict__.keys():
            exec(f"self.{attr}=other.{attr}")
        return self

#############################################################
# CONSTANTS
#############################################################


# NUMERICAL CONSTANTS
PI = np.pi
RAD = 1 / const.degree

# PHYSICAL CONSTANTS
MIN = const.minute  # s
HOUR = const.hour  # s
DAY = const.day  # s
YEAR = const.Julian_year  # s
KYEAR = 1e3 * YEAR
MYEAR = 1e6 * YEAR  # s
GYEAR = 1e9 * YEAR  # s
GCONST = 1.0  # 6.6740831e-11  # m^3 / kg s^2
sigmasb = const.Stefan_Boltzmann
c_light = const.c  # m / s

# ASTRONOMICAL CONSTANTS
AU = const.au  # m
MSUN = c.M_sun.value  # kg
RSUN = 6.957e8  # m
MENCEL = 1.08e20  # kg
RENCEL = 2574700  # m
MEARTH = c.M_earth.value  # kg
REARTH = 6378000  # m
LSUN = 3.846e26  # W

PLANETS = dict2obj(dict(
                   Jupiter=dict2obj(dict(
                                    M=1.898e27,  # kg
                                    R=6.9911e7,  # m
                                    P=29.5 * YEAR,  # s
                                    Prot=9.4 * HOUR,  # s
                                    alpha=0.126,
                                    beta=0.020)),
                   Saturn=dict2obj(dict(
                                   M=5.683e26,  # kg
                                   R=6.0268e7,  # m
                                   P=10.8 * YEAR,  # s
                                   Prot=10.656 * HOUR,  # s
                                   alpha=0.219,
                                   beta=0.196)),
                   Uranus=dict2obj(dict(
                                   M=86.8e24,  # kg
                                   R=2.5632e7,  # m
                                   P=84 * YEAR,  # s
                                   Prot=17.24 * HOUR,  # s
                                   alpha=0.30,
                                   beta=0.093)),
                   Neptune=dict2obj(dict(
                                    M=1.024e26,  # kg
                                    R=2.4622e7,  # m
                                    P=164.8 * YEAR,  # s
                                    Prot=16.11 * HOUR,  # s
                                    alpha=0.35,
                                    beta=0.131))
                   ))


# ############################################################
# Util Functions
# ############################################################

def fmt(x, pos):
    """Writes scientific notatino formater for plots

    Returns:
        string: Scientific notation of input string
    """
    a, b = f"{x:.1e}".split("e")
    b = int(b)
    return rf'${a} \times 10^{{{b}}}$'


def fmt_2(x):
    """Writes scientific notatino formater for plots

    Returns:
        string: Scientific notation of input string
    """
    a, b = f"{x:.2e}".split("e")
    b = int(b)
    return rf'{a} \times 10^{{{b}}}'


def create_meshgrid(minpar_1, maxpar_1,
                    minpar_2, maxpar_2,
                    grid_division=500):
    """Creates 2D meshgrid to be used in countour plots.

    Args:
        minpar_1 (float): Minimum for first parameter.
        maxpar_1 (float): Maximum for first parameter.
        minpar_2 (float): Minimum for second parameter.
        maxpar_2 (float): Maximum for second parameter.
        grid_division (int, optional): Description

    Returns:
        Meshgrid of two different parameter to be used for contours or 2D
        histograms.

    """
    # VECTORS FOR PARAMETERS
    first_parameter = np.linspace(minpar_1, maxpar_1, grid_division)
    second_parameter = np.linspace(minpar_2, maxpar_2, grid_division)

    first_par, second_par = np.meshgrid(first_parameter, second_parameter)

    Outputs = namedtuple("Outputs", "first_par second_par")

    return Outputs(first_par, second_par)


# ############################################################
# FIT CONSTANTS
# ############################################################
# Exponent of radius mass-scaling law
ALPHAR = 0.156
# Coefficients of radius time-scaling law
A = 3.964
B = -0.064
C = 3.364
t0 = 1e8

# Scaling constants for alpha and beta
KP = 90.0742985384
BP = 4.0
DP = -0.232

#############################################################
# GYRATION RADIUS
#############################################################
GR = 0.2  # None
K2QVE = 3E-5
