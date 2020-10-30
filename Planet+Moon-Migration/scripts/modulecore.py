#!/usr/bin/env python
# -*- coding:utf-8 -*-
from util import (PI, GCONST, PLANETS, sigmasb, DP, BP, KP,
                  YEAR, GYEAR, LSUN, MSUN, RSUN)
import numpy as np

#############################################################
# SPECIFIC STELLAR ROUTINES
#############################################################


def k2Q_star_envelope(alpha, beta, epsilon, **args):
    """
      Source: Mathis, 2015
      alpha = Rc/Rp
      beta = Mc/Mp
      epsilon = Omega/Omega_crit
      args = contains behavior
    """
    gamma = alpha**3. * (1 - beta) / (beta * (1 - alpha**3.))

    line1 = 100 * PI / 63 * epsilon**2 * (alpha**5. / (1 - alpha**5.)) * (1 - gamma)**2.
    line2 = ((1 - alpha)**4.0 * (1 + 2 * alpha + 3 * alpha**2. + 1.5 * alpha**3.)**2.0
             * (1 + (1 - gamma) / gamma * alpha**3.))
    line3 = (1 + 1.5 * gamma + 2.5 / gamma * (1 + 0.5 * gamma - 1.5 * gamma**2.)
             * alpha**3. - 9. / 4. * (1 - gamma) * alpha**5.)

    k2q1 = line1 * line2 / line3**2.0

    return k2q1


#############################################################
# SPECIFIC PLANET ROUTINES
#############################################################
def k2Q_planet_envelope(alpha, beta, epsilon, **args):
    """
      Source: Mathis, 2015
      alpha = Rc/Rp
      beta = Mc/Mp
      epsilon = Omega/Omega_crit
      args = contains behavior
    """
    # if args["qk2q"]==0:return args["k2q"]

    fac0 = alpha**3.0
    fac1 = alpha**5.0
    fac2 = fac1 / (1 - fac1)

    gamma = fac0 * (1 - beta) / (beta * (1 - fac0))
    fac3 = (1 - gamma) / gamma * fac0

    k2q = 100 * PI / 63 * epsilon**2 * fac2 * (1 + fac3) / (1 + 5. / 2 * fac3)**2

    return k2q


def k2Q_planet_core(G, alpha, beta, Mp, Rp):
    gamma = alpha**3.0 * (1 - beta) / (beta * (1 - alpha**3.0))

    AA = 1.0 + 2.5 * gamma**(-1.0) * alpha**3.0 * (1.0 - gamma)
    BB = alpha**(-5.0) * (1.0 - gamma)**(-2.0)
    CC = (38.0 * PI * (alpha * Rp)**4.0) / (3.0 * GCONST * (beta * Mp)**2.0)
    DD = (2.0 / 3.0) * AA * BB * (1.0 - gamma) * (1.0 + 1.5 * gamma) - 1.5

    num = PI * G * (3.0 + 2.0 * AA)**2.0 * BB * CC
    den = DD * (6.0 * DD + 4.0 * AA * BB * CC * G)
    k2qcore = num / den
    return k2qcore


# ############RODRIGUEZ 2011########################
def S(kQ1, Mp, Ms, Rs):
    return (9 * kQ1 * Mp * Rs**5.0) / (Ms * 4.0)


def p(kQ, Mp, Ms, Rp):
    return (9 * kQ * Ms * Rp**5.0) / (Mp * 2.0)


def D(pp, SS):
    return pp / (2 * SS)
# ############RODRIGUEZ 2011########################


def Mp2Rp(Mp, Mpo, Rpo):
    c = Rpo / Mpo**0.01  # Proportionality constant
    Rp = c * Mp**0.01
    return Rp


def mloss_atmo(t, Ls, a, Mp, Rp):
    #  Zuluaga et. al (2012)
    ti = 0.06 * GYEAR * (Ls / LSUN)**-0.65

    if t < ti:
        Lx = 6.3E-4 * Ls
    else:
        Lx = 1.8928E28 * t**(-1.55)
    # Sanz-forcada et. al (2011)
    Leuv = 10**(4.8 + 0.86 * np.log10(Lx))
    k_param = 1.0  # Sanz-forcada et. al (2011)

    lxuv = (Lx + Leuv) * 1E-7
    fxuv = lxuv / (4 * PI * a**2.0)

    num = PI * Rp**3.0 * fxuv
    deno = GCONST * Mp * k_param
    return num / deno


def mloss_dragging(a, Rp, Rs, Ms, oms):
    alpha_eff = 0.3  # Zendejas et. al (2010) Venus

    return (Rp / a)**2.0 * mloss_star(Rs, Ms, oms) * alpha_eff / 2.0


def mloss_star(Rs, Ms, oms):
    smlr_sun = 1.4E-14 * MSUN / YEAR  # Zendejas et. al (2010) - smlr sun
    oms_sun = 2.67E-6
    m_loss = (smlr_sun * (Rs / RSUN)**2.0
              * (oms / oms_sun)**1.33 * (Ms / MSUN)**-3.36)

    return m_loss


def omegaCritic(M, R):
    Oc = (GCONST * M / R**3)**0.5
    return Oc


def omegadt_braking(kappa, OS, OS_saturation, osini, dobbs=False):
    if dobbs:
        gam = 1.0
        tao = GYEAR
        odt_braking = -gam / 2 * (osini / tao) * (OS / osini)**3.0
        return odt_braking

    if isinstance(OS, np.ndarray):
        odt_braking = []
        for k, o in zip(kappa, OS):
            odtb = -k * o * min(o, OS_saturation)**2.0
            # for i in range(len(k)):
            # odtb.append(-k[i] * o[i] * min(o[i], OS_saturation)**2.0)
            odt_braking.append(np.array(odtb))
        return odt_braking
    odt_braking = -kappa * OS * min(OS, OS_saturation)**2.0

    return odt_braking


def kappa_braking(OS, stellar_age, skumanich=True, alpha=0.495):

    alpha_s = 0.5  # Skumanich (1972)
    kappa = OS**-2.0 / (2.0 * stellar_age)  # Weber-Davis

    if not skumanich:
        alpha_s = alpha  # Brown et. al (2011)
        kappa = OS**(-1.0 / alpha_s) / (stellar_age / alpha_s)  # Brown (2011)
        return kappa
    return kappa


def meanMotion(a, M):
    n = (GCONST * M / a**3.0)**0.5
    return n


def eff_temp(L, Ab, a):
    Teff = (L * (1 - Ab) / (16 * PI * sigmasb * a**2.0))**(1 / 4)
    return Teff


def equilibrium_temperature(T, R, a, Ab):
    T_eq = T * (R / (2 * a))**0.5 * (1 - Ab)**0.25
    return T_eq


def luminosity(R, T):
    L = 4 * PI * R**2.0 * sigmasb * T**4.0
    return L


def alpha2beta(Mp, alpha, **args):
    beta = KP * (Mp / PLANETS.Saturn.M)**DP * alpha**BP
    return beta


def semiMajorAxis(P, M):
    a = (GCONST * M * P**2.0 / (2.0 * PI)**2.0)**(1.0 / 3.0)
    return a


def mean2axis(N, M):
    return (GCONST * M / N**2.0)**(1.0 / 3.0)


# ###################DOBS-DIXON 2004#######################
def f1e(ee):
    numer = (1 + 3.75 * ee**2.0 + 1.875 * ee**4.0 + 0.078125 * ee**6.0)
    deno = (1 - ee**2.0)**6.5
    return numer / deno


def f2e(ee):
    numer = (1 + 1.5 * ee**2.0 + 0.125 * ee**4.0)
    deno = (1 - ee**2.0)**5.0
    return numer / deno


def f3e(ee):
    numer = (1 + 7.5 * ee**2.0 + 5.625 * ee**4.0 + 0.3125 * ee**6.0)
    deno = (1 - ee**2.0)**6.0
    return numer / deno


def f4e(ee):
    numer = (1 + 3 * ee**2.0 + 0.375 * ee**4.0)
    deno = (1 - ee**2.0)**4.5
    return numer / deno


def factorbet(ee, OM, OS, N, KQ, KQ1, MP, MS, RP, RS):
    fac1 = f1e(ee) - 0.611 * f2e(ee) * (OM / N)
    fac2 = f1e(ee) - 0.611 * f2e(ee) * (OS / N)
    lamb = (KQ / KQ1) * (MS / MP)**2.0 * (RP / RS)**5.0
    return 18.0 / 7.0 * (fac1 + fac2 / lamb)


def power(ee, aa, KQ, Ms, Rp):
    keys = (GCONST * Ms)**1.5 * ((2 * Ms * Rp**5.0 * ee**2.0 * KQ) / 3)
    coeff = 15.75 * aa**(-7.5)
    return coeff * keys
# ###################DOBS-DIXON 2004#######################


def mag_vec(a, b, c):
    mag = (a**2. + b**2. + c**2.)**0.5

    return mag


def state_vector_to_semimajor(M, m, r, v):
    mu = GCONST * (M + m)
    eps = v**2. / 2 - mu / r
    a = -mu / (2. * eps)
    # a = 1 / (2 / r - v**2. / mu)
    return a


def escape_velocity(M, r):
    ev = (2 * GCONST * M / r)**0.5

    return ev


def polar_mom_iner(M, R):
    return 2 * M * R**2. / 5.


#############################################################
# DIFFERENTIAL EQUATIONS
#############################################################
def dr1_dt(y, t, parameters):
    """
    Integration of the governing vector differential equation.
    d2r_dt2 = -(mu/R^3)*r with d2r_dt2 and r as vectors.
    Initial position and velocity are given.
    y[0:2] = position components
    y[3:] = velocity components
    """
    M1 = parameters["M1"]
    R1 = parameters["R1"]
    R2 = parameters["R2"]
    M2 = parameters["M2"]
    Ms = parameters["Ms"]
    c_light = parameters["c_light"]
    mu1 = GCONST * (Ms + M1)

    alpha_1 = parameters["alpha_1"]
    beta_1 = parameters["beta_1"]
    rigidity_1 = parameters["rigidity_1"]

    alpha_2 = parameters["alpha_2"]
    beta_2 = parameters["beta_2"]
    rigidity_2 = parameters["rigidity_2"]

    # Variables of body 2
    q2_x = parameters["q2_x"]
    q2_y = parameters["q2_y"]
    q2_z = parameters["q2_z"]
    q2_vx = parameters["q2_vx"]
    q2_vy = parameters["q2_vy"]
    q2_vz = parameters["q2_vz"]

    om1_x = parameters["om1_x"]
    om1_y = parameters["om1_y"]
    om1_z = parameters["om1_z"]

    om2_x = parameters["om2_x"]
    om2_y = parameters["om2_y"]
    om2_z = parameters["om2_z"]

    # Secondary properties body 1
    eps1 = mag_vec(om1_x, om1_y, om1_z) / omegaCritic(M1, R1)

    if parameters["key"] == 0:
        k2q1 = parameters["k2q1"]
    else:
        k2q1_e = k2Q_planet_envelope(alpha_1,
                                     beta_1,
                                     eps1)
        k2q1_c = k2Q_planet_core(rigidity_1,
                                 alpha_1,
                                 beta_1, M1, R1)
        k2q1 = k2q1_e + k2q1_c

    # Secondary properties body 2
    eps2 = mag_vec(om2_x, om2_y, om2_z) / omegaCritic(M2, R2)

    if parameters["key"] == 0:
        k2q2 = parameters["k2q2"]
    else:
        k2q2_e = k2Q_planet_envelope(alpha_2,
                                     beta_2,
                                     eps2)
        k2q2_c = k2Q_planet_core(rigidity_2,
                                 alpha_2,
                                 beta_2, M2, R2)
        k2q2 = k2q2_e + k2q2_c

    # Position vectors magnitude
    r1 = mag_vec(y[0], y[1], y[2])
    v1 = mag_vec(y[3], y[4], y[5])
    r2 = mag_vec(q2_x, q2_y, q2_z)
    v2 = mag_vec(q2_vx, q2_vy, q2_vz)

    a1 = state_vector_to_semimajor(Ms, M1, r1, v1)
    n1 = meanMotion(a1, Ms)

    a2 = state_vector_to_semimajor(Ms, M2, r2, v2)
    n2 = meanMotion(a2, Ms)

    # Force of central body on body 1
    f_01_x = -(mu1 / r1**3) * y[0]
    f_01_y = -(mu1 / r1**3) * y[1]
    f_01_z = -(mu1 / r1**3) * y[2]

    # Force of body 2 on body 1
    f_21_x = GCONST * M2 * ((q2_x - y[0]) / mag_vec((q2_x - y[0]), (q2_y - y[1]), (q2_z - y[2]))**3. - q2_x / r2**3.)
    f_21_y = GCONST * M2 * ((q2_y - y[1]) / mag_vec((q2_x - y[0]), (q2_y - y[1]), (q2_z - y[2]))**3. - q2_y / r2**3.)
    f_21_z = GCONST * M2 * ((q2_z - y[2]) / mag_vec((q2_x - y[0]), (q2_y - y[1]), (q2_z - y[2]))**3. - q2_z / r2**3.)

    # Relativity correction
    dot_product = y[0] * y[3] + y[1] * y[4] + y[2] * y[5]
    f_rel_x = GCONST * M1 * Ms / (c_light**2 * r1**3.) * ((4 * GCONST * Ms / r1 - v1**2.) * y[0] + 4 * dot_product * y[3])
    f_rel_y = GCONST * M1 * Ms / (c_light**2 * r1**3.) * ((4 * GCONST * Ms / r1 - v1**2.) * y[1] + 4 * dot_product * y[4])
    f_rel_z = GCONST * M1 * Ms / (c_light**2 * r1**3.) * ((4 * GCONST * Ms / r1 - v1**2.) * y[2] + 4 * dot_product * y[5])

    # Tidal force on body 1
    cross_pr_r_om1 = [y[1] * om1_z - y[2] * om1_y,
                      y[2] * om1_x - y[0] * om1_z,
                      y[0] * om1_y - y[1] * om1_x]

    f_tidal_x = -3 * k2q1 * GCONST * Ms**.2 * R1**5. / (n1 * r1**10.) * (2 * y[0] * dot_product + r1**2. * (cross_pr_r_om1[0] + y[3]))
    f_tidal_y = -3 * k2q1 * GCONST * Ms**.2 * R1**5. / (n1 * r1**10.) * (2 * y[1] * dot_product + r1**2. * (cross_pr_r_om1[1] + y[4]))
    f_tidal_z = -3 * k2q1 * GCONST * Ms**.2 * R1**5. / (n1 * r1**10.) * (2 * y[2] * dot_product + r1**2. * (cross_pr_r_om1[2] + y[5]))

    # Tidal force on body 2
    cross_pr_r_om2 = [q2_y * om2_z - q2_z * om2_y,
                      q2_z * om2_x - q2_x * om2_z,
                      q2_x * om2_y - q2_y * om2_x]

    g_tidal_x = -3 * k2q2 * GCONST * Ms**.2 * R2**5. / (n2 * r2**10.) * (2 * q2_x * dot_product + r2**2. * (cross_pr_r_om2[0] + q2_vx))
    g_tidal_y = -3 * k2q2 * GCONST * Ms**.2 * R2**5. / (n2 * r2**10.) * (2 * q2_y * dot_product + r2**2. * (cross_pr_r_om2[1] + q2_vy))
    g_tidal_z = -3 * k2q2 * GCONST * Ms**.2 * R2**5. / (n2 * r2**10.) * (2 * q2_z * dot_product + r2**2. * (cross_pr_r_om2[2] + q2_vz))

    dy0 = y[3]
    dy1 = y[4]
    dy2 = y[5]
    dy3 = f_01_x + f_21_x + (Ms + M1) * (f_tidal_x + f_rel_x) / (Ms * M1) + g_tidal_x / Ms
    dy4 = f_01_y + f_21_y + (Ms + M1) * (f_tidal_y + f_rel_y) / (Ms * M1) + g_tidal_y / Ms
    dy5 = f_01_z + f_21_z + (Ms + M1) * (f_tidal_z + f_rel_z) / (Ms * M1) + g_tidal_z / Ms

    return [dy0, dy1, dy2, dy3, dy4, dy5]


def dr2_dt(y, t, parameters):
    """
    Integration of the governing vector differential equation.
    d2r_dt2 = -(mu/R^3)*r with d2r_dt2 and r as vectors.
    Initial position and velocity are given.
    y[0:2] = position components
    y[3:] = velocity components
    """
    M1 = parameters["M1"]
    R1 = parameters["R1"]
    M2 = parameters["M2"]
    R2 = parameters["R2"]
    Ms = parameters["Ms"]
    c_light = parameters["c_light"]
    mu2 = GCONST * (Ms + M2)

    alpha_1 = parameters["alpha_1"]
    beta_1 = parameters["beta_1"]
    rigidity_1 = parameters["rigidity_1"]

    alpha_2 = parameters["alpha_2"]
    beta_2 = parameters["beta_2"]
    rigidity_2 = parameters["rigidity_2"]

    # Variables of body 1
    q1_x = parameters["q1_x"]
    q1_y = parameters["q1_y"]
    q1_z = parameters["q1_z"]
    q1_vx = parameters["q1_vx"]
    q1_vy = parameters["q1_vy"]
    q1_vz = parameters["q1_vz"]

    om1_x = parameters["om1_x"]
    om1_y = parameters["om1_y"]
    om1_z = parameters["om1_z"]

    om2_x = parameters["om2_x"]
    om2_y = parameters["om2_y"]
    om2_z = parameters["om2_z"]

    # Secondary properties body 1
    eps1 = mag_vec(om1_x, om1_y, om1_z) / omegaCritic(M1, R1)

    if parameters["key"] == 0:
        k2q1 = parameters["k2q2"]
    else:
        k2q1_e = k2Q_planet_envelope(alpha_1,
                                     beta_1,
                                     eps1)
        k2q1_c = k2Q_planet_core(rigidity_1,
                                 alpha_1,
                                 beta_1, M1, R1)
        k2q1 = k2q1_e + k2q1_c

    # Secondary properties body 2
    eps2 = mag_vec(om2_x, om2_y, om2_z) / omegaCritic(M2, R2)

    if parameters["key"] == 0:
        k2q2 = parameters["k2q2"]
    else:
        k2q2_e = k2Q_planet_envelope(alpha_2,
                                     beta_2,
                                     eps2)
        k2q2_c = k2Q_planet_core(rigidity_2,
                                 alpha_2,
                                 beta_2, M2, R2)
        k2q2 = k2q2_e + k2q2_c

    # Position vectors magnitude
    r1 = mag_vec(q1_x, q1_y, q1_z)
    v1 = mag_vec(q1_vx, q1_vy, q1_vz)
    r2 = mag_vec(y[0], y[1], y[2])
    v2 = mag_vec(y[3], y[4], y[5])

    a1 = state_vector_to_semimajor(Ms, M1, r1, v1)
    n1 = meanMotion(a1, Ms)

    a2 = state_vector_to_semimajor(Ms, M2, r2, v2)
    n2 = meanMotion(a2, Ms)

    # Force of central body on body 2
    f_02_x = -(mu2 / r2**3) * y[0]
    f_02_y = -(mu2 / r2**3) * y[1]
    f_02_z = -(mu2 / r2**3) * y[2]

    # Force of body 1 on body 2
    f_12_x = GCONST * M1 * ((q1_x - y[0]) / mag_vec((q1_x - y[0]), (q1_y - y[1]), (q1_z - y[2]))**3. - q1_x / r1**3.)
    f_12_y = GCONST * M1 * ((q1_y - y[1]) / mag_vec((q1_x - y[0]), (q1_y - y[1]), (q1_z - y[2]))**3. - q1_y / r1**3.)
    f_12_z = GCONST * M1 * ((q1_z - y[2]) / mag_vec((q1_x - y[0]), (q1_y - y[1]), (q1_z - y[2]))**3. - q1_z / r1**3.)

    # Relativity correction
    dot_product = q1_x * q1_vx + q1_y * q1_vy + q1_z * q1_vz
    f_rel_x = GCONST * M1 * Ms / (c_light**2 * r1**3.) * ((4 * GCONST * Ms / r1 - v1**2.) * q1_x + 4 * dot_product * q1_vx)
    f_rel_y = GCONST * M1 * Ms / (c_light**2 * r1**3.) * ((4 * GCONST * Ms / r1 - v1**2.) * q1_y + 4 * dot_product * q1_vy)
    f_rel_z = GCONST * M1 * Ms / (c_light**2 * r1**3.) * ((4 * GCONST * Ms / r1 - v1**2.) * q1_z + 4 * dot_product * q1_vz)

    # Tidal force on body 1
    cross_pr_r_om1 = [q1_y * om1_z - q1_z * om1_y,
                      q1_z * om1_x - q1_x * om1_z,
                      q1_x * om1_y - q1_y * om1_x]

    f_tidal_x = -3 * k2q1 * GCONST * Ms**.2 * R1**5. / (n1 * r1**10.) * (2 * q1_x * dot_product + r1**2. * (cross_pr_r_om1[0] + q1_vx))
    f_tidal_y = -3 * k2q1 * GCONST * Ms**.2 * R1**5. / (n1 * r1**10.) * (2 * q1_y * dot_product + r1**2. * (cross_pr_r_om1[1] + q1_vy))
    f_tidal_z = -3 * k2q1 * GCONST * Ms**.2 * R1**5. / (n1 * r1**10.) * (2 * q1_z * dot_product + r1**2. * (cross_pr_r_om1[2] + q1_vz))

    # Tidal force on body 2
    cross_pr_r_om2 = [y[1] * om2_z - y[2] * om2_y,
                      y[2] * om2_x - y[0] * om2_z,
                      y[0] * om2_y - y[1] * om2_x]

    g_tidal_x = -3 * k2q2 * GCONST * Ms**.2 * R2**5. / (n2 * r2**10.) * (2 * y[0] * dot_product + r2**2. * (cross_pr_r_om2[0] + y[3]))
    g_tidal_y = -3 * k2q2 * GCONST * Ms**.2 * R2**5. / (n2 * r2**10.) * (2 * y[1] * dot_product + r2**2. * (cross_pr_r_om2[1] + y[4]))
    g_tidal_z = -3 * k2q2 * GCONST * Ms**.2 * R2**5. / (n2 * r2**10.) * (2 * y[2] * dot_product + r2**2. * (cross_pr_r_om2[2] + y[5]))

    dy0 = y[3]
    dy1 = y[4]
    dy2 = y[5]
    dy3 = f_02_x + f_12_x + (Ms + M2) * g_tidal_x / (Ms * M2) + (f_rel_x + f_tidal_x) / Ms
    dy4 = f_02_y + f_12_y + (Ms + M2) * g_tidal_y / (Ms * M2) + (f_rel_y + f_tidal_y) / Ms
    dy5 = f_02_z + f_12_z + (Ms + M2) * g_tidal_z / (Ms * M2) + (f_rel_z + f_tidal_z) / Ms

    return [dy0, dy1, dy2, dy3, dy4, dy5]


def dom1_dt(y, t, parameters):
    """
    Integration of the angular velocity differential equation.
    """
    Ms = parameters["Ms"]
    M1 = parameters["M1"]
    R1 = parameters["R1"]

    alpha_1 = parameters["alpha_1"]
    beta_1 = parameters["beta_1"]
    rigidity_1 = parameters["rigidity_1"]

    om1_x = y[0]
    om1_y = y[1]
    om1_z = y[2]

    # Variables of body 1
    q1_x = parameters["q1_x"]
    q1_y = parameters["q1_y"]
    q1_z = parameters["q1_z"]
    q1_vx = parameters["q1_vx"]
    q1_vy = parameters["q1_vy"]
    q1_vz = parameters["q1_vz"]
    r1 = mag_vec(q1_x, q1_y, q1_z)
    v1 = mag_vec(q1_vx, q1_vy, q1_vz)

    a1 = state_vector_to_semimajor(Ms, M1, r1, v1)
    n1 = meanMotion(a1, Ms)
    mom1 = polar_mom_iner(M1, R1)

    # Secondary properties planet
    eps1 = mag_vec(om1_x, om1_y, om1_z) / omegaCritic(M1, R1)
    if parameters["key"] == 0:
        k2q1 = parameters["k2q1"]
    else:
        k2q1_e = k2Q_planet_envelope(alpha_1,
                                     beta_1,
                                     eps1)
        k2q1_c = k2Q_planet_core(rigidity_1,
                                 alpha_1,
                                 beta_1, M1, R1)
        k2q1 = k2q1_e + k2q1_c

    dot_pr_r_om = q1_x * om1_x + q1_y * om1_y + q1_z * om1_z
    cross_pr_r_v = [q1_y * q1_vz - q1_z * q1_vy,
                    q1_z * q1_vx - q1_x * q1_vz,
                    q1_x * q1_vy - q1_y * q1_vx]

    dom1dt0 = -3 * k2q1 * GCONST * Ms**.2 * R1**5. / (mom1 * n1 * r1**8.) * (dot_pr_r_om * q1_x - r1**2. * om1_x + cross_pr_r_v[0])
    dom1dt1 = -3 * k2q1 * GCONST * Ms**.2 * R1**5. / (mom1 * n1 * r1**8.) * (dot_pr_r_om * q1_y - r1**2. * om1_y + cross_pr_r_v[1])
    dom1dt2 = -3 * k2q1 * GCONST * Ms**.2 * R1**5. / (mom1 * n1 * r1**8.) * (dot_pr_r_om * q1_z - r1**2. * om1_z + cross_pr_r_v[2])

    return [dom1dt0, dom1dt1, dom1dt2]


def dom2_dt(y, t, parameters):
    """
    Integration of the angular velocity differential equation.
    """
    Ms = parameters["Ms"]
    M2 = parameters["M2"]
    R2 = parameters["R2"]

    alpha_2 = parameters["alpha_2"]
    beta_2 = parameters["beta_2"]
    rigidity_2 = parameters["rigidity_2"]

    om2_x = y[0]
    om2_y = y[1]
    om2_z = y[2]

    # Variables of body 1
    q2_x = parameters["q2_x"]
    q2_y = parameters["q2_y"]
    q2_z = parameters["q2_z"]
    q2_vx = parameters["q2_vx"]
    q2_vy = parameters["q2_vy"]
    q2_vz = parameters["q2_vz"]
    r2 = mag_vec(q2_x, q2_y, q2_z)
    v2 = mag_vec(q2_vx, q2_vy, q2_vz)

    a2 = state_vector_to_semimajor(Ms, M2, r2, v2)
    n2 = meanMotion(a2, Ms)

    mom2 = polar_mom_iner(M2, R2)
    # Secondary properties planet
    eps2 = mag_vec(om2_x, om2_y, om2_z) / omegaCritic(M2, R2)

    if parameters["key"] == 0:
        k2q2 = parameters["k2q2"]
    else:
        k2q2_e = k2Q_planet_envelope(alpha_2,
                                     beta_2,
                                     eps2)
        k2q2_c = k2Q_planet_core(rigidity_2,
                                 alpha_2,
                                 beta_2, M2, R2)
        k2q2 = k2q2_e + k2q2_c

    dot_pr_r_om = q2_x * om2_x + q2_y * om2_y + q2_z * om2_z
    cross_pr_r_v = [q2_y * q2_vz - q2_z * q2_vy,
                    q2_z * q2_vx - q2_x * q2_vz,
                    q2_x * q2_vy - q2_y * q2_vx]

    dom2dt0 = -3 * k2q2 * GCONST * Ms**.2 * R2**5. / (mom2 * n2 * r2**8.) * (dot_pr_r_om * q2_x - r2**2. * om2_x + cross_pr_r_v[0])
    dom2dt1 = -3 * k2q2 * GCONST * Ms**.2 * R2**5. / (mom2 * n2 * r2**8.) * (dot_pr_r_om * q2_y - r2**2. * om2_y + cross_pr_r_v[1])
    dom2dt2 = -3 * k2q2 * GCONST * Ms**.2 * R2**5. / (mom2 * n2 * r2**8.) * (dot_pr_r_om * q2_z - r2**2. * om2_z + cross_pr_r_v[2])

    return [dom2dt0, dom2dt1, dom2dt2]


#############################################################
# INTEGRATION OF THE WHOLE SYSTEM
#############################################################
def global_differential_equation(q, t, parameters):

    q1_x = q[0]
    q1_y = q[1]
    q1_z = q[2]
    q1_vx = q[3]
    q1_vy = q[4]
    q1_vz = q[5]
    om1_x = q[6]
    om1_y = q[7]
    om1_z = q[8]

    q2_x = q[9]
    q2_y = q[10]
    q2_z = q[11]
    q2_vx = q[12]
    q2_vy = q[13]
    q2_vz = q[14]
    om2_x = q[15]
    om2_y = q[16]
    om2_z = q[17]

    parameters["q1_x"] = q1_x
    parameters["q1_y"] = q1_y
    parameters["q1_z"] = q1_z
    parameters["q1_vx"] = q1_vx
    parameters["q1_vy"] = q1_vy
    parameters["q1_vz"] = q1_vz
    parameters["om1_x"] = om1_x
    parameters["om1_y"] = om1_y
    parameters["om1_z"] = om1_z

    parameters["q2_x"] = q2_x
    parameters["q2_y"] = q2_y
    parameters["q2_z"] = q2_z
    parameters["q2_vx"] = q2_vx
    parameters["q2_vy"] = q2_vy
    parameters["q2_vz"] = q2_vz
    parameters["om2_x"] = om2_x
    parameters["om2_y"] = om2_y
    parameters["om2_z"] = om2_z

    dr1dt_sol = dr1_dt([q1_x, q1_y, q1_z, q1_vx, q1_vy, q1_vz],
                       t, parameters)

    dr2dt_sol = dr2_dt([q2_x, q2_y, q2_z, q2_vx, q2_vy, q2_vz],
                       t, parameters)

    dom1dt_sol = dom1_dt([om1_x, om1_y, om1_z], t, parameters)

    dom2dt_sol = dom2_dt([om2_x, om2_y, om2_z], t, parameters)

    return dr1dt_sol + dom1dt_sol + dr2dt_sol + dom2dt_sol
