// Astronomical constants in SI
#define AU      1.495978707e11
#define RSUN    6.957e8
#define MSUN    1.989e30
#define MJUP    1.898e27
#define RJUP    6.9911e7
#define MSAT    5.683e26
#define RSAT    6.9911e7
#define MENCEL  1.08e20
#define RENCEL  2574700
#define REARTH  6378000
#define LSUN    3.846e26
#define c_light 2.99792458e8
// Constants of time in SI
#define MIN     60.0
#define HOUR    60.0 * MIN
#define DAY     24.0 * HOUR
#define MONTH	  30.0 * DAY
#define YEAR    365.25 * DAY
#define KYEAR   1e3 * YEAR
#define MYEAR   1e6 * YEAR

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libconfig.h>
#include "cunits.h"
// #include "parameters.h"
#include "functions.h"
