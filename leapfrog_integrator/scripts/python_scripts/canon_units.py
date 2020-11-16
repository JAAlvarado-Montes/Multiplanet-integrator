#!/usr/bin/env python

############################################################
# Calcula sistema de unidades canonicas

G = 6.6740831e-11


def units(**kwargs):
    if 'uM' in kwargs.keys() and 'uL' in kwargs.keys():
        uT = (kwargs['uL']**3 / (G * kwargs['uM']))**0.5
        return [kwargs.get('uM'), kwargs.get('uL'), uT]

    elif 'uM' in kwargs.keys() and 'uT' in kwargs.keys():
        uL = (G * kwargs['uM'] * kwargs['uT']**2)**(1.0 / 3.0)
        return [kwargs.get('uM'), uL, kwargs.get('uT')]

    elif 'uL' in kwargs.keys() and 'uT' in kwargs.keys():
        uM = (kwargs['uL']**3 / (kwargs['uT']**2 * G))
        return [uM, kwargs.get('uL'), kwargs.get('uT')]
