"""
For physical constants and tolerances.
"""

import numpy as np
from scipy.constants import (
    speed_of_light,
    mu_0,
    epsilon_0,
    elementary_charge,
    hbar,
    electron_mass,
    physical_constants,
    angstrom,
    Boltzmann,
)

# float epsilon.
M_ZERO = np.finfo(float).eps

# e = 1.602176634e-19 [C]
elem_charge_SI = elementary_charge

# m_e = 9.1093837015e-31 [kg]
elec_mass_SI = electron_mass

# hbar = 1.0545718176461565e-34 [J * s]
hbar_SI = hbar

# k_B = 1.380649e-23 [J / K]
k_B_SI = Boltzmann

# mu_B = 9.274010078362164e-24 [J / T]
bohr_magn_SI = elem_charge_SI * hbar_SI / (2 * elec_mass_SI)

# epsilon_0 = 8.8541878128e-12  [F / m]
eps0_SI = epsilon_0

# mu_0 = 1.25663706212e-06 [H / m]
mu0_SI = mu_0

# c = 299792458.0 [m / s]
speedlight_SI = speed_of_light

# Electron Volt in atomic units [E_h / eV]
eV_au = physical_constants["electron volt-hartree relationship"][0]

# Electron Volt in seconds []
eV_seconds = 6.582119e-16

# Rydberg constant times hc in eV = 13.605693122994 [eV]
Ry_eV = physical_constants["Rydberg constant times hc in eV"][0]

# Ang = 1e-10 [m]
Ang_SI = angstrom

# Bohr to Ang
bohr_angstrom_internal = physical_constants["Bohr radius"][0] / angstrom
bohr = physical_constants["Bohr radius"][0] / angstrom

# joule to electron volt [eV / J]
joul_to_eV = physical_constants["joule-electron volt relationship"][0]
