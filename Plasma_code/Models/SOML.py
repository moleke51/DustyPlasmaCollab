import numpy as np
import scipy.special as sps
from scipy.optimize import bisect


def get_name():
    return "SOML"


def colour():
    return "violet"


def get_info():
    assumptions_list = [
        "Model assumptions:\n",
        "Spherical symmetry\n",
        "No collisions\n",
        "No magnetic field\n",
        "No external electric field\n",
        "No electron emission of any kind\n",
        "Quasi-neutrality in bulk plasma\n",
        "Conservation of particle energy\n",
        "Conservation of particle angular momentum\n",
        "Limiting trajectory is the grazing incidence\n",
        "Flowing Maxwell-Boltzmann velocity distribution for ions\n",
    ]
    validity_list = [
        "Validity:\n",
        "Flowing plasma\n",
        "Small dust\n",
        "Any ion temperature\n",
    ]
    reference_list = [
        "References:\n",
        "C. T. N. Willis, “Dust in stationary and flowing plasmas,” Physics PhD Thesis, Imperial College London, March 2012\n",
        "D. M. Thomas, “Theory and simulation of the charging of dust in plasmas,” Physics PhD Thesis, Imperial College London, March 2016\n",
    ]
    string = (
        " ".join(assumptions_list) + " ".join(validity_list) + " ".join(reference_list)
    )
    return string


# SOML (Shifted OML) model for normalised dust surface potential - eqn 2.138 in Thomas' thesis
# Define SOML equation to solve
def SOML_function(Phi, Theta, mu, z, alpha, upsilon):  # gamma = 3 for flowing plasmas
    s_1 = ((np.sqrt(np.pi)) * (1 + 2 * (upsilon ** 2)) * sps.erf(upsilon)) / (
        4 * upsilon
    ) + 0.5 * np.exp(-(upsilon ** 2))
    s_2 = (np.sqrt(np.pi) * sps.erf(upsilon)) / (2 * upsilon)
    return (np.sqrt(Theta) / mu) * (s_1 - (s_2 * Phi) / Theta) - np.exp(Phi)


def potential_finder(dictionarylist):
    for _vardict in dictionarylist:
        if _vardict.get("Norm_var_name") is not None:
            if _vardict.get("Norm_var_name") == "alpha":
                alpha = _vardict.get("Norm_value")
            elif _vardict.get("Norm_var_name") == "z":
                z = _vardict.get("Norm_value")
            elif _vardict.get("Norm_var_name") == "mu":
                mu = _vardict.get("Norm_value")
            elif _vardict.get("Norm_var_name") == "upsilon":
                upsilon = _vardict.get("Norm_value")
            elif _vardict.get("Norm_var_name") == "Theta":
                Theta = _vardict.get("Norm_value")

    Phi = bisect(SOML_function, -10, 10, args=(Theta, mu, z, alpha, upsilon))
    return np.absolute(Phi)


def priority(dictionarylist):
    P = 0
    for _vardict in dictionarylist:
        if _vardict.get("Norm_var_name") is not None:
            if _vardict.get("Norm_var_name") == "alpha":
                alpha = _vardict.get("Norm_value")
            elif _vardict.get("Norm_var_name") == "z":
                P += 1
            elif _vardict.get("Norm_var_name") == "mu":
                P += 1
            elif _vardict.get("Norm_var_name") == "upsilon":
                upsilon = _vardict.get("Norm_value")
                if upsilon == 0:
                    return 0
                P += 1
            elif _vardict.get("Norm_var_name") == "Theta":
                Theta = _vardict.get("Norm_value")
                if Theta >= 1e-4:
                    P += 1
                else:
                    P += 0.5
            else:
                if _vardict.get("Norm_value") != _vardict.get("default value"):
                    return 0

    if alpha > 1.25 * Theta ** (0.4):
        return 0
    P += 1

    return P
