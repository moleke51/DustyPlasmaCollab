import numpy as np
import SOML as soml
import SMOML as smoml


def get_name():
    return "FNF"


def colour():
    return "purple"


def get_info():
    assumptions_list = [
        "Model assumptions:\n",
        "Spherical symmetry\n",
        "No collisions\n",
        "No magnetic field\n",
        "No external electric field\n",
        "No electron emission of any kind\n",
        "Quasi-neutrality in bulk plasma\n",
    ]
    validity_list = [
        "Validity:\n",
        "Flowing plasma\n",
        "Intermediate sized dust\n",
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


# Linear model for normalised dust surface potential - eqn 4.3 in Willis' thesis
def Linear_function(phi_SOML, phi_TS, alpha_SOML, alpha_TS, alpha):
    x = ((phi_TS - phi_SOML) / (np.log(alpha_TS) - np.log(alpha_SOML))) * np.log(
        (alpha) / (alpha_TS)
    ) + phi_TS
    return x


def potential_finder(dictionarylist):  # gamma = 3 for flowing plasmas

    for _vardict in dictionarylist:
        if _vardict.get("Norm_var_name") is not None:
            if _vardict.get("Norm_var_name") == "alpha":
                alpha = _vardict.get("Norm_value")
            elif _vardict.get("Norm_var_name") == "Theta":
                Theta = _vardict.get("Norm_value")
    alpha_OML = 1.25 * (Theta) ** 0.4  # Assume this is the same as the static case
    alpha_TS = 100
    Phi_SMOML = smoml.potential_finder(dictionarylist)
    Phi_SOML = soml.potential_finder(dictionarylist)
    Phi = Linear_function(Phi_SOML, Phi_SMOML, alpha_OML, alpha_TS, alpha)
    return Phi  # returned phi is positive


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
                P += 1
    if alpha < 100 and alpha > 1.25 * Theta ** (0.4):
        P += 1
    return P
