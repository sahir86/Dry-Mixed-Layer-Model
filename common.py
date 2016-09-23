"""
Leif Denby's attribute dictionary class that allows you to use the dot notation to access dictionary entries.
Must use as variable = AttrDict(constants).Then can access using variable.g for example
""" 

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

# Library of relevant constants

constants = {
    "Cp": 1005.46, # J/Kg.K
    "g": 9.80665, # m/s^2
    "D": 4e-6, #1/s
    "delta_f": 40.0, #W/m^2
    "delta_s": 12.5e3, #J/kg
    "V": 0.008, #m/s
    "T": 273.0 + 25,
    "z": 0.0
    }

ak_consts = {
    "Cp": 1005.46, # J/Kg.K
    "D": 4e-6, #1/s
    "beta": 0.42, # 0.42
    "delta_f": 35.0, # W/m^2?? - Check how this fits in with your equations.
    "Q_rad": -1.0, # K/day
    "V": 0.005, #m/s
    "gamma": 5.0, #K/km
    "theta_0_ft": 298.0, #K
    "theta_0": 301.0 #K
    } # note that these are related to Ann-Kristin's potential temperature values.

