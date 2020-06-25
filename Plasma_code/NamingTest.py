from pathlib import Path
import numpy as np
import os
import Model as mdl

Theta = 1
mu = 43
z = 1
alpha = 1
upsilon = 0

modellist = mdl.modelpicker('Plasma_code/Models/',Theta,mu,z,alpha,upsilon)
priority = 0
for model in modellist:
    __import__(model.get_name())
    if model.priority() > priority:
        
        priority = model.priority()
        modelindex = modellist.index(model)

print(modellist[modelindex], modellist[modelindex].potential_finder(modellist[modelindex].get_name()))



