from pathlib import Path
import numpy as np
import os
import Model as model

Theta = 1
mu = 43
z = 1
alpha = 1
upsilon = 0

FileList = os.listdir('Plasma_code/Models/')
modellist = []
for file in FileList:
    if ".py" in file:
        m = model.model(file,Theta,mu,z,alpha,upsilon)
        modellist.append(m)

priority = 0
for model in modellist:
    if model.priority() > priority:
        print(model)
        priority = model.priority()
        modelindex = modellist.index(model)

print(modellist[modelindex], priority)



