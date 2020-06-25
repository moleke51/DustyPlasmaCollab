import os
import Model 
FileList = os.listdir('Plasma_code/Models/')
modellist = []
for file in FileList:
    if ".py" in file:
        m = model.model(file,Theta,mu,z,alpha,upsilon)
        modellist.append(m)


__import__(self._name)