import os
import sys
sys.path.insert(1, 'Plasma_code/Models/')
class model:
    def __init__(self,filename,Theta,mu,z,alpha,upsilon):
        self._name = filename.strip('.py')
        self._Theta = float(Theta)
        self._mu = float(mu)
        self._z = float(z)
        self._alpha = float(alpha)
        self._upsilon = float(upsilon)
    def __repr__(self):
        return f'Model: {self._name}, at Theta = {self._Theta}, mu = {self._mu}, z = {self._z}, alpha = {self._alpha} and upsilon = {self._upsilon}'
    def get_name(self):
        return self._name
    def get_colour(self):
        return getattr(__import__(self._name),'colour')()
    def priority(self):
        return getattr(__import__(self._name), 'priority')(self._Theta,self._alpha,self._upsilon)
    def potential_finder(self):
        return getattr(__import__(self._name), 'potential_finder')(self._Theta,self._mu,self._z,self._alpha,self._upsilon)
    
def modelpicker(path,Theta,mu,z,alpha,upsilon):
    FileList = os.listdir(path)
    modellist = []
    for File in FileList:
        if ".py" in File:
            m = model(File,Theta,mu,z,alpha,upsilon)
            modellist.append(m)
    return modellist

