from inspect import getmembers, isfunction
import Testmodule

class model:
    def __init__(self,filename,Theta,mu,z,alpha,upsilon):
        self._name = filename.strip('.py')
        self._Theta = float(Theta)
        self._mu = float(mu)
        self._z = float(z)
        self._alpha = float(alpha)
        self._upsilon = float(upsilon)
    def get_name(self):
        return self._name
    def priority(self):
        return getattr(__import__(self._name), 'priority')(self._Theta,self._mu,self._z,self._alpha,self._upsilon)
    def potential_finder(self,filename):
        return getattr(__import__(self._name), 'potential_finder')(self._Theta,self._mu,self._z,self._alpha,self._upsilon)
    


