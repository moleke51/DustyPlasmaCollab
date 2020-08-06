#============================IMPORT PACKAGES============================#
import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
import periodictable as pt
from scipy import integrate
from termcolor import colored
import os
import sys
sys.path.insert(1, 'DustyPlasmaCollab/Plasma_code/Models/')
#==============================CONSTANTS================================#
#We define some of the universal constants
e = 1.60e-19 #[C] The charge on an electron
epsilon_0 = 8.85e-12 #[F][m^-1] The permitivity of free space
k_B = 1.38e-23 #[m^2][kg][s^-2][K^-1] The Boltzmann constant
m_e = 9.11e-31 #[kg] The mass of an electron
u = 1.66e-27 #[kg] The mass of a nucleon
numgraphpoints = 10000
colour = 'red'
#This allows the chemical symbols or element names to be entered in either upper or lower case.
def speciesinput():
    word = input("Enter the plasma ion species; ")
    if len(word)==1 : 
        word = word.upper()
    elif len(word)==2 :
        word = word.lower().capitalize()
    else:
        word = word.lower()
    return(word)
#Define the periodic table
elements = pt.core.default_table()
elementList = pt.core.define_elements(elements,globals())


def choosespecies():
    species = speciesinput()
    while (species in elementList) == False:
        print('This species does not exist')
        species = speciesinput()
    index = np.floor(elementList.index(species)/2)
    if index == 119: #Deuterium
        index = 1
        m_a = pt.deuterium.mass
    elif index == 120: #Tritium
        index = 1
        m_a = pt.tritium.mass
    else:   
        m_a = (elements[index].mass) #[kg]
    m_i = (m_a)*u #[kg]
    mu = np.sqrt(m_i/m_e) #Mu value
    z_max = elements[index].number
    return m_i,mu,z_max


class IsInteger:
    """Class for the condition the input is an integer"""
    def check(self,input_number):
        try:
            if input_number == int(input_number):
                return True
            return False
        except TypeError:
            return False
    def error_message(self):
        return 'an integer'

class GreaterThan:
    """Class for the condition greater than x"""
    def __init__(self,x):
        self._x = x
    
    def check(self,input_number):
        try:
            if self._x == None:
                return True
            elif input_number > self._x:
                return True
            return False
        except TypeError:
            return False
    def error_message(self):
        return f'greater than {self._x}'

class GreaterThanEqualTo:
    """Class for the condition greater than or equal to x"""
    def __init__(self,x):
        self._x = x
    
    def check(self,input_number):
        try:
            if self._x == None:
                return True
            elif input_number >= self._x:
                return True
            return False
        except TypeError:
            return False
    def error_message(self):
        return f'greater than or equal to {self._x}'

class LessThan:
    """Class for the condition less than x"""
    def __init__(self,x):
        self._x = x
    
    def check(self,input_number):
        try:
            if self._x == None:
                return True
            elif input_number < self._x:
                return True
            return False
        except TypeError:
            return False
    def error_message(self):
        return f'less than {self._x}'

class LessThanEqualTo:
    """Class for the condition less than or equal to x"""
    def __init__(self,x):
        self._x = x
    
    def check(self,input_number):
        try:
            if self._x == None:
                return True
            elif input_number <= self._x:
                return True
            return False
        except TypeError:
            return False
    def error_message(self):
        return f'less than or equal to {self._x}'

class Norm:
    def __init__(self,dictionary_list = None):
        self._dictlist = dictionary_list
    def update(self,dictionary_list):
        self._dictlist = dictionary_list
    def getvarvalue(self,varname):
        for _vardict in self._dictlist:
            if _vardict.get('var_name') == varname:
                return _vardict.get('Value')


class Normunity(Norm):
    def getnormfactor(self):
        return 1
class Norm_T_i(Norm):

    def getnormfactor(self):
        _T_e = self.getvarvalue('Electron temperature')
        return _T_e

class Norm_a(Norm):
    def getnormfactor(self):
        _T_e = self.getvarvalue('Electron temperature')
        _n_e = self.getvarvalue('Electron density')
        return np.sqrt(epsilon_0*k_B*_T_e/(_n_e*(e**2)))

class Norm_v(Norm):
    def getnormfactor(self):
        _T_i = self.getvarvalue('Ion temperature')
        _m_i = self.getvarvalue('Ion mass')
        return np.sqrt(2*k_B*_T_i/_m_i)

T_e_dict = {
    'var_name' : 'Electron temperature',
    'Requirements' : [GreaterThan(0)],
    'Unit' : 'Kelvin',
    'Unit_symbol' : 'K',
    'isvariable' : True
}
nfTi = Norm_T_i()
T_i_dict = {
    'var_name' : 'Ion temperature',
    'Requirements' : [GreaterThan(0)],
    'Norm_factor' : nfTi,
    'Norm_var_name' : 'Theta',
    'Unit' : 'Kelvin',
    'Unit_symbol' : 'K',
    'isvariable' : True
}
nfz = Normunity()
z_dict = {
    'var_name' : 'Relative ion charge',
    'Requirements' : [GreaterThan(0),IsInteger()],
    'Norm_var_name' : 'z',
    'Norm_factor' : nfz,
    'Unit' : None,
    'Unit_symbol' : 'e',
    'isvariable' : False
}

m_i_dict = {
    'var_name' : 'Ion mass',
    'Requirements' : [GreaterThan(0)],
    'Norm_var_name' : 'mu',
    'Unit_symbol' : 'kg',
    'isvariable' : False
}

n_e_dict = {
    'var_name' : 'Electron density',
    'Requirements' : [GreaterThan(0)],
    'Unit' : 'electrons per meter cubed',
    'Unit_symbol' : 'm^-3',
    'isvariable' : False
}
nfa = Norm_a()
a_dict = {
    'var_name' : 'Dust radius',
    'Requirements' : [GreaterThanEqualTo(0)],
    'Norm_factor' : nfa,
    'Norm_var_name' : 'alpha',
    'Unit' : 'meters',
    'Unit_symbol' : 'm',
    'isvariable' : True
}
nfv = Norm_v()
v_dict = {
    'var_name' : 'Flow speed',
    'Requirements' : [GreaterThanEqualTo(0)],
    'Norm_factor' : nfv,
    'Norm_var_name' : 'upsilon',
    'Unit' : 'meters per second',
    'Unit_symbol' : 'ms^-1',
    'isvariable' : True,
    'default value' : 0
}



dict_list = [T_i_dict,T_e_dict,z_dict,m_i_dict,n_e_dict,a_dict,v_dict]


def eval_input(x):
    x = x.replace('^','**')
    if '**' in x:
        x = x.split('**',)
        a = x[0].split('*')
        b = x[-1].split('*')
        A = 1
        for i in range(len(a)-1):
            A *= float(a[i])
        for i in range(1,len(b)):
            A *= float(b[i])
        if len(x) == 3:
            c,d = x[1].split('*')
            B = (float(a[-1])**float(c))*(float(d)**float(b[0])) 
        else:
            B = float(a[-1])**float(b[0])
        X = A*B
    else:
        x = x.split('*')
        X = 1
        for i in range(len(x)):
            X *= float(x[i])
    return X

def is_valid(name,requirements,var_counter,units = None,isvariable = True):
    if requirements == None:
        requirements = []
    units_message = '; '
    if units != None:
        units_message = f' with units of {units}; '
    user_input = input(f'Enter the {name.lower()} value{units_message.lower()}')
    if user_input.lower() == 'variable':
        if isvariable == True:
            if user_input.lower() == 'variable' and var_counter == False:
                maximum = None
                minimum = None
                while minimum == None:
                    minimum, placeholder = is_valid('minimum '+name,requirements,True,units)
                    requirements.append(GreaterThan(minimum))
                while maximum == None:
                    maximum, placeholder = is_valid('maximum '+name,requirements,True,units)
                return(np.linspace(minimum,maximum,numgraphpoints),True)
            elif var_counter == True: 
                print('The variable has already been selected')
                return is_valid(name,requirements,var_counter,units = None)
        else:
            print('This parameter cannot be varied')
            return is_valid(name,requirements,var_counter,units = None,isvariable = False)
    try:
        prefixes = {'Y': '*10**24', 'Z': '*10**21', 'E': '*10**18', 'P': '*10**15', 'T': '*10**12', 'G': '*10**9', 'M': '*10**6', 'k': '*10**3', 'm': '*10**-3', 'u': '*10**-6', 'n': '*10**-9', 'p': '*10**-12', 'f': '*10**-15', 'a': '*10**-18', 'z': '*10**-21', 'y': '*10**-24'}
        user_input = eval_input((user_input.replace('eV','ev').replace('ev','*11594')).translate(str.maketrans(prefixes)))
        errors = []
        for req in requirements:
            if req.check(user_input) == False:
                errors.append(req.error_message())
        if len(errors) == 0:
            return user_input, var_counter
        error_message = f'The {name.lower()} value should be '
        error_message += errors[0]
        if len(errors) > 1:
            for i in range(1,len(errors)-1):
                error_message += ', '+errors[i]
            error_message += ' and '+errors[-1]
        print(error_message)
        return is_valid(name,requirements,var_counter,units,isvariable)
        
    except ValueError:
        print(f'{name.lower()} must be a number.')
        return is_valid(name,requirements,var_counter,units,isvariable)
    except NameError:
        print(f'{name.lower()} must be a number.') 
        return is_valid(name,requirements,var_counter,units,isvariable)
    except TypeError:
        print(f'{name.lower()} must be a number.')
        return is_valid(name,requirements,var_counter,units,isvariable)


class Model:
    def __init__(self,filename,dictionarylist,dimensionless):
        self._name = filename[:-3]
        self._dictlist = dictionarylist
        self._dim = dimensionless
    def __repr__(self):
        self._rep = f'Model: {self._name}'
        if self._dim == False:
            for _vardict in self._dictlist:
                _name = _vardict.get('var_name')
                _val = _vardict.get('Value')
                _us = _vardict.get('Unit_symbol')
                self._rep += f', {_name} = {_val}{_us}'
        else:
            for _vardict in self._dictlist:
                _name = _vardict.get('Norm_var_name')
                _val = _vardict.get('Norm_value')
                if _name != None:
                    self._rep += f', {_name} = {_val}'
        return self._rep
    def get_name(self):
        return self._name
    def get_colour(self):
        return getattr(__import__(self._name),'colour')()
    def get_info(self):
        return getattr(__import__(self._name),'get_info')()
    def priority(self):
        return getattr(__import__(self._name), 'priority')(self._dictlist)
    def potential_finder(self):
        return getattr(__import__(self._name), 'potential_finder')(self._dictlist)
    
#Define the Debye-Huckle potential.
def DH_potential_to_charge(a,Phi_a,lambda_d):
    x = 4*np.pi*(epsilon_0)*a*Phi_a*np.exp((a)/(lambda_d))
    return x

#It should be noted that the Debye-Huckle potential reduces into the point charge potential when a<<lambda_d
def Spherical_potential_to_charge(a,Phi_a):
    x = 4*np.pi*(epsilon_0)*a*Phi_a
    return x
class PotentialCalculator:

    '''
    Class for tracking the variables involved in the calculation
    of the potential at the surface of a dust grain
    '''

    def __init__(self,dictionary_list):
        self._dictlist = dictionary_list
        self._variablecounter = False
        self._varnamelist = []
        for _vardict in self._dictlist:
            self._varnamelist.append(_vardict.get('var_name'))


    def check_variables(self):
        return self._varnamelist


    def initialise(self):
        _choice = None 
        while _choice != 'y' and _choice != 'n':
            _choice = input('Do you want to use dimensionless variables (y/n); ').lower()
        if _choice == 'y':
            self._dimensionless = True
            for _vardict in self._dictlist:
                if _vardict.get('Norm_var_name') != None:
                    _vardict['Norm_value'], self._variablecounter = is_valid(_vardict.get('Norm_var_name'),_vardict.get('Requirements'),self._variablecounter,None,_vardict.get('isvariable'))
            
            if self._variablecounter == True:
                for _vardict in self._dictlist:
                    if _vardict.get('Norm_var_name') != None:
                        if type(_vardict.get('Norm_value')) == np.ndarray:
                            self._variabletracker = _vardict.get('Norm_var_name')
                        else:
                            _vardict.update({'Norm_value' : np.ones(numgraphpoints)*_vardict.get('Norm_value')}) 
        else:
            self._dimensionless = False
            m_i,mu, z_max = choosespecies()
             
            for _vardict in self._dictlist:
                if _vardict.get('var_name') == 'Relative ion charge':
                    _vardict.update({'Requirements' : [GreaterThan(0),IsInteger(),LessThanEqualTo(z_max)]})
                if _vardict.get('var_name') == 'Ion mass':
                    _vardict['Value'] = m_i
                    _vardict['Norm_value'] = mu
                    
                else:
                    _vardict['Value'], self._variablecounter = is_valid(_vardict.get('var_name'),_vardict.get('Requirements'),self._variablecounter,_vardict.get('Unit'),_vardict.get('isvariable'))
            if self._variablecounter == True:
                for _vardict in self._dictlist:
                    if type(_vardict.get('Value')) == np.ndarray:
                        self._variabletracker = _vardict.get('var_name')
                    elif _vardict.get('var_name') == 'Ion mass':
                        _vardict.update({'Norm_value' : np.ones(numgraphpoints)*_vardict.get('Norm_value')})
                    else:
                        _vardict.update({'Value' : np.ones(numgraphpoints)*_vardict.get('Value')})  
            for _vardict in self._dictlist:
                if _vardict.get('Norm_factor') != None:
                    _vardict.get('Norm_factor').update(self._dictlist)
                    _vardict.update({'Norm_factor' : _vardict.get('Norm_factor').getnormfactor()})
                    _vardict['Norm_value'] = _vardict.get('Value')/_vardict.get('Norm_factor')
                #print(_vardict.get('Norm_value'))
        
                               

    def get_Norm_Potential(self):
        self.initialise()
        FileList = os.listdir('DustyPlasmaCollab/Plasma_code/Models/')
        if self._variablecounter == False:
            modellist = []
            for File in FileList:
                if ".py" in File:
                    m = Model(File,self._dictlist,self._dimensionless)
                    modellist.append(m)
            
            priority = 0
            for model in modellist:
                __import__(model.get_name())
                if model.priority() > priority:
                    
                    priority = model.priority()
                    modelindex = modellist.index(model)
        
        
            print(colored(modellist[modelindex].__repr__(),colour))
            _Phi = modellist[modelindex].potential_finder()
            if self._dimensionless == True:
                print(f'The normalised potential is: {_Phi}')
            else:
                for _vardict in self._dictlist:
                    if _vardict.get('var_name') == 'Dust radius':
                        _a =  _vardict.get('Value')
                        _alpha = _vardict.get('Norm_value')
                        _lambda_d = _a/_alpha
                    if _vardict.get('var_name') == 'Electron temperature':
                        _T_e =  _vardict.get('Value')
                        _k_B = k_B
                        _e = e
                        _phi = (_Phi * k_B * _T_e)/(e) 
                _Q = DH_potential_to_charge(_a,_phi,_lambda_d)
            
                print(f'The dust grain surface potential is {_phi} V')
                print(f'The charge on the dust grain is {_Q} C')
                print(f'The relative charge of the dust grain is {_Q/_e}')

            _info = modellist[modelindex].get_info()
            print(colored(info,colour))
        else:
            if self._dimensionless == True:
                val = 'Norm_value'
            else:
                val = 'Value'
            _Phiarray = np.ones(len(self._dictlist[0].get(val)))
            _modelname = []
            for i in range(len(self._dictlist[0].get(val))):
                _newdictlist = []
                for _vardict in self._dictlist:
                    _newdict = _vardict.copy()
                    if type(_vardict.get('Norm_value')) != type(None):
                        _newdict.update({'Norm_value' : _vardict.get('Norm_value')[i]})
                        _newdictlist.append(_newdict)
                
                modellist = []
                for File in FileList:
                    if ".py" in File:
                        m = Model(File,_newdictlist,self._dimensionless)
                        modellist.append(m)
                
                priority = 0
                for model in modellist:
                    __import__(model.get_name())
                    if model.priority() > priority:
                        
                        priority = model.priority()
                        modelindex = modellist.index(model)
            
                #print(modellist)
                #print(colored(modellist[modelindex].__repr__(),colour))
                _Phiarray[i] =  modellist[modelindex].potential_finder()#,modellist[modelindex].get_info()
                _modelname.append(modellist[modelindex].get_name())
            #print(_Phiarray)

            
            plt.xlabel(self._variabletracker)
            if self._dimensionless == True:
                for _vardict in self._dictlist:
                    if _vardict.get('Norm_var_name') == self._variabletracker:
                        _X = _vardict.get('Norm_value')
                plt.plot(_X,_Phiarray)
                plt.ylabel('Normalised potential')
                
            else:
                for _vardict in self._dictlist:
                    if _vardict.get('var_name') == self._variabletracker:
                        _X = _vardict.get('Value')
                    if _vardict.get('var_name') == 'Dust radius':
                        _a =  _vardict.get('Value')
                        _alpha = _vardict.get('Norm_value')
                        _lambda_d = _a/_alpha
                    if _vardict.get('var_name') == 'Electron temperature':
                        _T_e =  _vardict.get('Value')
                        _k_B = k_B*np.ones(len(_T_e))
                        _e = e*np.ones(len(_T_e))
                        _phiarray = (_Phiarray * k_B * _T_e)/(e) 
                _Q = DH_potential_to_charge(_a,_phiarray,_lambda_d)
                plt.plot(_X,_phiarray)
                plt.ylabel('Potential [V]')
            plt.grid()
            plt.show()
                
            #return _Phiarray,_modelname           
        
       
pc = PotentialCalculator(dict_list)
pc.get_Norm_Potential()
#print(colored(f'The normalised potential is: {Phi}',colour))
#print(colored(info,colour))


