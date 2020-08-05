#============================IMPORT PACKAGES============================#
import numpy as np
import scipy as sp
import scipy.special as sps
import matplotlib.pyplot as plt
import periodictable as pt
from scipy import integrate
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
    return mu,z_max


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
        return is_valid(name,requirements,var_counter,units = None)
        
    except ValueError:
        print(f'{name.lower()} must be a number.')
        return is_valid(name,requirements,var_counter,units = None)
    except NameError:
        print(f'{name.lower()} must be a number.') 
        return is_valid(name,requirements,var_counter,units = None)
    except TypeError:
        print(f'{name.lower()} must be a number.')
        return is_valid(name,requirements,var_counter,units = None)

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



T_e_dict = {
    'var_name' : 'Electron temperature',
    'Requirements' : [GreaterThan(0)],
    'Unit' : 'Kelvin',
    'isvariable' : True
}
nfTi = Norm_T_i()
T_i_dict = {
    'var_name' : 'Ion temperature',
    'Requirements' : [GreaterThan(0)],
    'Norm_factor' : nfTi,
    'Norm_var_name' : 'Theta',
    'Unit' : 'Kelvin',
    'isvariable' : True
}

z_dict = {
    'var_name' : 'Relative ion charge',
    'Requirements' : [GreaterThan(0),IsInteger()],
    'Norm_var_name' : 'z',
    'Unit' : None,
    'isvariable' : False
}
nfmi = Normunity()
m_i_dict = {
    'var_name' : 'Ion mass',
    'Norm_factor' : nfmi,
    'Norm_var_name' : 'mu',
    'isvariable' : False
}


dict_list = [T_e_dict,T_i_dict,z_dict,m_i_dict]


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
        else:
            self._dimensionless = False
            mu, z_max = choosespecies()
             
            for _vardict in self._dictlist:
                if _vardict.get('var_name') == 'Relative ion charge':
                    _vardict.update({'Requirements' : [GreaterThan(0),IsInteger(),LessThanEqualTo(z_max)]})
                if _vardict.get('var_name') == 'Ion mass':
                    _vardict['Value'] = mu
                    
                else:
                    _var, self._variablecounter = is_valid(_vardict.get('var_name'),_vardict.get('Requirements'),self._variablecounter,_vardict.get('Unit'),_vardict.get('isvariable'))
                    _vardict['Value'] = _var
            if self._variablecounter == True:
                for _vardict in self._dictlist:
                    if type(_vardict.get('Value')) == np.ndarray:
                        self._variabletracker = _vardict.get('var_name')
                    else:
                        _vardict.update({'Value' : np.ones(numgraphpoints)*_vardict.get('Value')})  
            for _vardict in self._dictlist:
                if _vardict.get('Norm_factor') != None:
                    _vardict.get('Norm_factor').update(self._dictlist)
                    _vardict.update({'Norm_factor' : _vardict.get('Norm_factor').getnormfactor()})
                    _vardict['Norm_value'] = _vardict.get('Value')/_vardict.get('Norm_factor')
                print(_vardict.get('Norm_value'))
                               
        
       
pc = PotentialCalculator(dict_list)
pc.initialise()


