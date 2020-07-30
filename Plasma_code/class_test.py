import numpy as np

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
            if input_number > self._x:
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
            if input_number >= self._x:
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
            if input_number < self._x:
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
            if input_number <= self._x:
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

def is_valid(name,requirements,var_counter,units = None):
    units_message = '; '
    if units != None:
        units_message = f' with units of {units}; '
    user_input = input(f'Enter the {name} value{units_message}')
    if user_input.lower() == 'variable':
        if user_input.lower() == 'variable' and var_counter == 0:
            maximum = None
            minimum = None
            while minimum == None:
                minimum, placeholder = is_valid('minimum '+name,requirements,1,units)
                requirements.append(GreaterThan(minimum))
            while maximum == None:
                maximum, placeholder = is_valid('maximum '+name,requirements,1,units)
            return(np.linspace(minimum,maximum,100000).tolist(),var_counter+1)
        elif var_counter == 1: 
            print('The variable has already been selected')
            return is_valid(name,requirements,var_counter,units = None)
        else:
            print('This parameter cannot be varied')
            return is_valid(name,requirements,var_counter,units = None)
    try:
        prefixes = {'Y': '*10**24', 'Z': '*10**21', 'E': '*10**18', 'P': '*10**15', 'T': '*10**12', 'G': '*10**9', 'M': '*10**6', 'k': '*10**3', 'm': '*10**-3', 'u': '*10**-6', 'n': '*10**-9', 'p': '*10**-12', 'f': '*10**-15', 'a': '*10**-18', 'z': '*10**-21', 'y': '*10**-24'}
        user_input = eval_input((user_input.replace('eV','ev').replace('ev','*11594')).translate(str.maketrans(prefixes)))
        errors = []
        for req in requirements:
            if req.check(user_input) == False:
                errors.append(req.error_message())
        if len(errors) == 0:
            return user_input, var_counter
        error_message = f'The {name} value should be '
        error_message += errors[0]
        if len(errors) > 1:
            for i in range(1,len(errors)-1):
                error_message += ', '+errors[i]
            error_message += ' and '+errors[-1]
        print(error_message)
        return is_valid(name,requirements,var_counter,units = None)
        
    except ValueError:
        print(f'{name} must be a number.')
        return is_valid(name,requirements,var_counter,units = None)
    except NameError:
        print(f'{name} must be a number.') 
        return is_valid(name,requirements,var_counter,units = None)
    except TypeError:
        print(f'{name} must be a number.')
        return is_valid(name,requirements,var_counter,units = None)


print(is_valid('Test',[],1,'Smeckles'))
