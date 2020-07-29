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
    prefixes = {'Y': '*10**24', 'Z': '*10**21', 'E': '*10**18', 'P': '*10**15', 'T': '*10**12', 'G': '*10**9', 'M': '*10**6', 'k': '*10**3', 'm': '*10**-3', 'u': '*10**-6', 'n': '*10**-9', 'p': '*10**-12', 'f': '*10**-15', 'a': '*10**-18', 'z': '*10**-21', 'y': '*10**-24'}
    check = True
    message = ''
    if user_input.lower() == 'variable':
        if user_input.lower() == 'variable' and var_counter == 0:
            maximum = None
            minimum = None
            while maximum == None:
                maximum, a = is_valid('maximum '+name,requirements,1,units)
            requirements.append(f'<{maximum}')
            while minimum == None:
                minimum, a = is_valid('minimum '+name,requirements,1,units)
            return(np.linspace(minimum,maximum,100000).tolist(),var_counter+1)
        elif var_counter == 1: 
            print('The variable has already been selected')
        else:
            print('This parameter cannot be varied')
    if 'is_num' in requirements:
        try:
            str_number = eval_input((user_input.replace('eV','ev').replace('ev','*11594')).translate(str.maketrans(prefixes)))
            number = float(str_number)
            if 'is_int' in requirements:
                if '.' in user_input:
                    check = False
                    if message == '':
                        message = f'{name} must be an integer'
                    else:
                        message += f', {name} must be an integer'
            if '>0' in requirements:
                if number <= 0:
                    check = False
                    if message == '':
                        message = f'{name} must be greater than zero'
                    else:
                        message += f', {name} must be greater than zero'
            if '>=0' in requirements: 
                if number < 0:
                    check = False
                    if message == '':
                        message = f'{name} must be greater than or equal to zero'
                    else:
                        message += f', {name} must be greater than or equal to zero' 
            for req in requirements:
                if '<=' in req:
                    max_num = float(req.strip('<='))
                    if number > max_num:
                        check = False
                        if message == '':
                            message = f'{name} must be smaller than {max_num}'
                        else:
                            message += f', {name} must be smaller than {max_num}'
            for req in requirements:
                if '<' in req and '<=' not in req:
                    max_num = float(req.strip('<'))
                    if number >= max_num:
                        check = False
                        if message == '':
                            message = f'{name} must be smaller than {max_num}'
                        else:
                            message += f', {name} must be smaller than {max_num}'
        except ValueError:
            check = False
            message = (f'{name} must be a number.')
            
        except NameError:
            check = False
            message = (f'{name} must be a number.')
            
        except TypeError:
            check = False
            message = (f'{name} must be a number.')
            
        if check == True:
            return(str_number, var_counter)
        else:
            print(message)
            return(None,var_counter)
    else:
        return(user_input,var_counter)