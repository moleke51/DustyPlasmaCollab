N = 7
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

#class Get_T_e:
    

T_e_dict = {
    'var_name' : 'Electron temperature',
    'Requirements' : [GreaterThan(0)],
    'Norm_factor' : 1
}

T_i_dict = {
    'var_name' : 'Ion temperature',
    'Requirements' : [GreaterThan(0)],
    'Norm_factor' : Get_T_e
}

