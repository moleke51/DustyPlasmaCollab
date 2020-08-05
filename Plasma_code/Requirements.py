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