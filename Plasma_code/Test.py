def eval_input(x):
    x = x.replace('^','**')
    if '**' in x:
        x = x.split('**')
        a = x[0].split('*')
        b = x[1].split('*')
        A = 1
        for i in range(len(a)-1):
            A *= float(a[i])
        for i in range(1,len(b)):
            A *= float(b[i])
        B = float(a[-1])**float(b[0])
        X = A*B
    else:
        x = x.split('*')
        X = 1
        for i in range(len(x)):
            X *= float(x[i])
    return X
def is_valid(user_input,name,requirements):
    check = True
    message = ''
    if 'is_num' in requirements:
        try:
            number = float(eval_input(user_input))
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
            return eval_input(user_input)
        else:
            print(message)
    else:
        return user_input

reqs = ['is_num','is_int','>=0']
N_max = 400
reqs.append('<='+str(N_max))
A = 0
#while A != '1':
    #A = is_valid(input('Input a number: '),'number',reqs)
    #if A != None:
        #print(A)

string = '1Y'

prefixes = {'Y': '*10**24', 'Z': '*10**21', 'E': '*10**18', 'P': '*10**15', 'T': '*10**12', 'G': '*10**9', 'M': '*10**6', 'k': '*10**3', 'm': '*10**-3', 'u': '*10**-6', 'n': '*10**-9', 'p': '*10**-12', 'f': '*10**-15', 'a': '*10**-18', 'z': '*10**-21', 'y': '*10**-24'}
print(string.translate(str.maketrans(prefixes)))