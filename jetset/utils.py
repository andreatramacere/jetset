__author__ = "Andrea Tramacere"

import re
import inspect
import sys
import warnings
import os
import json



__all__=['check_frame','unexpected_behaviour']



def check_frame(frame):
    allowed=['obs','src','blob']
    if frame not in allowed:
        raise RuntimeError('rest frame', frame, 'not allowed',allowed)

def unexpected_behaviour():
    raise RuntimeError('the code reached a condition that should never happen!')


def clean_var_name(s):
    _s = s
    s.replace('-', '_')
    s.replace(' ', '_')

    # Remove invalid characters
    s = re.sub('[^0-9a-zA-Z_]', '', s)

    # Remove leading characters until we find a letter or underscore
    s = re.sub('^[^a-zA-Z_]+', '', s)

    if len(s) == 0:
        raise RuntimeError('the string ', _s, 'is not valid, please change it (do not use leading numbers or math symbols)')
    return s


class NoTraceBackWithLineNumber(Exception):
    def __init__(self, msg):
        try:
            ln = sys.exc_info()[-1].tb_lineno
        except AttributeError:
            ln = inspect.currentframe().f_back.f_lineno
        self.args = "{0.__name__} (line {1}): {2}".format(type(self), ln, msg),
        sys.exit(self)

def new_version_warning():
    m = '\n\n' + '*'*80 + '\n'
    m+= 'Something wrong has happened. Please, look at the exception message.\n'
    m+= '*' * 80 + '\n'
    warnings.warn(m)


def parameters_warning():
    pass


def old_model_warning():
    m = '\n\n' + '*'*80 + '\n'
    m+= 'you are loading a model supported for version<1.1.0, starting from version 1.1.0 \n'
    m+= 'the saved model has changed,  plase update to the new model the new format, \n'
    m += 'by saving it with this version\n'
    m+= '*' * 80 + '\n'
    warnings.warn(m)

class JetkerneltException(Exception):

    def __init__(self, message='Jeset  exception', debug_message=''):
        super(JetkerneltException, self).__init__(message)
        self.message=message
        self.debug_message=debug_message


#def str_hook(pairs):
#    new_pairs = []
#    for key, value in pairs:
#        if isinstance(value, unicode):
#            value = value.encode('utf-8')
#        if isinstance(key, unicode):
#            key = key.encode('utf-8')
#        new_pairs.append((key, value))
#    return dict(new_pairs)



def safe_run(func):

    def func_wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
           message =  'the jetkernel failed\n'
           message += '\n exception message: \n'
           message += '%s'%str(e)
           new_version_warning()

           raise JetkerneltException(message=message)


    return func_wrapper


def set_str_attr(obj,name,val):
    #print('set obj', obj,'name',name ,'to', val)
    try:

        try:
            setattr(obj, name,val)
        except:
            setattr(obj, name, val.encode('ascii'))
    except Exception as e:
        raise RuntimeError('error setting attr',name,'execption:',e)



def get_info():
    with open(os.path.dirname(__file__) + '/pkg_info.json') as fp:
        _info = json.load(fp)

    return _info
