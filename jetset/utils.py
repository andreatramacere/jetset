__author__ = "Andrea Tramacere"

import re
import inspect
import sys
import warnings
import os
import json



__all__=['check_frame','unexpected_behaviour','get_nested_attr','set_nested_attr']



def check_frame(frame):
    """Check frame.
    
    Parameters
    ----------
    frame : object
        Reference frame for data/model values.
    """
    allowed=['obs','src','blob']
    if frame not in allowed:
        raise RuntimeError('rest frame', frame, 'not allowed',allowed)

def unexpected_behaviour():
    """Unexpected behaviour."""
    raise RuntimeError('the code reached a condition that should never happen!')


def get_nested_attr(obj, name):
    #print("==> get",obj,name)
    """Return nested attr.
    
    Parameters
    ----------
    obj : object
        Parameter controlling obj.
    name : object
        Name identifier.
    
    Returns
    -------
    object
        Requested value.
    """
    parts = name.split('.')
    for part in parts:
        obj = getattr(obj, part)
    return obj


def set_nested_attr(obj, name, val):
    """Set nested attr.
    
    Parameters
    ----------
    obj : object
        Parameter controlling obj.
    name : object
        Name identifier.
    val : object
        Value to assign.
    """
    parts = name.split('.')
    target = obj
    for part in parts[:-1]:
        target = getattr(target, part)
    setattr(target, parts[-1], val)


def clean_var_name(s):
    """Clean var name.
    
    Parameters
    ----------
    s : object
        Parameter controlling s.
    
    Returns
    -------
    object
        Computed result.
    """
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
    """Exception helper that reports message with originating line number."""
    def __init__(self, msg):
        """Create a new `NoTraceBackWithLineNumber` instance.
        
        Parameters
        ----------
        msg : object
            Parameter controlling msg.
        """
        try:
            ln = sys.exc_info()[-1].tb_lineno
        except AttributeError:
            ln = inspect.currentframe().f_back.f_lineno
        self.args = "{0.__name__} (line {1}): {2}".format(type(self), ln, msg),
        sys.exit(self)

def new_version_warning():
    """New version warning."""
    m = '\n\n' + '*'*80 + '\n'
    m+= 'Something wrong has happened. Please, look at the exception message.\n'
    m+= '*' * 80 + '\n'
    warnings.warn(m)


def parameters_warning():
    """Parameters warning."""
    pass


def old_model_warning():
    """Old model warning."""
    m = '\n\n' + '*'*80 + '\n'
    m+= 'you are loading a model supported for version<1.1.0, starting from version 1.1.0 \n'
    m+= 'the saved model has changed,  please update to the new model the new format, \n'
    m += 'by saving it with this version\n'
    m+= '*' * 80 + '\n'
    warnings.warn(m)

class JetkerneltException(Exception):
    """Custom exception for jetkernel runtime failures and wrappers."""

    def __init__(self, message='Jeset  exception', debug_message=''):
        """Create a new `JetkerneltException` instance.
        
        Parameters
        ----------
        message : str, optional
            Parameter controlling message.
        debug_message : str, optional
            Parameter controlling debug message.
        """
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

    """Safe run.
    
    Parameters
    ----------
    func : object
        Parameter controlling func.
    
    Returns
    -------
    object
        Computed result.
    """
    def func_wrapper(*args, **kwargs):
        """Func wrapper.
        
        Parameters
        ----------
        *args : tuple
            Additional positional arguments.
        **kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        object
            Computed result.
        """
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
    """Set str attr.
    
    Parameters
    ----------
    obj : object
        Parameter controlling obj.
    name : object
        Name identifier.
    val : object
        Value to assign.
    """
    try:
        if '.' in name:
            try:
                set_nested_attr(obj, name, val)
            except Exception:
                set_nested_attr(obj, name, val.encode('ascii'))
        else:
            try:
                setattr(obj, name,val)
            except Exception:
                setattr(obj, name, val.encode('ascii'))
    except Exception as e:
        raise RuntimeError('error setting attr',name,'execption:',e)


def set_particle_attr(obj, emitters_type):
    """Set ``core.PARTICLE`` using enum codes while keeping string fallback.

    Parameters
    ----------
    obj : object
        Backend blob-like object.
    emitters_type : str
        Particle selector (``'electrons'`` or ``'protons'``).
    """
    try:
        from .jetkernel import jetkernel as BlazarSED
        particle_map = {
            'electrons': BlazarSED.PARTICLE_ELECTRONS,
            'protons': BlazarSED.PARTICLE_PROTONS,
            'secondaries_el': BlazarSED.PARTICLE_SECONDARIES_EL,
            'primaries_el': BlazarSED.PARTICLE_PRIMARIES_EL,
        }
    except Exception:
        particle_map = {
            'electrons': 0,
            'protons': 1,
            'secondaries_el': 2,
            'primaries_el': 3,
        }

    if emitters_type not in particle_map:
        raise RuntimeError('emitters type', emitters_type, 'not valid', tuple(particle_map.keys()))

    try:
        set_nested_attr(obj, 'core.PARTICLE', particle_map[emitters_type])
    except Exception:
        # Compatibility with pre-enum builds where PARTICLE is still char[].
        set_str_attr(obj, 'core.PARTICLE', emitters_type)



def get_info():
    """Return info.
    
    Returns
    -------
    object
        Requested value.
    """
    with open(os.path.dirname(__file__) + '/pkg_info.json') as fp:
        _info = json.load(fp)

    return _info
