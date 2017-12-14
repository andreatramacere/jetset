#!/usr/bin/env python
import os,sys

try:
    import IPython
except ImportError:
    print 'IPython not installed, requires IPython'
    sys.exit(3)

import BlazarSEDFit 
init_script=os.path.dirname(BlazarSEDFit.__file__)
init_script=init_script+'/'+'interactive_init.py'

from optparse import OptionParser

op = OptionParser()
op.usage = '%prog [options]'
op.add_option("--qt", action="store_true", dest="qt")


ops,args = op.parse_args()

print ops
print args


command="""-c 'run %s'"""%init_script

 
    
if ops.qt is None:
    sys.argv=['ipython ','-i' , command] 
    print sys.argv
else:
    sys.argv=['ipython ', 'qtconsole','--pylab=qt',command] 
    print sys.argv
    



from IPython.frontend.terminal.ipapp import launch_new_instance
launch_new_instance()
        

