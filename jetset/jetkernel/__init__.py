import os
import jetkernel

#package_dir=os.path.dirname(BlazarSED.__file__)

os.environ['BLAZARSED']=os.path.dirname(os.path.abspath(__file__))
print 'package_dir',os.path.dirname(os.path.abspath(__file__))