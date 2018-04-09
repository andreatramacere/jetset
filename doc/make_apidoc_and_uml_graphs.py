from __future__ import absolute_import, division, print_function,unicode_literals

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)


import sys,os



# !in this way docs source module, not the installed one
sys.path.insert(0, os.path.abspath('../../'))


import jetset
import pkgutil

package = jetset

print( package.__path__)

def make_module_uml(modname):
    cmd="pyreverse  -AS -mn -k jetset.%s -o png -p jetset"%modname
    os.system(cmd)
    cmd="mv classes_jetset.png  source/api/classes_%s.png"%modname
    os.system(cmd)


def  make_apidoc_automod(mod_list):
    f = open('modules.rst', 'w')
    text = """

Reference/API
=============
        """
    print(text,file=f)

    for modname in mod_list:


        text = """.. automodapi:: jetset.%s """ % (modname)
        print(text, file=f)
    f.close()
    cmd = "mv modules.rst  source/api/"
    os.system(cmd)

def make_apidoc(mod_list):
    f=open('modules.rst','w')
    text="""

Modules
=============

In the following the package modules are listed.
 
.. toctree::
   :maxdepth: 2
    """
    print(text, file=f)
    
    for modname in mod_list:
        print("   %s    <%s.rst>"%(modname,modname),file=f)
        
        f1=open('source/api/%s.rst'%modname,'w')

        under_title="-"*len(modname)

        text=""".. automodapi:: jetset.%s"""%(modname)

        print(text, file=f1)

        f1.close()
    
    f.close()

    cmd="mv modules.rst  source/api/"
    os.system(cmd)

def main():
    mod_list=[]
    for importer, modname, ispkg in pkgutil.iter_modules(package.__path__):
        if ispkg==False:
            print("Found submodule %s " % (modname))
            print("generating classes uml graph")
            mod_list.append(modname)
            
    print (mod_list)
    
    make_apidoc(mod_list)
    
    for modname in mod_list:
            
            
            #make_module_uml(modname)
            pass


if __name__ == "__main__":
    main()
    