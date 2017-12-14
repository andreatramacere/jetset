import sys,os

# !in this way docs source module, not the installed one
sys.path.insert(0, os.path.abspath('../../'))


import BlazarSEDFit
import pkgutil

package = BlazarSEDFit

print package.__path__

def make_module_uml(modname):
    cmd="pyreverse  -AS -mn -k BlazarSEDFit.%s -o png -p BlazarSEDFit"%modname
    os.system(cmd)
    cmd="mv classes_BlazarSEDFit.png  modules_doc/classes_%s.png"%modname
    os.system(cmd)

def make_apidoc(mod_list):
    f=open('code.rst','w')
    text="""

Modules
=============

 In the following the package modules are listed.
 

=============

.. toctree::
    :maxdepth: 2
    """
    print>>f,text
    
    for modname in mod_list:
        print>>f,"    %s    <modules_doc/%s.rst>"%(modname,modname)
        
        f1=open('modules_doc/%s.rst'%modname,'w')
        
        text="""
.. automodule:: BlazarSEDFit.%s
   :members:
   :private-members:
   :undoc-members:
   :show-inheritance:
   
   """%modname
        
        print>>f1,text

        f1.close()
    
    f.close()
    
def main():
    mod_list=[]
    for importer, modname, ispkg in pkgutil.iter_modules(package.__path__):
        if ispkg==False:
            print "Found submodule %s " % (modname)
            print "generating classes uml graph"
            mod_list.append(modname)
            
    print mod_list
    
    make_apidoc(mod_list)
    
    for modname in mod_list:
            
            
            make_module_uml(modname)
            print


if __name__ == "__main__":
    main()
    