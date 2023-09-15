import os
import platform
import shutil
name=os.getenv('CONDABUILDJETSET')

if platform.system()=='Darwin':
    system_str = 'macos'
else:
    system_str = 'linux'

new_name=name.replace('.tar.bz2','_%s.tar.bz2'%system_str)

shutil.move(name,new_name)

os.makedirs('conda-binary')

shutil.move(new_name,'conda-binary')
