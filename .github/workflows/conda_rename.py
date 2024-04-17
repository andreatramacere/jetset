import os
import platform
import shutil
name=os.getenv('CONDABUILDJETSET')

if platform.system()=='Darwin':
    system_str = 'macos'
else:
    system_str = 'linux'

if 'arm64' in platform.platform():
    system_str= '%s_arm64'%system_str

new_name=name.replace('.tar.bz2','_conda_%s.tar.bz2'%system_str)

shutil.move(name,new_name)

os.makedirs('conda-binary')

shutil.move(new_name,'conda-binary')
