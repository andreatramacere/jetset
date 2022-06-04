import sys
import io
import pip
from contextlib import redirect_stdout
from packaging import version
import subprocess

class capture(redirect_stdout):

    def __init__(self):
        self.f = io.StringIO()
        self._new_target = self.f
        self._old_targets = []  # verbatim from parent class

    def __enter__(self):
        self._old_targets.append(getattr(sys, self._stream))  # verbatim from parent class
        setattr(sys, self._stream, self._new_target)  # verbatim from parent class
        return self  # instead of self._new_target in the parent class

    def __repr__(self):
        return self.f.getvalue()


def check_version():
    try:
        p = subprocess.run(['python', '-m','pip', 'index', '--retries', '1', 'versions', 'jetset'], capture_output=True, text=True,timeout=1)
        message = p.stdout
    except:
        p.kill()
        outs, errs = p.communicate()

    lines=str(message).split('\n')
    installed= None
    latest = None
    for line in lines:
        if 'INSTALLED:' in line:
            installed=line.split(':')[1].strip()
        if 'LATEST:' in line:
            latest=line.split(':')[1].strip()
    if installed is not None and latest is not None:
        if  version.parse(installed) < version.parse(latest):
            o_message=('there is a new jetset version %s\n '%latest)
            o_message+=('there installed one is %s',installed)
        else:
            o_message='jetset version is up to date v:%s'%installed
    else:
        o_message='not able to get latest version info from pip'
    
    print(o_message)

def run_version_checking():
    try:
        check_version()
    except Exception as e:
        print("can not establish internet connection, I can't verify if jetset version is up to date")
