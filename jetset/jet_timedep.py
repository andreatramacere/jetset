from __future__ import absolute_import, division, print_function

from builtins import (str, open, super, range,
                      object, map)




__author__ = "Andrea Tramacere"



import os

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd == True:
    try:
        from jetkernel import jetkernel as BlazarSED
    except ImportError:
        from .mock import jetkernel as BlazarSED
else:
    from jetkernel import jetkernel as BlazarSED


from .output import makedir,WorkPlace


from .utils import safe_run,set_str_attr, old_model_warning


__all__=['TempEvol']

class TempEvol(object):

    def __init__(self,out_dir='temp_ev',flag='tests',clean_work_dir=True):

        self.build_TempEv()

        self.set_path(out_dir, clean_work_dir=clean_work_dir)

        self.set_flag(flag)

    def build_TempEv(self,duration=1E5,
                     TStart_Acc=0.,
                     TStop_Acc=0.,
                     TStart_Inj=0.,
                     TStop_Inj=0.,
                     T_esc_Coeff=1E60,
                     Esc_Index=0.,
                     Acc_Index=1.,
                     Diff_Index=2.,
                     T_SIZE=5000,
                     NUM_SET=50,
                     Lambda_max_Turb=1E30):




        self._temp_ev = BlazarSED.MakeTempEv()

        self._temp_ev.duration = duration
        self._temp_ev.TStart_Acc = TStart_Acc
        self._temp_ev.TStop_Acc = TStop_Acc
        self._temp_ev.TStart_Inj = TStart_Inj
        self._temp_ev.TStop_Inj = TStop_Inj
        self._temp_ev.T_esc_Coeff = T_esc_Coeff
        self._temp_ev.Esc_Index = Esc_Index
        self._temp_ev.Acc_Index = Acc_Index
        self._temp_ev.Diff_Index = Diff_Index
        self._temp_ev.T_SIZE = T_SIZE
        self._temp_ev.NUM_SET = NUM_SET
        self._temp_ev.Lambda_max_Turb=Lambda_max_Turb



    def set_path(self, path, clean_work_dir=True):
        if path.endswith('/'):
            pass
        else:
            path += '/'

        set_str_attr(self._temp_ev, 'path', path)
        # set_str(self._blob.path,path)
        makedir(path, clean_work_dir=clean_work_dir)



    def set_flag(self,flag):
        self._temp_ev.STEM=flag


    def run(self,jet):
        BlazarSED.Run_temp_evolution(jet._blob, self._temp_ev)