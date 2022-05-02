import pytest
import inspect
class TestBase:

    def _all(self,**kwargs):
        m_list=inspect.getmembers(self, inspect.ismethod)
        for m in m_list:
            if m[0].startswith('test_') and  m[0]!='test_all':
                m[1](kwargs)    

    def integration_suite(self):
        raise RuntimeError('this method has to be implemented in each subclass')