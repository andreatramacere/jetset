import  pytest
from .base_class import TestBase

class TestData(TestBase):

    def integration_suite(self,plot=False):
        self._all(plot=plot)

    def test_data_loader(self,plot=True,sed_number=2):
        from jetset.data_loader import ObsData, Data
        from jetset.test_data_helper import test_SEDs

        data = Data.from_file(test_SEDs[0])
        sed_data = ObsData(data_table=data)
        sed_data.show_data_sets()
        sed_data.filter_data_set('-1', exclude=True)
        sed_data.show_data_sets()
        sed_data.reset_data()
        sed_data.show_data_sets()
        sed_data.filter_data_set('-1', exclude=False)
        sed_data.show_data_sets()
        sed_data.reset_data()

        data=Data.from_file(test_SEDs[sed_number])
        sed_data=ObsData(data_table=data)
        sed_data.group_data(bin_width=0.2)
        print('sed file->',test_SEDs[sed_number])
        sed_data.add_systematics(0.1,[10.**6,10.**29])
        if plot is True:
            p=sed_data.plot_sed()

        return sed_data
    

    def test_data_from_numpy(self,plot=False):
        pass

    def test_data_from_table(self,plot=False):
        pass