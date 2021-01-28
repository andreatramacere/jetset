from jetset.jet_model import Jet
def test_my_foo():
    from jetset.plot_sedfit import  plt
    plt.ioff()
    j = Jet()
    j.eval()
    j.energetic_report()
