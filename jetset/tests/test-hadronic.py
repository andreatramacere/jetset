from jetset.jet_model import  Jetpp
j=Jetpp(proton_distribution='plc')
j.set_verbosity(0)
j.parameters.gmin.val=2
j.parameters.gmax.val=1E8
j.parameters.NH_pp.val=1E10
j.parameters.N.val=1E1
j.parameters.B.val=80

j.parameters.p.val=2.5
j.eval()
j.show_model()
j.plot_model()
