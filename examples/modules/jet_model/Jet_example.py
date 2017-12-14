import BlazarSEDFit as SEDFit

"""
create a Jet  object, named "test" , with a "lppl" electron  distribution
"""
myJet=SEDFit.Jet('test','lppl')

myJet.parameters.show_pars()

myJet.eval()

myPlot=SEDFit.Plot()

myPlot.add_model_plot(myJet,autoscale=True)

myPlot.save('jet.png')

"""
chante the value of the 'gamma0_log_parab' parameter
"""
myJet.set_par('gamma0_log_parab',val=1.0E5)

myJet.eval()

myPlot.add_model_plot(myJet,label='gamma0_log_parab=1.0E5')

myPlot.save('jet1.png')