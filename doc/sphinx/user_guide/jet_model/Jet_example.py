
# coding: utf-8

# # jet_model module user guide


In this section we describe how to use the module  :mod:`.jet_model`
to buil a model of jet able to reproduce SSC/EC emission processes. 


=============================================The :mod:`.jet_model`  allows to build a jet  model  providing an interface 
to call the BlazarSED code. The BlazarSED code is a numerical 
accurate C code, to evaluate SSC/EC emission processes in a relativistic jet. 
The python wrappper is  built using SWIG. 

A jet can be built using the  the :class:`.Jet` class, istanciating a jet object.

# In[41]:

import BlazarSEDFit as SEDFit
myJet=SEDFit.Jet('test','lppl')

This instruction will create a ``Jet`` object with ``name`` **test**,
and for an electron distribution **lppl** that is a log-parabola with
a low-energy power-law branch.
The parameters of the model arre accessible throug the instruction
# In[42]:

myJet.parameters.show_pars()


Each parameter has default values. All the parameters listed are handled by
:class:`.ModelParameterArray`, and each parameter is an instance of the the
:class:`.JetParameter`.
# 
# At this point one can evaluate the SSC/EC emission for this jet model using the instruction

# In[43]:

myJet.eval()


# and plot the corresponding SED:
# 

# In[44]:

myPlot=SEDFit.Plot()

myPlot.add_model_plot(myJet,autoscale=True)
myPlot.save('jet.png')


# To change one of the parameter in the model:
# 

# In[45]:


myPlot=SEDFit.Plot()

myPlot.add_model_plot(myJet,autoscale=True)
myJet.set_par('gamma0_log_parab',val=1.0E5)

myJet.eval()

myPlot.add_model_plot(myJet,label='gamma0_log_parab=1E5',autoscale=True)

myPlot.save('jet1.png')


# To plot all the components

# In[46]:

myPlot=SEDFit.Plot()
for c in myJet.spectral_components: myPlot.add_model_plot(c.SED,autoscale=True)


# In[ ]:




# In[ ]:




# In[ ]:



