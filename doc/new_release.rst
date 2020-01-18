What's new in version 1.1.1
===========================

This release fixes some minor bugs, and fixes a problem with swig for mac os mojave for the Jet.energetic_report method

The QuickStart notebook has been improved.


The conda-forge channel has ben removed becase was installing on mac os mojave a bugged version of matplotlib.

In the following are listed the main new features common to 1.1.0 and 1.1.1


#. Astropy Tables with units have been implemented for most of the products.



#. Improved model serialization.



#. Improved External Compton implementation with possibility to use a double approach

   * transformation of the external  fields to the blob rest frame :cite:`Dermer2000`

   *  transformation of the electron emitting distribution from the blob restframe to
      disk/BH restframe :cite:`Dermer95` :cite:`GKM01`


#. Implementation of Broad Line Region radiative field using the approach of :cite:`Donea2003`



.. important::
    Starting from version 1.1.0, the `R` parameter of the :class:`.jet_model.Jet` class, as default *is linear and not logarithmic*, please update your old scripts
    setting `R` with linear values (see :ref:`jet_physical_guide` for more details).

.. important::
    starting from version 1.1.0 the saved model format has changed, if you have models saved vith version<1.1.0,
    plase update them the new models by loading the old models with the :meth:`.jet_model.Jet.load_old_model`
    and then saving them again (see :ref:`jet_physical_guide` for more details).


.. bibliography:: references.bib
