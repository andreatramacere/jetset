leastsqbound
============

What is leastsqbound
--------------------

leastsqbound is a enhanced version of scipy's optimize.leastsq function which
allows users to include min, max bounds for each fit parameter. Constraints are 
enforced by using an unconstrained internal parameter list which is
transformed into a constrained parameter list using non-linear functions.
