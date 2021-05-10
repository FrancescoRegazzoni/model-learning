.. _conventions:

====================================
Notes about the conventions
====================================


Function arguments
---------------------------------

Functions typically accept a list of mandatory arguments, briefly documented in the function help. Optional arguments are typically passed by a structure array ``opt``. In the first lines of the function, default values are typically assigned to undefined optional arguments. These lines thus provide a list of the possible optional arguments.

Variable names
---------------------------------

In the code the names of the variables concerning the models refelct the convention used in `[1] <https://doi.org/10.1016/j.jcp.2019.07.050>`_. Specifically, models are denoted as follows:

.. math::
	\frac{d \mathbf{x}}{dt} &= \mathbf{f}(\mathbf{x}(t),\mathbf{u}(t),\boldsymbol{\alpha}) \\
	\mathbf{y}(t) &= \mathbf{g}(\mathbf{x}(t))

where :math:`t` is time, :math:`\mathbf{x}(t)` is the internal state, :math:`\mathbf{u}(t)` (if present) is the time-dependent input, :math:`\boldsymbol{\alpha}` (if present) is a constant parameter, :math:`\mathbf{y}(t)` is the output.