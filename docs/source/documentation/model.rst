.. _model:

.. highlight:: matlab

.. index:: ! model

============================================================
The struct ``model``
============================================================

What is a model?
----------------

A :ref:`model<model>` can be regarded as a map from a time-dependent input to a time-dependent output, according to the definition of Sec. 2.1 of `[1] <https://doi.org/10.1016/j.jcp.2019.07.050>`_. Each model is associated to a :ref:`problem<problem>`. Notice that a model comprises both the mathematical model and its discretization: it corresponds thus to the concept of **numerical model**.

How is a model represented?
---------------------------

A model is represented as a `MATLAB struct <https://www.mathworks.com/help/matlab/ref/struct.html>`_. The struct must contain a field ``model.problem``, containing the associated :ref:`problem struct<problem>`.

Models can be of two types: black-box models or explicitly-defined models.

Black-box models
^^^^^^^^^^^^^^^^
Black-box models are characterized by the field ``model.blackbox = 1``. The input-output map is defined through a function handle ``model.forward_function``, with the following signature: ::

	output = forward_function(test, options)

This function receives an input :ref:`test struct<test>` and returns an output :ref:`test struct<test>`. An example of definition of such function is contained in ``\examples\tutorial\NTL1_getmodel_blackbox.m``.

Black-box models are useful when one wants to build a wrapper to an external software.


Explicitly-defined models
^^^^^^^^^^^^^^^^^^^^^^^^^
Explicitly-defined models contain the explicit definition of the properties of the model. In this case, the model struct must contain the following fields:

**Basic fields**

- ``model.nX``: number of internal states

- ``model.x0``: initial state

- ``model.dt``: integration time step

**Model dynamics**

The model dynamics can be defined in different ways, according to the field ``model.advance_type``, that takes the following values:

- ``'solveonestep'``: the most generic case. It requires a function ``x = model.solveonestep(x,u,dt)`` that advances the state of one time step.

- ``'solveonestep_linear'``: like ``'solveonestep'``, but this entails the solution of a linear system :math:`Kx = r`. It requires a function ``[K,r] = model.solveonestep(x,u,dt)``.

- ``'nonlinear_explicit'``: explicit Euler advance, with generic right-hand side, defined by the function ``rhs = model.f(x,u)``.

- ``'linear_CN_uaffine'``: generic linear model, with Crank-Nicolson advance, with affine dependence on :math:`\mathbf{u}(t)`.

- ``'linear_advance'``: generic linear model.

- ``'linear_advance_timeexplicit'``: generic linear model, with explicit dependence on :math:`\Delta t`.

- ``'linear_advance_uaffine'``: generic linear model, with affine dependence on :math:`\mathbf{u}(t)`.

- ``'linear_advance_timeexplicit_uaffine'``: generic linear model, with explicit dependence on :math:`\Delta t` and with affine dependence on :math:`\mathbf{u}(t)`.

For the definition of the fields required according to the field ``model.advance_type``, see the function ``core/model_solve.m``

**Model output**

Also the output of the model can be defined in different ways, according to the field ``model.output_type``, that takes the following values:

- ``'insidestate'``: the output is given by the first ``nY`` internal states (i.e. ``y = x(1:nY)``).

- ``'nonlinear'``: the output is a non-linear function of the state: ``y = model.g(x)``.

- ``'linear'``: the output is a non-linear function of the state: ``y = model.G * x + model.g0``.


.. note::

	The first case corresponds to the *input-inside-the-state approach* of `[1] <https://doi.org/10.1016/j.jcp.2019.07.050>`_, while the other two correspond to the *input-outside-the-state approach*.


How to create or load a model struct?
-------------------------------------

The best-practice is that of creating models by means of ad-hoc functions. Examples are contained in ``\examples\tutorial\NTL1_getmodel_blackbox.m`` (black-box case) and ``\examples\tutorial\NTL1_getmodel.m`` (explicitly-defined case). 

Such model constructors, that receive in input the :ref:`problem struct<problem>`, can be used in the ``.ini`` file defining a problem to define the high-fidelity model handler. In this case, the high-fidelity model associated with a problem can be obtained, for example, with the folowing commands: ::

   problem = problem_get('tutorial','NTL1.ini');
   HFmod = problem.get_model(problem);

Also reduced models based on ANNs are **models**, in the sense of the definition of this page. To load a trained model, you can use the following function: ::

	ANNmod = read_model_fromfile(problem,'FOLDER_OF_TRAINED_MODEL');


What can I do with a model struct?
-------------------------------------

Models structs can be used to perform simulations. With the following command, for instance, we employ the model defined in the struct ``model`` to obtain the output :ref:`test struct<test>` ``test_output`` associated with the input :ref:`test struct<test>` ``test_input`` and we plot the solution: ::

   figure();
   test_output = model_solve(test_input, model, struct('do_plot',1));