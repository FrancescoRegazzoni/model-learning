.. _test:

.. highlight:: matlab

.. index:: ! test

============================================================
The struct ``test``
============================================================


What is a test?
-------------------------------

A **test** represents a time-dependent input :math:`\mathbf{u}(t)` in an interval :math:`[0,T]` and, optionally, the corresponding output :math:`\mathbf{y}(t)` and state :math:`\mathbf{x}(t)`. 


How is a test represented?
-------------------------------

A test is represented as a `MATLAB struct <https://www.mathworks.com/help/matlab/ref/struct.html>`_. The time is defined in the field ``test.tt``, while the input, output and state are defined in the fields ``test.uu``, ``test.yy``, ``test.xx``, respectively.

A test can be defined analytically, as in the following example: ::

	test.tt = [0, 10];
	test.uu = @(t) sin(t);

Otherwise, the test can be defined by a set of discrete values, as in the following example: ::

	test.tt = linspace(1,10,200);
	test.uu = sin(test.tt);

In case of vector-valued input (i.e. ``nU > 1`` in the :ref:`problem<problem>` definition), each entries of the input corresponds to a row of ``test.uu`` (the same convention if employed also for ``test.yy`` and ``test.xx``), that is: ::

	test.uu = @(t) [sin(t); cos(t)];

or ::

	test.uu = [sin(test.tt); cos(test.tt)];

In case the value of ``test.uu`` for a time step not contained in ``test.tt`` is needed, a linear interpolation is performed.

.. note::

	Both the *input* and the *output* of a simulation are represented as test structs. For example, with the following line one can obtain the solution corresponding the input :math:`\mathbf{u}(t) = \sin(t)` in the interval :math:`t \in [0,10]` for the model represented by the :ref:`struct<model>` ``model``: ::

		test_input.tt = linspace(1,10,200);
		test_input.uu = sin(test.tt);
		test_output = model_solve(test_input, model);

In this example, both ``test_input`` and ``test_output`` are test structs. However, while the former contains only the fields ``test_input.tt`` and ``test_input.uu``, the latter contains the fields ``test_output.tt``, ``test_output.uu`` and ``test_output.yy``. Additionally, by specifying the option ``save_x`` as follows: ::

	test_output = model_solve(test_input, model, struct('save_x', 1));

the struct ``test_output`` also contains the field ``test_output.xx``.
