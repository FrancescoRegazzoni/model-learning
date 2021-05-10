.. highlight:: matlab

============================================================
Learning models with inter-individual variability 
============================================================


In the example considered in the previous section, all the indivuals in the training set shared the same dynamics and the same initial state. We now consider the case of **inter-individuals variability**, that is when the individuals may differ for either the dynamics (more precisely, we suppose that the dynamics follows a common law, dependent on parameters, which may be be different for different individuals) or the intial state, or both.

We consider in particular the following exponential model:

.. math::
	\frac{dx}{dt} &= \alpha x(t) \\
	x(0) &= x_0

where the output is :math:`x` itself and the initial state :math:`x_0` and the parameter :math:`\alpha` may be different for different individuals.

Notice that such object is not a **model** according to the definition of Sec. 2.1 of `[1] <https://doi.org/10.1016/j.jcp.2019.07.050>`_, since it does not uniquely defines a map from the space of inputs (which is empty in this case) to the space of outputs. Thus, we refer to it as a **metamodel** since it defines rather a family of models.

.. index:: problem

Definig the problem
-------------------------------

The problem is defined in ``\examples\tutorial\expmod.ini``. With respect to the previous example, we have two new fields:

- ``fixed_x0 = 0``: the individuals do not share the same initial state

- ``samples_variability = 1``: the law of evolution may depend on parameters, different for each individual.


.. index:: model

Definig the high-fidelity model
-------------------------------

The model (or better, the metamodel) is defined in ``\examples\tutorial\expmod_getmodel.m``. With respect to the previous example, we notice the following differences:

- Also the number of parameters (``nA``) and their bounds are defined

- The right-hand side depends on the parameter, thus it is assigned into ``f_alpha``

- The derivatives of the right-hand side w.r.t. the state and the parameters are provided. This is not mandatory, unless the model will be later used for data assimilation purposes.

We can now load the model::

	problem = problem_get('tutorial','expmod.ini');
	HFmod = problem.get_model(problem);

Notice that we cannot test directly the model, since it is just a metamodel. First, we have to **particularize** it, i.e. assigning initial condtion (e.g. ``x_0 = 1.1``) and parameter (e.g. ``a = 0.5``)::

	HFmod_particularized = metamodel_particularize(HFmod, 1.1, .5);

Now, we can solve the model as before::

	figure()
	model_solve(struct('tt',[0 1]), HFmod_particularized, struct('do_plot',1))

Otherwise, with the following command, we assign a random initial state and a random parameter, and show the result::

	model_show_example(HFmod);


.. index:: dataset

Generating training datasets
-------------------------------

As before, with the following commands we generate a training dataset. In this case, initial state and parameters will be randomly generated.::

	rng('default')

	opt_gen.do_plot = 1;
	opt_gen.do_save = 1;
	opt_gen.T = 1;
	opt_gen.outFile = 'samples.mat';
	dataset_generate_random(HFmod,100,opt_gen);

Training the ANN
-------------------------------

As before, we define the learning specifications in the file ``\examples\tutorial\expmod_opt.ini``. The main differences are:

- We must specify the number of parameters in the learned model, by ``Model\N_alpha``

- We add some noise, by ``Problem\noise_y``

To train the network, run::

	model_learn('expmod_opt.ini')

Using the trained model
-------------------------------

As before, the learned model can be loaded by::

	ANNmod = read_model_fromfile(problem,'test_int_N1_hlayF3_dof13_2019-02-28_17-51-59');

Notice that the learned model is also a metamodel! To use it, it should first be particularized. Otherwise, one can show a possible evolution of an individual with the command::

	model_show_example(ANNmod);

When the parameter is just one, the original and learned parameters can be compared by::

	model_alpha_plot(ANNmod);

Using data assimilation on the learned model
--------------------------------------------------------------

To test how Extended Kalman Filter performs on a given metamodel and for a given level of noise, the command ``da_test`` can be used. It can be used both with the high-fidelity and learned model, since they are metamodels::

	noise = 1e-3;
	da_test(HFmod, noise);
	da_test(ANNmod, noise);

The following command instead, performs the following steps:

- It randomly generates an indidual and simulates its evolution by means of the high-fidelity model;

- It adds the prescribed amount of noise;

- For the first half of the time span of the simulation, it performs data assimilation through the learned model, to estimate the parameter of the individual and its state;

- It uses those two estimated values to predict the future evolution of the individual and computes the error.

::

	opt_ep.mod_HF = HFmod;
	opt_ep.obs_err = 1e-3;
	opt_ep.pause_each_test = 1;
	opt_ep.do_plot = 1;
	da_estimate_predict(problem,ANNmod,opt_ep);

By pressing ``ENTER``, the process is reapeated for other individuals and the mean error is computed.

With the following command instead, the same process is repeated for different noise levels (without stopping at each individual), and the convergence plot is shown.::

	opt_ep.mod_HF = HFmod;
	opt_ep.obs_err = 10.^-(1:8);
	opt_ep.pause_each_test = 0;
	opt_ep.do_plot = 0;
	opt_ep.n_tests = 20;
	da_estimate_predict(problem,ANNmod,opt_ep);
