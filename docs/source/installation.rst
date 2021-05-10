.. _installation:

============
Installation
============

To install the library, follow these steps.

- In the main library folder copy the file ``options_example.ini`` into ``options.ini`` (not git tracked) and customize the field ``datapath`` (path where generated data will be stored). By default, data are stored in the library folder.
- Run ``addpaths.m``.
- Enjoy.


.. note::

	A collection of pre-trained ANNs is available in the repo `model-learning_data <https://github.com/FrancescoRegazzoni/model-learning_data>`_. To use them, clone or download the repo and make the ``datapath`` field point to the path of the repo (or alternatively copy the content of the repo into the ``datapath`` folder).


==============
Python wrapper
==============

To use the Python wrapper, add to your ``.bashrc`` file the following line: ::

	export PYTHONPATH="${PYTHONPATH}:/path/to/model-learning"

Then, you can import the python wrapper as ::

	import pyModelLearning
