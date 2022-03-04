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

The library provides a Python interface to deploy the trained model. The module can be loaded as: ::

	import pyModelLearning
	
To install it, you can choose among the following options.	


Option 1: Install the library by `setuptools <https://setuptools.readthedocs.io/>`_
------------------------------------------------------------------------------------------

From the main folder of this repository run: ::

    $ pip install . --use-feature=in-tree-build


Option 2: Add the library path to Python paths
------------------------------------------------------------------------------------------

Linux / macOS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add to your ``.bashrc``: ::

    export PYTHONPATH="${PYTHONPATH}:/path/to/model-learning"

Windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add the path of ``model-learning`` to the ``PYTHONPATH`` environment variable.
Alternatively, you can write the path of ``model-learning`` inside a ``pth`` file, as described `here <https://docs.python.org/3/using/windows.html#finding-modules>`_.

Option 3: Add the library path within your python script
------------------------------------------------------------------------------------------

Put the following lines at the beginning of your Python script: ::

    import sys
    sys.path.append("/path/to/model-learning")
    import pyModelLearning

