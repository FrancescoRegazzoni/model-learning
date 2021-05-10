.. _example:

.. highlight:: matlab

.. index:: ! example

============================================================
The folder ``example``
============================================================

What is an example?
-------------------------------

In this library, an **example** is simply a logical container that groups the files corresponding to a apecific application of the library. It corresponds to a folder contained in ``/examples``. For instance, each of the test cases considered in `[1] <https://doi.org/10.1016/j.jcp.2019.07.050>`_ has its own example folder. 

How to create a new example?
-------------------------------

If you want to develop your own application with the library ``model-learning``, you simply have to create a folder inside ``/examples``. 

.. note::

	Inside ``/examples``, folders starting with the prefix ``app_`` (the standard prefix for custom applications) are automatically gitignored.