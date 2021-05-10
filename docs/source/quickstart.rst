====================================
Quickstart
====================================

The repository contains several examples of usage of the library, some of which are documented and described in `[1] <https://doi.org/10.1016/j.jcp.2019.07.050>`_. 
Each example is contained into a dedicated folder in ``\examples`` (e.g. ``\examples\wave_oned``). The standard structure for an :ref:`example<example>` envisages:

- The :ref:`problem<problem>` specifications, through an ``.ini`` file (e.g. ``\examples\wave_oned\wave1D.ini``).

- The high-fidelity :ref:`model<model>` definition, through a ``.m`` file (e.g. ``\examples\wave_oned\wave1D_getmodel.ini``).

- A script that generates training and testing :ref:`datasets<dataset>` (e.g. ``\examples\wave_oned\wave1D_generate_tests.ini``) and stores them into the ``datapath`` folder.

- One or more configuration files to train a reduced model, with the format ``opt_*.ini`` (e.g. ``\examples\wave_oned\opt_wave1D.ini``). To use them, execute the command::

	model_learn('opt_wave1D.ini')

- One or more scripts to analyze and postprocess the results obtained with the trained models and the HF model.

To create a new :ref:`example<example>`, create a new folder in the folder ``\examples``. 

This documentation also contains a :ref:`Tutorial<tutorial>` to get acquainted with the library. 