:orphan:

Configparsing Subpackage
==============================

The purpose of this subpackage is to create usable code ready for simulation from human-readable configuration files.

This package intends to facilitate the creation and modification of ODE models for different diseases.

If you want to learn how to use this subpackage jump directly to the :ref:`Usage <configparsing usage>` section.



Modules
^^^^^^^

.. image:: imgs/salida_configparsing.png

Call graph showing the interrelation of the modules in the package.


Main module
------------------------------------
.. automodule:: victoriaepi.configparsing.config
   :members: ModelConfig
   :undoc-members:
   :show-inheritance:
   :special-members: __repr__
   :autosummary: 

  

Call module
----------------------------------------

.. automodule:: victoriaepi.configparsing.modelcalls
   :members: ModelCall, Rrates, Zone, rrate, strevaldate
   :undoc-members:
   :show-inheritance:
   :special-members: __repr__
   :autosummary:
    

Class module
----------------------------------------

.. automodule:: victoriaepi.configparsing.modelclass
   :members: EFractions, ModelMatrix, functionCollection, functiondef, plottingMethods, rawfuncbody, solvingMethods
   :undoc-members:
   :show-inheritance:
   :special-members: __repr__
   :autosummary:
    

configparsing.generic module
-------------------------------------

.. automodule:: victoriaepi.configparsing.generic
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __repr__
   :autosummary:





