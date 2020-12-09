
.. _configparsing usage:

Using the configuration parser
******************************
This part explains how to process configuration files using Configparsing Subpackage. We strongly recommend getting acquainted with the :ref:`Extended modules <Extended modules>` section before attempting to create your own configuration file.


.. _Getting started:
Getting started
===============

Three basic configuration files are provided as example, each one corresponding to recreate the  :ref:`extension modules <Extended modules>`:

* :download:`SEIDconfig.hjson  <../victoriaepi/examples/SEIDconfig.hjson>`
* :download:`AMAconfig.hjson <../victoriaepi/examples/AMAconfig.hjson>`
* :download:`AMA_GRPconfig.hjson <../victoriaepi/examples/GRPconfig.hjson>`

The file :download:`exampleconfig.py <../victoriaepi/examples/exampleconfig.py>` shows how to use the configuration files.

.. literalinclude:: ../victoriaepi/examples/exampleconfig.py

After each execution of the `writecode()` method :meth:`victoriaepi.configparsing.config.ModelConfig.writecode` two files will be created as per the `moduleName` directive located inside the configuration file :any:`victoriaepi.configparsing.config.ModelConfig.moduleName`.

For example, for the SEID configuration file, these could be ``autoseid.py`` and ``call_autoseid.py``; where the former includes the model's class and the latter can be executed via ``python call_autoseid.py`` to get the results.

Please note that you still need to provide the data file as per the ``data_fnam`` directive in the configuration file, and you can expect the results in the supplied ``workdir``.



Config file structure
=====================
The configuration files are written in the `Hjson` file format. Hjson is a syntax extension to JSON and it is intended to be used as a user interface for humans, to read and edit before passing the JSON data to machine [#f2]_. 

In the next sections, we will analyze the structure of a config file with the highest grade of complexity which is the ``grp`` see the :ref:`Complexity <Complexity>` section.

General structure
-----------------
The config  files are hierarchical, and the parameters found in the root structure are directly passed to :any:`victoriaepi.configparsing.config.ModelConfig`.

Let's start with the following configuration. We have trimmed some child members to facilitate viewing:

.. code-block:: python


    {
    
    complexity: grp
    className: ama2_grp
    moduleName: autoGRP
    #version control
    baseModelVersion: "1.2.3"
    ModelMatrix:{
            #defining forward map
            prnDefault:False
            names: "S E I^A I^S I^B H^1 H^2 H^3 U^1 U^2 D R"
            BaseVar: S
            Exits:{}
            ExitFractions:{}
        }
    num_pars : 3 # Number of parameters to be inferred
    factor_foi: 0.1*20
    #model call
    ModelCall:{  }
    odeint: "integrate.odeint( rhs, X0, t_quad, args)"

    solvingMethods:{
        rhs:{}
        solve_plain1:{}
        solve_plain2:{ }
        solve:{}
        llikelihood:{}
        lprior:{}
        support:{}
        sim_init:{}
        }

    plottingMethods:{
        PlotEvolution:{}
        PlotPanelEvolution:{}
        PlotFit:{}
        PlotStateVar:{}
        StateSlice:{}
        PlotStateSlice:{}
        PlotOmega:{}
        PlotPanelEvolution:{}
        }

    additionalMethods:{
    #Mainly used for overriding
        SetTime:{}
        GetMask:{}
        PrintConfigParameters:{}
        }

    additionalCallMethods:{}


From here you can see that all the parameters directly translate to the main attributes from the class which are:

* :any:`victoriaepi.configparsing.config.ModelConfig.complexity`
* :any:`victoriaepi.configparsing.config.ModelConfig.className`
* :any:`victoriaepi.configparsing.config.ModelConfig.moduleName`
* :any:`victoriaepi.configparsing.config.ModelConfig.baseModelVersion`
* :any:`victoriaepi.configparsing.config.ModelConfig.ModelMatrix`
* :any:`victoriaepi.configparsing.config.ModelConfig.num_pars`
* :any:`victoriaepi.configparsing.config.ModelConfig.factor_foi`
* :any:`victoriaepi.configparsing.config.ModelConfig.ModelCall`
* :any:`victoriaepi.configparsing.config.ModelConfig.odeint`
* :any:`victoriaepi.configparsing.config.ModelConfig.solvingMethods`
* :any:`victoriaepi.configparsing.config.ModelConfig.plottingMethods`
* :any:`victoriaepi.configparsing.config.ModelConfig.additionalMethods`
* :any:`victoriaepi.configparsing.config.ModelConfig.additionalCallMethods`

In the next sections we will review some important aspects from these attributes.


.. _Complexity:
Complexity
----------
If you look at the above-mentioned config files (:ref:`Getting started <Getting started>`) you will find that they have a different value of a parameter called ``complexity``. This is an internal parameter that helps the parser to locate 
the required templates, as each one of these models has different requirements for plotting and data handling.
The three built-in complexity values are:

* ``basic`` for simple ODE models with MCMC Bayesian parameter inference.
* ``intervention`` includes ``basic`` but adds the possibility of adding intervention dates.
* ``grp`` all of the above plus it is able to handle age groups differently.

Of course, it is possible to use the ``grp`` value for simpler tasks but, using the right `complexity` value will output a leaner code, better suited to that particular problem.

You can define your own ``complexity`` value by modifying the files found under the ``calls`` and ``models`` directories and giving them a new name. 


.. code-block:: text


    victoriaepi/configparsing
    ├── calls
    │   ├── basic.py
    │   ├── grp.py
    │   └── intervention.py
    └── models
        ├── buildingblocks
        │   ├── basic_ModelMatrix.py
        │   ├── grp_ModelMatrix.py
        │   └── intervention_ModelMatrix.py
        ├── basic.py
        ├── grp.py
        └── intervention.py


ModelMatrix
-----------


Defines the forward map and has mainly information that goes to writing the class file.

.. code-block:: python

    ModelMatrix:{
            #defining forward map
            prnDefault:False
            names: "S E I^A I^S I^B H^1 H^2 H^3 U^1 U^2 D R"
            BaseVar: S
            Exits:{
              #AuxMatrix.Exit( 'S', 'E')
              "S" : {"E":{}}

              #AuxMatrix.SplitExit( 'E', 'I^A', 'I^S', 1-f, prob_symb=['1-f','f'])
              "E" : {
                  "I^A":{
                      "prob":1-f
                      "prob_symb":'1-f'
                      }
                  "I^S":{
                      "prob":f
                      "prob_symb":'f'}
                      }

              #AuxMatrix.Exit( 'I^A', 'R')
              "I^A" : {"R":{}}

              #AuxMatrix.SplitExit( 'I^S', 'H^1', 'I^B', g, prob_symb=['g', '1-g'])
              "I^S" : {
                  "H^1":{
                      "prob":g
                      "prob_symb":'g'}
                  "I^B":{
                      "prob":1-g
                      "prob_symb":'1-g'
                      }
                      }

              #AuxMatrix.Exit( 'I^B', 'R')
              "I^B" : {"R":{}}

              #AuxMatrix.SplitExit( 'H^1', 'U^1', 'H^2', h, prob_symb=['h', '1-h'])
              "H^1" : {
              "U^1":{"prob":h
                      "prob_symb":'h'}
                      "H^2":{"prob":1-h
                      "prob_symb":'1-h'}}

              #AuxMatrix.Exit( 'H^2', 'R')
              "H^2" : {"R":{}}

              #AuxMatrix.Exit( 'U^1', 'U^2')
              "U^1" : {"U^2":{}}

              #AuxMatrix.SplitExit( 'U^2', 'D', 'H^3', i, prob_symb=['i', '1-i'])
              "U^2" : {
              "D":{"prob":i
                      "prob_symb":'i'}
                      "H^3":{"prob":1-i
                      "prob_symb":'1-i'}}

              #AuxMatrix.Exit( 'H^3', 'R')
              "H^3" : {"R":{}}

              #AuxMatrix.NoExit('R')
              "R" : {}
              #AuxMatrix.NoExit('D')
              "D" : {}
            }
            ExitFractions:{
                #Each can be an expression or a List of expressions
                # fraction of severe infections
                f:[0.05, 0.05, 0.05, 0.05]
                # fraction of severe infections that require hospitalization
                g:[0.004     , 0.0348    , 0.12333333, 0.2396    ]
                # fraction of hospitations that require ICU
                h:[0.05      , 0.0552    , 0.17266667, 0.5112    ]
                # fraction of ICU patients who die
                i:[0.5, 0.5, 0.5, 0.5]
            }
        }


Is processed by :any:`victoriaepi.configparsing.modelclass.ModelMatrix`. You can check individual parameters:

* :any:`victoriaepi.configparsing.modelclass.ModelMatrix.prnDefault`
* :any:`victoriaepi.configparsing.modelclass.ModelMatrix.names`
* :any:`victoriaepi.configparsing.modelclass.ModelMatrix.BaseVar`
* :any:`victoriaepi.configparsing.modelclass.ModelMatrix.Exits`
* :any:`victoriaepi.configparsing.modelclass.ModelMatrix.ExitFractions`

The Exit Fractions :any:`victoriaepi.configparsing.modelclass.ModelMatrix.ExitFractions` will later processed as :class:`victoriaepi.configparsing.modelclass.EFractions`.
See also :any:`victoriaepi.victoria.AuxMatrix` for more information on how these values are used.


ModelCall
---------
Is processed by :any:`victoriaepi.configparsing.modelcalls.ModelCall`.

Mainly parameters to create the file that makes the instance and executes the simulation:

.. code-block:: python

    ModelCall:{
            Region:"CDMX COVID-19"
            #instanceName: seidinstance
            # population size
            N:261
            burnin:1000
            # Data file name, should be in workdir + 'data/' + data_fnam
            data_fnam:"9-01.csv"
            out_fnam:"auto_9-01_grp" # MCMC output file name, without .txt, workdir/output/out_fnam + '.txt'
            init_index:0 # day number from data start, where to start the plot
            init_date:"date(2020, 4, 1)"# date of init_index
            trim:0,# how many data points to trim
            pred: 200
            plotpred: 28
            workdir:"./../amafromconfig/"
            T:10000
            plot_fit:False
            # Contact matrix
            Int_M: np.ones((ngrp,ngrp))/ngrp #Uniform
            #probability of recording an infection
            Pobs_I:0.85
            # probability of recording a death
            Pobs_D:0.95
            #Time to symptoms (days)
            daysdelay: 4
            # Residence rates 1/day, names and Erlang series
            R_rates:{
                'E':{
                    residenceRate:1/5
                    name:\sigma_1
                    erlang:4
                    }
                'I^S':{
                    residenceRate:1/4
                    name:\sigma_2
                    erlang:3
                   }
                'H^1':{
                    residenceRate:1/1
                    name:\sigma_3
                    erlang:1
                    }
                'U^1':{
                    residenceRate:1/5
                    name:\sigma_4
                    erlang:3
                    }
                'U^2':{
                    residenceRate:1/1
                    name:\sigma_5
                    erlang:1
                    }
                'I^A':{
                    residenceRate:1/7
                    name:\gamma_1
                    erlang:3
                    }
                'I^B':{
                    residenceRate:1/7
                    name:\gamma_2
                    erlang:3
                    }
                'H^2':{
                    residenceRate:1/6
                    name:\gamma_3
                    erlang:3
                    }
                'H^3':{
                    residenceRate:1/3.3
                    name:\gamma_4
                    erlang:1
                    }
            }
            # General information for the metro zone or region to be analyzed
            Zone:{
                id: "clave9-01"
                name: "Mexico city"
                num_relax_days: 2
                population: 21.942666e6
                init_date: date(2020, 2, 27)
                intervention_n_relax_dates: [ "date(2020, 3, 22)", "date(2020, 4, 3)", "date(2020, 5, 10)", "date(2020, 6, 3)"]
            }
            plottingQuantiles: [10,25,50,75,90]
            ngrp: 4
            age_prop: [0.45213552, 0.36863423, 0.11800457, 0.06122568]
            age_groups: [0, 25, 50, 65, 100]
        }


You can check individual parameters:

* :any:`victoriaepi.configparsing.modelcalls.ModelCall.Region` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.N` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.burnin` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.data_fnam` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.out_fnam` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.init_index` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.init_date` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.trim` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.pred` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.plotpred` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.workdir` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.T` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.plot_fit` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.Int_M` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.Pobs_I` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.Pobs_D` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.daysdelay` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.R_rates` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.Zone` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.plottingQuantiles` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.ngrp` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.age_prop` 
* :any:`victoriaepi.configparsing.modelcalls.ModelCall.age_groups` 



Method Groups
--------------

These are dictionaries filled with method definitions:

* ``solvingMethods`` which is processed by :any:`victoriaepi.configparsing.config.ModelConfig.solvingMethods` and expects some specific methods to be defined.
* ``plottingMethods`` which is processed by :any:`victoriaepi.configparsing.config.ModelConfig.plottingMethods` and expects some specific methods to be defined.
* ``additionalMethods`` which is processed by :any:`victoriaepi.configparsing.modelclass.functionCollection` optional, this is useful for overriding :any:`victoriaepi.victoria.mcmc`'s methods. 
* ``additionalCallMethods`` which is processed by :any:`victoriaepi.configparsing.modelclass.functionCollection` optional, to add additional methods to the call file.

A single method definition gets processed by :class:`victoriaepi.configparsing.modelclass.functiondef` and has the following structure:

.. code-block:: python

    rhs:{
            defn:'rhs( self, x, t, p)'
            body:
            '''
            beta1 = p[3 + np.where(t < self.intervention_day)[0][0]]

            # total number of asymptomatic infections in each group
            I_A = np.array([np.sum(x[self.sel_[grp]] * self.mask_IAs_flat) for grp in range(self.ngrp)])
            I_S = np.array([np.sum(x[self.sel_[grp]] * self.mask_ISs_flat) for grp in range(self.ngrp)])
            # array of length self.ngrps

            #force of infection in each group
            foi = beta1/self.N * (self.Int_M @ I_A + self.factor_foi*self.Int_M @ I_S)

            for grp in range(self.ngrp):
                self.par[self.Ts[grp].ConvertVar('S')] = foi[grp]

            for grp in range(self.ngrp):
                self.rt[self.sel_[grp]] = self.Ts[grp].M @ (x[self.sel_[grp]] * (self.Ts[grp].par_mask @ self.par))

            return self.rt
            '''
            }





.. [#f2] https://hjson.github.io/ 
