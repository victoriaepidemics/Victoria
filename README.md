![Logo](https://www.cimat.mx/sites/all/themes/danland/logo.png)
![Logo](https://www.matem.unam.mx/++theme++im-theme-blue/images/logo_imunam.png)

# Victoria Epidemics

## Overview

The ongoing COVID-19 pandemic has posed a major challenge to public health systems of all countries with the imminent risk of saturated hospitals, resulting in patients not receiving proper medical care. However, having the possibility of estimating the excess hospital care demand in advance, it is possible to mitigate this consequences, improving thus patients' health conditions and saving more lives.  

Inspired by this, we develop the Victoria Epidemic suite in order to forecast hospital occupancy in metropolitan areas during, not only the current COVID-19 outbreak, but also any other infectious disease that can follow a SEIRD type model. 

A SEIRD model is an epidemiological compartmental model which divides the population into groups, called Susceptible, Exposed, Infected, Recovered and Dead, whose evolution in time and relations are described by an Ordinary Differential Equations System. It is useful for describing hospital dynamics and, in our case, Bayesian inference is used to calibrate some key features of the model given some measured data: using both hospital admittance confirmed cases and deaths we infer the contact rate and the initial conditions of the dynamical system. 

Having a basic background in mathematical statistics, and even without any knowledge of object-oriented programming, you can use our tool to model the dynamics of some infectious diseases or to simply change some parameters to explore possible scenarios. It is important to mention that our model has been used by the federal government of Mexico during the current COVID-19 pandemic to assist public policy, and has been applied for the analysis of more than 70 metropolitan areas and the 32 states in the country.


## Installing

You can install the package via Git:

```shell
python -m pip install git+<this repo>
```
or you can download the tarball [victoriaepi-0.0.1.tar.gz](dist/victoriaepi-0.0.1.tar.gz) and run:

```shell
python -m pip install victoriaepi-0.0.1.tar.gz
```
or you can simply download the code, `cd` into the directory containing `setup.py` and then run:

```shell
pip install .
```

## Getting started

After installing the package, download the file [AnalyzeEyam.py](victoriaepi/examples/AnalyzeEyam.py) or copy the following code to a file:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 09:27:13 2020

@author: jac

SEID model analysis of the Bubonic plague outbreak in Eyam, England

Data from the second black plague outbreak in the village of Eyam,
England, from 19 June to 1 November 1666 (Stojkova and D. A. Campbell 2017).
"""

import os
from datetime import date
from victoriaepi.seid import SEID
import matplotlib.pyplot as plt

workdir = "./"

# Residence rates 1/day, names and Erlang series
R_rates = \
        {'E':[1/5       , r'\sigma_1',  4],
         'I':[1/4      , r'\sigma_2',  3]}

trim=-30
eyam = SEID( Region="Eyam, England, Bubonic plague outbreak 1,666",\
                N=261, data_fnam="Eyam.csv", out_fnam="Eyam",\
                init_index=0, init=date( 1666, 6, 1), R_rates=R_rates, trim=trim, workdir= workdir)

os.makedirs(workdir+"output", exist_ok=True)

T=2000
pred=-trim+7*3
if T>0:
    eyam.RunMCMC( T=T, burnin=1000, pred=pred, plot_fit=False)

eyam.PlotEvolution(pred=pred, cumm=True, right_axis=False)
plt.show()
```
Create a subfolder named `data`, place [Eyam.csv](victoriaepi/examples/data/Eyam.csv) into it and then
run:

```shell
 $ python AnalyzeEyam.py
```

If everything worked right, you should see something like:


```text
$ python AnalyzeEyam.py
File with mcmc samples does not exist, run RunMCMC first.
S --> E
E --> I
I --> D
pytwalk: Running the twalk with 2000 iterations .  Sun, 04 Oct 2020, 13:57:25.
           Finish in approx. 1 min and 15 sec.
pytwalk: finished, Wed, 25 Nov 2020, 13:57:52.
AutoMaxlag: maxlag= 72.

Effective sample size: 41
Sampling 42 model solutions.

Saving files in  ./output/Eyam_*.pkl
```
### AMA model

If you just managed to run the previous example, you can try this one.


Download the file [AnalysisMZ.py](victoriaepi/examples/AnalysisMZ.py),
create a subfolder named `data` and place [9-01.csv](victoriaepi/examples/data/9-01.csv)
and [9-01_DinHosp.csv](victoriaepi/examples/data/hosp/9-01_DinHosp.csv) as follows:

```text
.
├── AnalysisMZ.py
└── data
    ├── 9-01.csv
    └── hosp
        └── 9-01_DinHosp.csv
```

Then, run:

```shell
 $ python AnalysisMZ.py
```

If everything worked right, you should have this new structure:


```text
.
├── AnalysisMZ.py
├── csv
│   ├── 9-01_D_long.csv
│   ├── 9-01_D_short.csv
│   ├── 9-01_Hs.csv
│   ├── 9-01_I_long.csv
│   ├── 9-01_I_short.csv
│   └── 9-01_U1.csv
├── data
│   ├── 9-01.csv
│   └── hosp
│       └── 9-01_DinHosp.csv
├── figs
│   ├── 9-01_D_long.png
│   ├── 9-01_D_short.png
│   ├── 9-01_Hs.png
│   ├── 9-01_I_short.png
│   ├── 9-01_Omega_f.png
│   ├── 9-01_Omega.png
│   ├── 9-01_R.png
│   └── 9-01_U1.png
└── output
    ├── 9-01_samples.pkl
    ├── 9-01_solns.pkl
    └── 9-01_solns_plain.pkl
```
 


## Features

Project highlights:
* A SEIRD model to set up and change the dynamics of some infectious diseases exploring different scenarios
* A Bayesian inference approach to calibrate some key features of the model given some measured data
* The possibility to handle intervention dates
* The possibility to handle age groups
* Usage of configuration files

## Contributors

This model has been developed by:
* [Andrés Christen](https://coronavirus.conacyt.mx/investigadores/jacg.html), CIMAT-CONACYT
* [Marcos Capistrán](https://coronavirus.conacyt.mx/investigadores/maco.html), CIMAT-CONACYT
* [Antonio Capella](https://coronavirus.conacyt.mx/investigadores/ack.html), Institute of Mathematics-UNAM 
* [Judith Esquivel](https://coronavirus.conacyt.mx/investigadores/jev.html), CIMAT-CONACYT
* [Oscar Miguel González](https://coronavirus.conacyt.mx/investigadores/jev.html), CIMAT-CONACYT


## Contributing

If you'd like to contribute, please fork the repository and use a feature
branch. Pull requests are warmly welcome.

## Links

* Project homepage: https://coronavirus.conacyt.mx/proyectos/ama.html
* Paper homepage: https://arxiv.org/abs/2006.01873
* Additional information: https://www.cimat.mx/~jac/twalk/

## Licensing

The code in this project is licensed under 3-Clause BSD License.
