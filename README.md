# ISOMIP+ with MOM6
The goal of this repository is to centralize the tests, issues/difficulties and steps needed during the implementation of ISOMIP+ using MOM6 and other model components from GFDL. 

## TESTS

Quiet runs:

* **[2D without ice shelf](https://github.com/gustavo-marques/ISOMIP/blob/master/ipynb/quiet_2D_noIS.ipynb)** 

* **[2D with ice shelf](https://github.com/gustavo-marques/ISOMIP/blob/master/ipynb/quiet_2D_yesIS.ipynb)** 

* **[3D without ice shelf](https://github.com/gustavo-marques/ISOMIP/blob/master/ipynb/quiet_3D_noIS.ipynb)**

* **[3D with ice shelf](https://github.com/gustavo-marques/ISOMIP/blob/master/ipynb/quiet_3D_yesIS.ipynb)**  

Sponge layer:

* **[2D, cold interior and cold forcing](https://github.com/gustavo-marques/ISOMIP/blob/master/ipynb/)** 

* **[2D, cold interior and warm forcing](https://github.com/gustavo-marques/ISOMIP/blob/master/ipynb/)** 

Passive tracers:

Coming soon!

## ISSUES

* **[Spurious velocities in "quiet experiments"](https://github.com/gustavo-marques/ISOMIP/blob/master/ipynb/spurious_velocities_quiet_experiments.ipynb)** 

* **[Horizontal pressure gradient errors due to the ice shelf](https://github.com/gustavo-marques/ISOMIP/blob/master/ipynb/horizontal_pressure_gradient_errors_due_to_the_ice_shelf.ipynb)** 

* **[Inconsistency in ALE sponge layer](https://github.com/gustavo-marques/ISOMIP/blob/master/ipynb/inconsistency_ale_sponge_layer.ipynb)**

## TO DO LIST
* **[Parameterization of vertical mixing under the ice shelf](ipynb/vertical_mixing_parameterization.ipynb)**
 
* **[~~Implementing sponge layer in layer mode~~](https://github.com/gustavo-marques/ISOMIP/blob/master/ipynb/sponge_layer_in_layer_mode.ipynb)**

* **[Impose an evaporative flux in the open ocean region to compensate the melt water](https://github.com/gustavo-marques/ISOMIP/blob/master/ipynb/evaporative_flux.ipynb)**

* **[Create a script to generate the requested output](ipynb/requested_output.ipynb)**
