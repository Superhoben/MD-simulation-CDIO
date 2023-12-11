.. Computational Physics Project documentation master file, created by
   sphinx-quickstart on Mon Dec  4 13:06:56 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Computational Physics Project's documentation!
=========================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

This is the documentation page for the Computational Pysics Project in the course TFYA99 at Linköping University.

Run the program
---------------
| To start the molecular dynamics program, run the following command in the terminal: 

```
python3 User_interface/user_interface.py
```

Tests
=====
In order to verify the simulations, as well as all implemented properties, a number of tests were performed. These tests include unit-tests which were performed on every individual function that calculates a specific property of the system, as well as larger system tests, which ran simulations and tested all implemented properties. In order to verify that the results were reasonable, the calculated values were compared to reference values from Physics Handbook and other online sources. It was made sure that the simulations reached equillibrium. Below are the results from the larger system tests.

NVT Simulations
-----------
Tests were run with the following simulation parameters

.. list-table:: NVT settings
   :widths: 50 50
   :header-rows: 1

   * - Setting name
     - Chosen setting
   * - Ensemble
     - NVT
   * - Temperature
     - 500
   * - Potential
     - EMT
   * - Step number
     - 10000
   * - Time step
     - 1
   * - Friction
     - 0.005
   * - Basic properties interval
     - 50
   * - Displacement properties interval
     - 50
   * - Basic properties interval
     - 50
   * - Elastic properties interval
     - 500
   * - Configuration interval
     - 500
   * - Optimal scaling interval
     - 50


For copper, Cu, the following results were obtained

.. list-table:: NVT Copper properties
   :widths: 33 33 33
   :header-rows: 1

   * - Property
     - Calculated value
     - Reference value
   * - Average pressure [GPa]
     - -0.005
     - 0
   * - Bulk modulus [GPa]
     - 130
     - 137
   * - Shear modulus [GPa]
     - 63
     - 48.3
   * - Youngs modulus [GPa]
     - 163
     - 129.8
   * - Poisson ratio
     - 0.29
     - 0.35
   * - Mean square displacement
     - 0.12
     - Low
   * - Lindemann coefficient
     - 0.135
     - Low (Below 0.1 preferably)
   * - Diffusion coefficient
     - 5e-5
     - Low
   * - Lattice constant [Å]
     - 3.47
     - 3.61

As we can see below, the total energy is pretty much constant.

.. image:: images/NVT_Cu_validation_test_energy.png
  :width: 500
  :align: center

The temperature for the NVT simulation is aroud 500 K throughout the simulation, which it should be.

.. image:: images/NVT_Cu_validation_test_temperature.png
  :width: 500
  :align: center

The internal preassure during the simulation is around 0 GPa, which is the case for optimized volume.

.. image:: images/NVT_Cu_validation_test_pressure.png
  :width: 500
  :align: center


For silver, Ag, the following results were obtained

.. list-table:: NVT Silver properties
   :widths: 33 33 33
   :header-rows: 1

   * - Property
     - Calculated value
     - Reference value
   * - Average pressure [GPa]
     - -0.006
     - 0
   * - Bulk modulus [GPa]
     - 94.8
     - 103
   * - Shear modulus [GPa]
     - 38.9
     - 30.3
   * - Youngs modulus [GPa]
     - 102.3
     - 82.7
   * - Poisson ratio
     - 0.319
     - 0.43
   * - Mean square displacement
     - 0.25
     - Low
   * - Lindemann coefficient
     - 0.171
     - Low (Below 0.1 preferably)
   * - Diffusion coefficient
     - 0.0001
     - Low
   * - Lattice constant [Å]
     - 4.05
     - 4.08


NVE Simulations
-----------
Tests were run with the following simulation parameters

.. list-table:: NVE Copper settings
   :widths: 50 50
   :header-rows: 1

   * - Setting name
     - Chosen setting
   * - Ensemble
     - NVE
   * - Temperature
     - 500
   * - Potential
     - EMT
   * - Step number
     - 10000
   * - Time step
     - 1
   * - Friction
     - 0.005
   * - Basic properties interval
     - 50
   * - Displacement properties interval
     - 50
   * - Basic properties interval
     - 50
   * - Elastic properties interval
     - 500
   * - Configuration interval
     - 500
   * - Optimal scaling interval
     - 10


For copper, Cu, the following results were obtained

.. list-table:: NVE Copper properties
   :widths: 33 33 33
   :header-rows: 1

   * - Property
     - Calculated value
     - Reference value
   * - Average pressure [GPa]
     - 0
     - 0
   * - Bulk modulus [GPa]
     - 130
     - 151
   * - Shear modulus [GPa]
     - 63
     - 57
   * - Youngs modulus [GPa]
     - 163
     - 110
   * - Poisson ratio
     - 0.29
     - 0.35
   * - Mean square displacement
     - 0.12
     - Low
   * - Lindemann coefficient
     - 0.135
     - Low (Below 0.1 preferably)
   * - Diffusion coefficient
     - 5e-5
     - Low
   * - Lattice constant [Å]
     - 3.47
     - 3.61


For silver, Ag, the following results were obtained

.. list-table:: NVT Silver properties
   :widths: 33 33 33
   :header-rows: 1

   * - Property
     - Calculated value
     - Reference value
   * - Average pressure [GPa]
     - 0.002
     - 0
   * - Bulk modulus [GPa]
     - 163.4
     - 217
   * - Shear modulus [GPa]
     - 34
     - 27
   * - Youngs modulus [GPa]
     - 94
     - 78.5
   * - Poisson ratio
     - 0.4
     - 0.52
   * - Mean square displacement
     - 0.208
     - Low
   * - Lindemann coefficient
     - 0.157
     - Low (Below 0.1 preferably)
   * - Diffusion coefficient
     - 8.8e-5
     - Low
   * - Lattice constant [Å]
     - 3.921
     - 4.078

As we can see below, the total energy in the NVE simulation is close to constant.

.. image:: images/NVE_Au_validation_test_energy.png
  :width: 500
  :align: center

The temperature is aroud 500 K throughout the simulation.

.. image:: images/NVE_Au_validation_test_temperature.png
  :width: 500
  :align: center

The internal preassure during the simulation is around 0 GPa, which is the case for optimized volume.

.. image:: images/NVE_Au_validation_test_pressure.png
  :width: 500
  :align: center

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
