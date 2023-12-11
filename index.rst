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
python User_interface/user_interface.py
```

Tests
=====
In order to verify the simulations, as well as all implemented properties, a number of tests were performed. These tests include unit-tests which were performed on every individual function that calculates a specific property of the system, as well as larger system tests, which ran simulations and tested all implemented properties. In order to verify that the results were reasonable, the calculated values were compared to reference values from Physics Handbook and other online sources. Below are the results from the larger system tests.

NVT Copper
-----------
This test was run with the following simulation parameters

.. list-table:: Title
   :widths: 50 50
   :header-rows: 5

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

The following results were obtained

.. list-table:: Title
   :widths: 33 33
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


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
