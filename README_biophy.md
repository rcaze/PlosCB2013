Author: Romain Cazé 
Date: 26/07/2013

===============
About this code
===============

This work is licensed under the Creative Commons Attribution Licence (CCAL). 
Under the CCAL, the author allows anyone to download, reuse, modify this code.
vecstim.mod has not been written by Romain Cazé but the author is unknown.

This code originates from Cazé, R. D., Humphries, M., & Gutkin, B. (2013). Passive Dendrites Enable Single Neurons to Compute Linearly Non-separable Functions. 
(O. Sporns, Ed.)PLoS Computational Biology, 9(2), e1002867. doi:10.1371/journal.pcbi.1002867"
It has been modified and dramatically compressed after publication.
The aim is to increase clarity and reusability, it was the first step to generate plots on Figure 3.C and Figure.4.
Only the core of the script has been kept which is the simulation class.
The BipolarNeuron class can be replaced by an arbitrarly complex neuron model.

This code snippet demonstrates under noiseless condition that scattered synaptic are more likely to make a neuron fires than clustered inputs.
This is what we call the scattering sensitivity, our article demonstrate that this sensitivity can be used to compute linearly non-separable functions.

More generally, this snippet enables to stimulate and record the response of a given neuron model using arbitrarly complex 
spike times list and where the position of the synapses can be precisely controlled, all within PYTHON (see instert_signal method).

=====
Files
=====

plos2013_biopĥy.py
vecstim.mod

============
Requirements
============

This snippet works with Python 2.6, Matplotlib 1.1.1, and Numpy 1.3 or greater.
It also requires NEURON simulator (http://www.neuron.yale.edu) version 7.2 or greater.

======
How to 
======

First run nrnivmodl in a terminal, within the folder containing vecstim.mod, to compile it.
Then run the plos2013_biophy.py. 
The default protocol is showing that a passive model is more sensitive to scattered than clustered inputs. 

===========
Contribute!
===========

Please report any bugs, make recommendations, and add your own source or help
enhance this library by contacting romain`dot`caze`at`gmail`dot`com. 



