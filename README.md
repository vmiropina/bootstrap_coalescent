# bootstrap_coalescent

This repository contains the code used in the artilce "A non-parametric estimator of the coagulation measure of Lambda-coalescents" [1], by Verónica Miró Pina, Emilien Joly and Arno Siri-Jégousse.

It uses msprime [2,3] to simulate the Site Frequency Spectrum of a population under a beta-coalescent model. 
Then, the coagulation measure Lambda is estimated using the techniques suggested in [1].

The code runs in python3 and the dependencies are: numpy, scipy, matplotlib and msprime. To install msprime see: https://tskit.dev/msprime/docs/stable/installation.html



References: 

[1] A non-parametric estimator of the coagulation measure of Lambda-coalescents, Verónica Miró Pina, Emilien Joly and Arno Siri-Jégousse. To appear

[2] Franz Baumdicker, Gertjan Bisschop, Daniel Goldstein, Graham Gower, Aaron P Ragsdale, Georgia Tsambos, Sha Zhu, Bjarki Eldon, E Castedo Ellerman, Jared G Galloway, Ariella L Gladstein, Gregor Gorjanc, Bing Guo, Ben Jeffery, Warren W Kretzschumar, Konrad Lohse, Michael Matschiner, Dominic Nelson, Nathaniel S Pope, Consuelo D Quinto-Cortés, Murillo F Rodrigues, Kumar Saunack, Thibaut Sellinger, Kevin Thornton, Hugo van Kemenade, Anthony W Wohns, Yan Wong, Simon Gravel, Andrew D Kern, Jere Koskela, Peter L Ralph and Jerome Kelleher (2022), Efficient ancestry and mutation simulation with msprime 1.0, Genetics, Volume 220, Issue 3. http://doi.org/10.1371/journal.pcbi.1004842

[3] Jerome Kelleher, Alison M Etheridge and Gilean McVean (2016), Efficient Coalescent Simulation and Genealogical Analysis for Large Sample Sizes, PLOS Comput Biol 12(5): e1004842. doi: 10.1371/journal.pcbi.1004842
