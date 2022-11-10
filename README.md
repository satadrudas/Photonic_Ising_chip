# Photonic_Ising_chip

A time-multiplexed Coherent Ising Machine chip is designed based of Silicon Photonics. Basically here I implement the table-top Ising machine design given in [1] on a simple chip. Though with some changes. The nonlinear function used for bifurcation in [1] is cosine squared since the spin information is encoded in the intensity of light. Here we use just cosine function for bifurcation and encode the spin information in the amplitude of the light. For the matrix multiplication part, and method photonic matrix multiplication using homodyne-detector given in [2] is used since this system is a time multiplex system.



References: <br />
1. [A poor man’s coherent Ising machine based on opto-electronic feedback systems for solving optimization problems](https://www.nature.com/articles/s41467-019-11484-3)<br />
2. [Large-Scale Optical Neural Networks Based on Photoelectric Multiplication](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.9.021032)<br />
