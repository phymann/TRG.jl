# TRG.jl

Classical two-dimensional Ising model in a longitudinal field, given by the following Hamiltonian,
$$H = -J \sum_{\langle i,j \rangle} \sigma_i^z \sigma_j^z - h \sum_i \sigma_i^z$$

is solved numerically using the TRG algorithm proposed in [PRL **99**, 120601 (2007)](https://link.aps.org/doi/10.1103/PhysRevLett.99.120601).

## Usage

Start from `runTRG.jl`