Dump various model Hamiltonians as FCIDUMP files.

  *  print help:

            dumpham -h

  *  rewrite FCIDUMP removing symmetry and adding all zero elements back:
  
            dumpham -d <FCIDUMP> [<FCIDUMP.NEW>]
        

  *  rewrite FCIDUMP removing symmetry and adding all zero elements back + create a file with orbital coefficients:
  

            dumpham -do <FCIDUMP> [<FCIDUMP.NEW> [<ORBDUMP>]]
        

  *  use [dumpham-input file](#input-file) to specify the model Hamiltonian, the output files and other parameters:
  

            dumpham <input.dh>


  *  the input parameters can be set from the command line (will be overwritten by the input file), e.g.,
    
  
            dumpham -i "hubbard,U=2.0" -i "ham,out=hubbard.fcidump" hubbard.dh


## Input file


The input file is a text file with the following format:

    % comment line
    set1,name1=value1, name2=value2,...
    set2,name1=value1, name2=value2,...
    ...
    \bham
      <hamiltonian specification>
    \eham

`set`s, `name`s, and default values together with descriptions can be found in the `params.reg` file.

The hamiltonian specification can either be the name of the FCIDUMP file,

    \bham
    H2O.FCIDUMP
    \eham
or several FCIDUMP files,
      
    \bham
    \output{adddump.FCIDUMP}
    \input{H2O.FCIDUMP} + \input{H2O_ADD.FCIDUMP}
    \eham     
or the specification of the [model Hamiltonian](#model-hamiltonians), e.g.,

    \bham
    \output{hubbard.FCIDUMP}
    \geom{dimension=10, pbc=1}
    \hubbard{U=1.0, t={0.5,0.25}}
    \eham

See `*.dh` files in the `test` directory for various examples of input files.

## Model Hamiltonians

For all model Hamiltonians one has to specify the geometry of the system using the `\geom` command. The geometry is specified by the dimension of the system and the periodic boundary conditions. The dimension is the number of sites in each direction. The periodic boundary conditions are specified by a list of 0s and 1s, where 0 means open boundary condition and 1 means periodic boundary condition. The length of the list must be equal to the dimensionality of the system. For example, the following command specifies a 2D system with 10 sites and the periodic boundary condition in the `x` direction, and 2 sites and open boundary condition in the `y` direction:

    \bham  
    \geom{dimension={10,2}, pbc={1,0}}
    \eham

More complicated geometries can be specified by defining coordinates of all sites within a unit-cell and lattice vectors. For example, the following command specifies fused benzene rings system with the periodic boundary condition in `y` direction (and the supercell of 3 rings):

    \geom{dimension={1,3}, pbc={0,1},
          lat={{5.2,0}{0,4.503332}},
          ucell={{1.3,0}{3.9,0}{0,2.251666}{5.2,2.251666}{1.3,4.503332}{3.9,4.503332}} }

The supercell can be stored in a `xyz` file, e.g., for visualization,

    periodic,xyzout=model.xyz
    ...
    \bham
    ...
    \eham

### Hubbard model

The Hubbard model 
```math
\hat H = -t \sum_{\langle i,j \rangle, \sigma} \hat a^\dagger_{i,\sigma} \hat a_{j,\sigma} + U \sum_{i,\sigma,\rho} \hat a^\dagger_{i\sigma} \hat a^\dagger_{i\rho} \hat a_{i\rho} \hat a_{i\sigma}
```
is specified by the `\hubbard` command. The following parameters can be specified:

    \bham
    \geom{...}
    \hubbard{U=2.0, t=1.0}
    \eham

* `U` - the on-site Coulomb repulsion
* `t` - the hopping parameters. The hopping parameters can be specified as a single number (nearest-neighbour hopping) or as a list of numbers (to include next-nearest neighbour hoppings etc)

### Heisenberg model

The Heisenberg model
```math
\hat H = J \sum_{\langle i,j \rangle} \hat S_i \hat S_j
```
is specified by the `\heisenberg` command. The following parameters can be specified:

    \bham
    \geom{...}
    \heisenberg{j=2.0, k=-10.0, norbs=2}
    \eham

* `j` - the nearest-neighbour exchange coupling
* `k` - on-site penalty exchange term
* `norbs` - the number of orbitals per site (= local spin times two)

### Pariser-Parr-Pople model

The Pariser-Parr-Pople model
```math
\hat H = -t \sum_{\langle i,j \rangle, \sigma} \hat a^\dagger_{i,\sigma} \hat a_{j,\sigma} + U \sum_{i,\sigma,\rho} \hat a^\dagger_{i\sigma} \hat a^\dagger_{i\rho} \hat a_{i\rho} \hat a_{i\sigma}
 + \sum_{i\lt j,\sigma,\rho} \frac{U}{\sqrt{1+ar_{ij}^2}} \hat a^\dagger_{i\sigma} \hat a^\dagger_{j\rho} \hat a_{j\rho} \hat a_{i\sigma}
```
is specified by the `\ppp` command. The following parameters can be specified:

    \bham
    \geom{...}
    \ppp{U=2.0, t=1.0, a=0.5}
    \eham

* `U` - the on-site Coulomb repulsion
* `t` - the hopping parameters. The hopping parameters can be specified as a single number (nearest-neighbour hopping) or as a list of numbers (to include next-nearest neighbour hoppings etc)
* `a` - the parameter of the long-range Coulomb repulsion
