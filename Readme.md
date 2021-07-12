# Log(v) 3LPF

<img src="https://github.com/Ignacio-Losada/Log-v-3LPF/blob/main/log(v)3LPF-logo.png" align="right" width="500" alt="Log(v) 3LPF logo">

Log(v)  3LPF is  a  linear power  flow  solver  for  unbalanced  three-phase  distribution  systems.  Log(v)  3LPF  uses  a  logarithmic  transform  of  the  voltage phasor  to  linearize  the  AC  power  flow  equations  around  the balanced   case.   It includes   the   modeling   of   ZIP   loads,transformers, capacitor banks, switches and their corresponding controls.  Log(v) 3LPF  uses  the  Sherman-Morrison-Woodbury  identity  for  an efficient computation of the inverse of a rank-k corrected system admittance matrix matrix, Ybus. In general, the computation of the inverse using Sherman-Morrison Woodbury is more efficient than traditional LU decomposition  methods  in  terms  of  FLOPS.  This tool has been validated against OpenDSS on different  network  sizes,  ranging  from  tens to  thousands  of  nodes.




## Installation (Python)
Use Git to install Log(v) 3LPF in your current python environment
```bash

```

## Installation with Git (Conda)
Use Git to install a conda environment for pysoda
```bash

```

## Basic Usage
Once Log(v) 3LPF is installed....


## Citing Log(V) 3LPF

If you find Log(v) 3LPF useful in your work, we kindly request that you cite the following [publication]():
```

```


## License

