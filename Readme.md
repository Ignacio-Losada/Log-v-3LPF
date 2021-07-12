

<img src="https://github.com/Ignacio-Losada/Log-v-3LPF/blob/main/log(v)3LPFlogo.png" align="justified" width="1000" alt="Log(v) 3LPF logo">

Log(v)  3LPF is  a  linear power  flow  solver  for  unbalanced  three-phase  distribution  systems.  Log(v)  3LPF  uses  a  logarithmic  transform  of  the  voltage phasor  to  linearize  the  AC  power  flow  equations  around  the balanced   case.   It includes   the   modeling   of   ZIP   loads,transformers, capacitor banks, switches and their corresponding controls.  Log(v) 3LPF  uses  the  Sherman-Morrison-Woodbury  identity  for  an efficient computation of the inverse of a rank-k corrected system admittance matrix matrix, Ybus. In general, the computation of the inverse using Sherman-Morrison Woodbury is more efficient than traditional LU decomposition  methods  in  terms  of  FLOPS.  This tool has been validated against OpenDSS on different  network  sizes,  ranging  from  tens to  thousands  of  nodes.




## Installation (Python)
Use Git to install Log(v) 3LPF in your current python environment
```bash
git clone https://github.com/Ignacio-Losada/Log-v-3LPF
cd .\Log-v-3LPF\
pip3 install -r requirements.txt

```

## Installation with Git (Conda)
Use Git to install a conda environment for Log(v) 3LPF
```bash
git clone https://github.com/Ignacio-Losada/Log-v-3LPF
cd SoDa
conda env create -f environment.yml

```

## Basic Usage
Once Log(v) 3LPF is installed, a power flow that uses any .dss test case within openDSS or any network with a similar format can be run as follows




## Citing Log(V) 3LPF

If you find Log(v) 3LPF useful in your work, we kindly request that you cite the following [publication]():
```
@article{carreno2021log,
  title={Log (v) 3LPF: A Linear Power Flow Formulation for Unbalanced Three-Phase Distribution Systems},
  author={Carre{\~n}o, Ignacio Losada and Saha, Shammya and Scaglione, Anna and Arnold, Daniel and Sy-Toan, Ngo and Roberts, Ciaran},
  year={2021},
  publisher={TechRxiv}
}
```


## License

