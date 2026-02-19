# GENIE-INCL user manual

## Installation

*Note: This installation is based on `/cvmfs/larsoft.opensciencegrid.org` at Fermilab, you need to have fnal `gpvm` account. (Please let me know if you would like to install GENIE-INCL in another place.)*

- Clone feature branch from github and source setup scripts
    ```
    $ git clone -b feature/incl-tune-beta git@github.com:LiangLiu212/Generator.git
    $ source setup.sh # This command may take some time, as it needs to load modules via Spack.
    $ source do_configure.sh
- Running `make` to compile
    ```
    $ make -j <N>
- Now you are ready to run GENIE-INCL simulation
    ```
    $ goto parent folder and list
    $ cd  ../
    $ ls 
    Generator  inclxx_genie  runarea  setup_env.sh
    $ cd runarea
    $ ./run_gevgen.sh INCL26_01a_00_000
- After compiling, you only need to `source setup_env.sh` the next time you log in.
## GENIE-INCL Event Record

### Here is an example of GENIE-INCL event record

```
|------------------------------------------------------------------------------------------------------------------|
|GENIE GHEP Event Record [print level:   3]                                                                        |
|------------------------------------------------------------------------------------------------------------------|
| Idx |          Name | Ist |        PDG |   Mother  | Daughter  |      Px |      Py |      Pz |       E |      m  |
|------------------------------------------------------------------------------------------------------------------|
|   0 |         nu_mu |   0 |         14 |  -1 |  -1 |   4 |   4 |   0.000 |   0.000 |   2.000 |   2.000 |   0.000 |
|   1 |           C12 |   0 | 1000060120 |  -1 |  -1 |   2 |   3 |   0.000 |   0.000 |   0.000 |  11.175 |  11.175 |
|   2 |       neutron |  11 |       2112 |   1 |  -1 |   5 |   5 |  -0.208 |  -0.077 |   0.062 |   0.920 | **0.940 | M = 0.891
|   3 |           C11 |   2 | 1000060110 |   1 |  -1 |  19 |  19 |   0.208 |   0.077 |  -0.062 |  10.255 |  10.254 |
|   4 |           mu- |   1 |         13 |   0 |  -1 |  -1 |  -1 |  -0.188 |   0.610 |   1.489 |   1.623 |   0.106 | P = (0.116,-0.377,-0.919)
|   5 |        proton |  14 |       2212 |   2 |  -1 |   7 |   8 |  -0.020 |  -0.687 |   0.573 |   1.297 |   0.938 | FSI = 1
|   6 |       neutron |  19 |       2112 |  -1 |  -1 |   7 |   8 |   0.091 |  -0.164 |  -0.071 |   0.960 |   0.940 | FSI = 1
|   7 |       neutron |  14 |       2112 |   5 |   6 |  16 |  16 |  -0.142 |  -0.148 |  -0.007 |   0.961 |   0.940 | FSI = 0
|   8 |        proton |  14 |       2212 |   5 |   6 |  10 |  11 |   0.213 |  -0.713 |   0.518 |   1.305 |   0.938 | FSI = 1
|   9 |        proton |  19 |       2212 |  -1 |  -1 |  10 |  11 |  -0.020 |   0.035 |   0.168 |   0.954 |   0.938 | FSI = 1
|  10 |        proton |  14 |       2212 |   8 |   9 |  13 |  14 |   0.068 |  -0.462 |   0.745 |   1.286 |   0.938 | FSI = 1
|  11 |        proton |  14 |       2212 |   8 |   9 |  17 |  17 |   0.126 |  -0.216 |  -0.060 |   0.973 |   0.938 | FSI = 0
|  12 |        proton |  19 |       2212 |  -1 |  -1 |  13 |  14 |   0.085 |   0.048 |   0.057 |   0.945 |   0.938 | FSI = 1
|  13 |        proton |  14 |       2212 |  10 |  12 |  15 |  15 |   0.347 |  -0.203 |   0.710 |   1.243 |   0.938 | FSI = 0
|  14 |        proton |  14 |       2212 |  10 |  12 |  18 |  18 |  -0.194 |  -0.210 |   0.092 |   0.985 |   0.938 | FSI = 0
|  15 |        proton |   1 |       2212 |  13 |  -1 |  -1 |  -1 |   0.351 |  -0.206 |   0.718 |   1.250 |   0.938 |
|  16 |       neutron |  14 |       2112 |   7 |  -1 |  -1 |  -1 |   0.195 |  -0.060 |   0.027 |   0.961 |   0.940 |
|  17 |        proton |  14 |       2212 |  11 |  -1 |  -1 |  -1 |   0.049 |   0.237 |  -0.087 |   0.973 |   0.938 |
|  18 |        proton |   1 |       2212 |  14 |  -1 |  -1 |  -1 |  -0.045 |  -0.049 |   0.021 |   0.941 |   0.938 |
|  19 |           B10 |  17 | 1000050100 |   3 |  -1 |  20 |  24 |  -0.118 |  -0.356 |  -0.229 |   9.403 | **9.324 | M = 9.393
|  20 |        proton |   1 |       2212 |  19 |  -1 |  -1 |  -1 |   0.077 |   0.014 |  -0.031 |   0.942 |   0.938 |
|  21 |            H2 |   1 | 1000010020 |  19 |  -1 |  -1 |  -1 |  -0.027 |   0.029 |  -0.084 |   1.878 |   1.876 |
|  22 |            H2 |   1 | 1000010020 |  19 |  -1 |  -1 |  -1 |   0.179 |  -0.024 |  -0.111 |   1.887 |   1.876 |
|  23 |           He4 |   1 | 1000020040 |  19 |  -1 |  -1 |  -1 |  -0.282 |  -0.332 |   0.014 |   3.753 |   3.727 |
|  24 |       neutron |   1 |       2112 |  19 |  -1 |  -1 |  -1 |  -0.064 |  -0.044 |  -0.018 |   0.943 |   0.940 |
|------------------------------------------------------------------------------------------------------------------|
|       Fin-Init:                                                |   0.000 |  -0.000 |  -0.000 |   0.042 |         |
|------------------------------------------------------------------------------------------------------------------|
```


The `Ist` is the status of particle, you can find the code in `$GENIE/src/Framework/GHEP/GHepStatus.h`

| Ist   | enum             | explanation             |
|-------|------------------|-------------------------|
| 0     | kIStInitialState | initial state           |
| 1     | kIStStableFinalState | final state           |
| 2     | kIStIntermediateState | intermediate state           |
| 14     | kIStHadronInTheNucleus | hadron propagate inside nucleus           |
| 15     | kIStFinalStateNuclearRemnant | final state nuclear remnant (A>4)           |
| 17     | kIStPreDeExNuclearRemnant | excited nuclear remnant           |
| 19     | kIStSpectator | the spectator particle in FSI collision           |

When a particle propagates inside nucleus, it will collide with a spectator. So, after collision, particles inside nucleus will have two parents.

