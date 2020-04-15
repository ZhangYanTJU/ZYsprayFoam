# ZYsprayFoam

[![OpenFOAM version](https://img.shields.io/badge/OpenFOAM-7-brightgreen)](https://github.com/OpenFOAM/OpenFOAM-7)

sprayFoam with liquid and vapor penetration output.

Original code is from [OpenFOAM-dev sprayFoam](https://github.com/OpenFOAM/OpenFOAM-dev/tree/master/applications/solvers/lagrangian/sprayFoam).

```
functions
{
    #includeFunc Qdot
    #includeFunc writeObjects(RR.C7H16)
    #includeFunc writeObjects(h)
}
```