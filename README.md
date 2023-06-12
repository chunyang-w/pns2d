# 🌊 PNS-2D: A paralleled 2-D Navier-Stoke Solver

## Compile the code

```shell
cd build
cmake ..
```

## Execute testing cases

```shell
cd build
ctest
```

## Execute the code

```shell
cd build
mpiexec -n <num_of_processors> ./bin/NS_Solver
```