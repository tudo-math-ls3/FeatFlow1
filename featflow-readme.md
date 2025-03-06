# FeatFlow 1.3

This repository contains the legacy version of FeatFlow 1.

## System Requirements

The code has been confirmed to compile and run with:
- GFortran 11.4

## Build Configuration

### Platform
The only confirmed working build configuration in the Makefile is:
- `pc64-x86-linux`

### Compiler Flags
The following compiler flags are necessary for successful compilation and execution:
```
-fallow-argument-mismath
-ffpe-trap=invalid,zero,overflow
```

### Optimization
Optimization must be disabled:
```
-O0
```

Please note that compatibility with other compiler versions, platforms, or build configurations has not been verified at this time.

# How to run FF1 on Lido3
```
git clone git@github.com:tudo-math-ls3/FeatFlow1.git
cd FeatFlow1/
module load gcc/11.5.0
chmod 750 bin/*
make -j 4 libs
cd applications/cc2d
./cc2d
make
OR
make -j 4 apps (build all apps)


```
