# Library: AutoIO
[//]: # (Badges)
[![Anaconda-Server Badge](https://anaconda.org/auto-mech/autoio/badges/version.svg)](https://anaconda.org/auto-mech/autoio)

Andreas V. Copan, Kevin B. Moore III, Sarah N. Elliott, and Stephen J. Klippenstein

Input writing, output parsing, and job submission tools created as part of the AutoMech package. This repository contains the following python packages: elstruct, autoparse, autoio, and autorun.

<hr size=20>

## Package: elstruct

### Installation
```python
>>> conda install elstruct -c auto-mech
```
### Description
The elstruct package handles input writing and output reading for electronic structure jobs. Current modules available:
- Single point calculation: CFOUR(v2), Gaussian(v09, v16), Molpro(v2015), Orca(v4), PSI4
- Geometry calculation: CFOUR(v2), Gaussian(v09, v16), Molpro(v2015), Orca(v4), PSI4
- Gradients calculation: CFOUR(v2), Gaussian(v09, v16), Molpro(v2015), Orca(v4), PSI4
- Frequency calculation: CFOUR(v2), Gaussian(v09, v16), Molpro(v2015), Orca(v4), PSI4
- VPT2 anharmonic analysis: Gaussian(v09, v16)
- IRC: Gaussian(v09, v16), PSI4


### Usage
Our pytest tests serve as an example for building filesystems

<hr>

## Package: autoparse

### Installation
```python
>>> conda install autoparse -c auto-mech
```
### Description
Set of wrappers for python.re that allow for more intelligible regular expressions syntax.

### Usage
Our pytest tests serve as an example for building filesystems

<hr>

## Package: autoio

### Installation
```python
>>> conda install autoio -c auto-mech
```
### Description
AutoIO generates the input and reads the output for  C++ and Fortran codes used by the AutoMech suite. 
Interfaces for codes distributed alongside AutoMech:
- [MESS](https://github.com/Auto-Mech/MESS): Master Equation System Solver
    The scope of this master equation solver code is described in [Reformulation and Solution of the Master Equation for Multiple-Well Chemical Reactions](https://pubs.acs.org/doi/10.1021/jp4060704).
    ```python
    >>> conda install mess -c auto-mech
    ```
- [ProjRot](https://github.com/Auto-Mech/ProjRot): Rotational Projections
    Projects rotational frequencies from the hessian and generates the resulting projected frequencies
    ```python
    >>> conda install projrot -c auto-mech
    ```
- [ThermP](https://github.com/Auto-Mech/ThermP): Thermochemistry Properties
    Converts partition function to thermochemical proprties
    ```python
    >>> conda install thermps -c auto-mech
    ```
- [PAC99](https://github.com/Auto-Mech/PAC99): Nasa Polynomials
    Converts thermochemical properties into a polynomial representation
    ```python
    >>> conda install pac99 -c auto-mech
    ```
- [PAC99](https://github.com/Auto-Mech/PAC99): Nasa Polynomials
    Converts thermochemical properties into a polynomial representation
    ```python
    >>> conda install pac99 -c auto-mech
    ```
- [OneDMin](https://github.com/Auto-Mech/OneDMin): 
- [VareCof](https://github.com/Auto-Mech/OneDMin): 
- Polyrate: 
- 
External Codes
- CHEMKIN: i

### Usage
Our pytest tests serve as an example for building filesystems

<hr>

## Package: autorun

### Installation
```python
>>> conda install autorun -c auto-mech
```
### Description
AutoRUN submits jobs and extracts relevant output through calls to autoio for  C++ and Fortran codes used by the AutoMech suite.
Interfaces for codes distributed alongside AutoMech:
- [MESS](https://github.com/Auto-Mech/MESS): Master Equation System Solver
    The scope of this master equation solver code is described in [Reformulation and Solution of the Master Equation for Multiple-Well Chemical Reactions](https://pubs.acs.org/doi/10.1021/jp4060704).
    ```python
    >>> conda install mess -c auto-mech
    ```
- [ProjRot](https://github.com/Auto-Mech/ProjRot): Rotational Projections
    Projects rotational frequencies from the hessian and generates the resulting projected frequencies
    ```python
    >>> conda install projrot -c auto-mech
    ```
- [ThermP](https://github.com/Auto-Mech/ThermP): Thermochemistry Properties
    Converts partition function to thermochemical proprties
    ```python
    >>> conda install thermps -c auto-mech
    ```
- [PAC99](https://github.com/Auto-Mech/PAC99): Nasa Polynomials
    Converts thermochemical properties into a polynomial representation
    ```python
    >>> conda install pac99 -c auto-mech
    ```
- [PAC99](https://github.com/Auto-Mech/PAC99): Nasa Polynomials
    Converts thermochemical properties into a polynomial representation
    ```python
    >>> conda install pac99 -c auto-mech
    ```
- [OneDMin](https://github.com/Auto-Mech/OneDMin): 
- [VareCof](https://github.com/Auto-Mech/OneDMin): 
- 
### Usage
Our pytest tests serve as an example for building filesystems.
