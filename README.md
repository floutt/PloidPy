# PloidPy

## Introduction
PloidPy is a program written in Python designed to infer ploidy from next-generation reads aligned to a haploid reference genome. The program makes use of a stastical model representing the distribution of specific nucleotide counts and selects the most probable ploidy on the basis of a minimum AIC.

## Installation
Installation of PloidPy is relatively simple and can be done easily using pip. The most recent version can be installed using the following command:
```
pip install git+git://github.com/floutt/PloidPy
```
### Dependencies
In order to run PloidPy, the following dependencies are required:
- `Python 3.6+`
- `NumPy`
- `SciPy`
- `Statsmodels`
- `matplotlib`
- `seaborn`
