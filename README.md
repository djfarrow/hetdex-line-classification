# Line Classification Tool

## Description

A tool to classify detected emission lines as OII or Lyman-alpha. Designed
for the HETDEX survey. Based on an idea and code from Leung et al 2017,
i.e. [this paper](http://adsabs.harvard.edu/abs/2017ApJ...843..130L). The
method described in that paper has been reimplemented in `classification_prob_leung`. 
A new approach, that deals with errors on the equivalent width and involves fewer integrals
is also included in `classification_prob` and will be presented in Farrow et al (in prep).

## Setup

To install, including packages needed for tests

```
pip install .[test]
```

if you want to run the tests then do

```
pytest
```

if you do not want to run tests you can just do

```
pip install .
```

## Known issues

1. Currently Python 3 is not supported. 
2. Colors aren't used in the probability yet.  

## References

If you use this code in a paper, please cite

[Leung et al 2017](http://adsabs.harvard.edu/abs/2017ApJ...843..130L)

Farrow et al (in prep)

If you use the configuration file in the test directory you'll also need to cite the following. It uses 
relative line strengths from

[Anders et al 2003](http://adsabs.harvard.edu/abs/2003A%26A...401.1063A)

and luminosity/equivalent width functions from

[Ciardullo et al 2013](http://adsabs.harvard.edu/abs/2013ApJ...769...83C)
Gronwall et al (in prep)

and cosmology from

[Planck 2013](https://ui.adsabs.harvard.edu/#abs/arXiv:1303.5076)



