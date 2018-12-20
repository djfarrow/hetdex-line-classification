# Line Classification Tool

## Description

A tool to classify detected emission lines as OII or Lyman-alpha. Designed
for the HETDEX survey. Based on an idea and code from Leung et al 2017,
i.e. [this paper](http://adsabs.harvard.edu/abs/2017ApJ...843..130L). The
method described in that paper has been reimplemented in classification\_prob\_leung. 
A new approach, that deals with errors on the equivalent width and involves fewer integrals
is also included classification\_prob and will be presented in Farrow et al (in prep).

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

Currently Python 3 is not supported. The new stuff in ```classification_prob.py`` tested.

## References

Farrow et al (in prep)
[Leung et al 2017](http://adsabs.harvard.edu/abs/2017ApJ...843..130L)
