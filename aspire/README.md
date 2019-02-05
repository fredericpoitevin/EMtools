# ASPIRE
this repository stores our current implementation of ASPIRE, so we are all on the same page.

## Installation:
1. download and install [ASPIRE](http://spr.math.princeton.edu/), add to path by running the initpath script:
`matlab -nodesktop -nosplash -r “run initpath.m; run install.m”`
2. download and install [MANOPT](http://manopt.org/), add to path by running `importmanop.m`
3. download and extract [SHT](https://www.mathworks.com/matlabcentral/fileexchange/43856-real-complex-spherical-harmonic-transform--gaunt-coefficients-and-rotations?requestedDomain=www.mathworks.com), add to `initpath.m`
4. download and extract [CWF](https://github.com/PrincetonUniversity/cwf_denoise), add to `initpath.m`
#
```
ln -s polarch-Spherical-Harmonic-Transform-4f4cf39 SHT
ln -s ASPIREv0.13 aspire
```
