# qdTCSPC
A simple program to help analyze TCSPC data from picoquant Microtime200 in python, specifically with QD analysis in mind.

The big requirements to get this to run that you may not have are:
1. pqreader from phconvert found here: https://phconvert.readthedocs.io/en/latest/
  This is a great library I found that reads picoquant's .ptu file format (T3 only!)
2. pycorrelate to do some FCS/Antibunching, found here: https://pypi.org/project/pycorrelate/
3. Numpy, matplotlib, scipy, pretty standard scientific python things as far as I can tell.


The idea to use this is to create a "ptu" object, that will hold all your photon stream data. Once you load your TTTR stream
in with pqreader you can then do lifetime plotting, blinking, antibunching, and some other analysis. The resulting data is 
usually saved to the "ptu" object so you can plot it again later. Still a work in progress! Need to write some good documentation!

Usage:
1. Make sure all dependencies are fulfilled
2. Run the qdtcspc.py script in your console/notebook, etc.
