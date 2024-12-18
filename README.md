The objective of this code is to use a deterministic mathematical model to computationally predict the appropriate loading of drugs in intravaginal rings (IVRs) for animal models that best reproduce drug pharmacokinetics (PK) in macaque and sheep.

Dependencies: MATLAB 2022 or above. 

main.m is the main file that calls or the required functions for the analysis. It's set up like a Jupyter notebook, so please run each section independently ("Run Section") in MATLAB to produce the results step by step.

main_parametric.m is a cleaned up version set up purely for running the parametric analysis in the cluster. It does no plotting for compatibility with the cluster, but plotting for these results can be done in post-processing using the various plotting functions available.

