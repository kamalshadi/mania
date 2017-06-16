Usage:

Installation:

-> pip install mania

In order to start using the package import the library:

-> import mania

The package exposes three functions.

———————————————————

Running mania on a subject probtrackx results:

-> mania_on_subject(study_path, number_of_streamlines_per_seed)

———————————————————

Calculating confidence of network edges:

-> conf(study_path)


———————————————————

Generating synthetic data (see the mania paper):

-> synth_probabilistic_anatomy(Number_of_nodes,density,mu1,mu2)


Note: <study_path> must point to a subject folder. Inside the subject folder,
there should be a subfolder called "probtrackx" in which the results of probtrackx
resides - see the sample_subject folder as an example.


