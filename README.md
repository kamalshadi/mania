# MANIA
This python package implements algorithms introduced in the paper:

*Shadi K, Bakhshi S, Gutman DA, Mayberg HS, Dovrolis C. A Symmetry-Based Method to Infer Structural Brain Networks from Probabilistic Tractography Data. Frontiers in Neuroinformatics. 2016.*

For technical details, please consult to the paper above. 

## Usage 

### Installation
Install the latest release with
```
-> pip install pymania
```
There are no hard dependencies other than the Python standard library and numpy. MANIA runs with Python 2.7.

### Importing the package in python environment:

```
-> import mania
```

### API

The package exposes four functions listed below.

* Running mania for a single subject's probtrackx results:

```
-> mania_on_subject(study_path, number_of_streamlines_per_seed)
```
> **Note:** <study_path> must point to a subject folder. Inside the subject folder,
there should be a subfolder called "probtrackx" in which the results of probtrackx
resides - see the sample_subject folder as an example.

* Calculating the confidence metric for the edges of a subject network:
```
-> conf(study_path)
```
* Running mania at group level:
```
-> group_mania(subject_list,output_folder)
```
subject_list is the list of study folders. each element of the list is a study folder for a subject.

* Generating synthetic data (see the mania paper):
```
-> synth_probabilistic_anatomy(Number_of_nodes,density,mu1,mu2)
```


