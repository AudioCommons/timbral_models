## AudioCommons Timbral Models
The timbral models were devleoped by the [Institute of Sound Recording (IoSR)](http://www.iosr.uk/projects/AudioCommons/) at the University of Surrey, and was completed as part of the [AudioCommons project](https://www.audiocommons.org).  

The current distribution contains python scripts for predicting eight timbral characteristics: *hardness*, *depth*, *brightness*, *roughness*, *warmth*, *sharpness*, *booming*, and *reverberation*.  

More detailed explanations of how the models function can be found in Deliverable D5.8: Release of timbral characterisation tools for semantically annotating non-musical content, available: http://www.audiocommons.org/materials/


## Installing the package
The timbral_models package can be installed using the pip command.  This will handle installation of all dependencies.  In the update to version 0.4, the dependency to essentia was removed and only pip installable packages are required.
```
pip install timbral_models
```

Please note that during testing, pip was unable to install some of the dependencies and produced an error.  In these cases, either rerun the `pip install timbral_models` command or install the offending dependency directly, e.g. `pip install numpy`. 

The package can also be installed locally and be made editable.  To do this, clone the repository, navigate to the folder, and run the pip command `pip install -e .`.  In this method, dependencies will not be installed.
  

## Dependencies
The script can also be downloaded manually from the github repository (https://github.com/AudioCommons/timbral_models).  If doing this, dependencies will need to be manually installed.  The timbral models rely on several other easily accessible python packages: `numpy`, `soundfile`, `librosa`, `sklearn`, and `scipy`.  These are all easily installed using the `pip install` command.  e.g.
```
$ pip install numpy
$ pip install soundfile
$ pip install librosa
$ pip install scipy
$ pip install sklearn
$ pip install six
$ pip install pyloudnorm   
```


## Using the models
The models are formatted in a python package that can be simply imported into a Python script.
The timbral extractor can be used to extract all timbral attributes with a single function call.

To calculate the timbral attributes, pass the timbral extractor function a string of the filename.  The method will then read in the audio file internally and return all timbral characteristics.
```
import timbral_models
fname = '/Documents/Music/TestAudio.wav'
timbre = timbral_models.timbral_extractor(fname)
```
In this example, `timbre` will be a python dictionary containing the predicted *hardness*, *depth*, *brightness*, *roughness*, *warmth*, *sharpness*, *booming*, and *reverberation* of the specified audio file.  


### Single attribute calculation

Alternative, each timbral attribute can be calculated individually by calling the specific timbral function, e.g. `timbral_hardness(fname)`.
These are named `timbral_xxx(fname)`, where `xxx` represents the timbral model, and also require a string of the filename to be analysed.
```
import timbral_models
fname = '/Documents/Music/TestAudio.wav'
timbre = timbral_models.timbral_hardness(fname)
```


## Model output
The *hardness*, *depth*, *brightness*, *roughness*, *warmth*, *sharpness*, and *booming* are regression based models, trained on subjective ratings ranging from 0 to 100.  However, the output may be beyond these ranges.   
The `clip_output` optional parameter can be  used to contrain the outputs between 0 and 100.
```
timbre = timbral_models.timbral_extractor(fname, clip_output=True)
```   
For additional optional parameters, please see Deliverable D5.8.

The *reverb* attribute is a classification model, returning 1 or 0, indicating if the file "sounds reverberant" or "does not sound reverberant" respectively.


## MATLAB Reverb model
Also contained in this repository is a full version of the timbral reverb model.  For instruction on installing and using this, please see Deliverable D5.8.

## Version History
This section documents the version history of the timbral models.  To download a specific version of the model that relate to a specific deliverable, please check this section and download the most recent version from that date.

2019/01/24 - Version 0.4 of timbral_models, relates to Audio Commons Deliverable D5.8.  This version of the repository relates to the software version 0.4 on PyPI.

2018/12/14 - Version 0.3 of timbral models, relates to Audio Commons Deliverable D5.7. This version of the repository relates to the software version 0.3 on PyPI.

2018/07/26 - Version 0.2 of timbral models, relates to Audio Commons Deliverable D5.6.  This version of the repository relates to the software version 0.2 on PyPI. 

2017/09/05 - Version 0.1 of timbral models, relates to Audio Commons Deliverable D5.3.  This version of the repository relates to the software version 0.1 on PyPI.

2017/04/27 - Version 0.0 of the timbral models, relates to Audio Commons Deliverable D5.2. 


## Citation
For refencing these models, please reference Deliverable D5.8, available: http://www.audiocommons.org/materials/