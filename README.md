# AudioCommons Timbral Models
This project contains Python scripts developed for extracting timbral attributes of audio files.

More detailed explanations of how the models function can be found in Deliverable D5.2: First prototype of timbral characterisation tools for semantically annotating non-musical content, available: http://www.audiocommons.org/materials/


# Installing the package
The timbral_models package can be installed using the pip command.  This will handle all dependencies.
```
pip install timbral_models
```
Note that during test installing the package with only basic python installed, an error occurred when installing dependencies.  This can be overcome by first installing numpy, followed by timbral_models.
```
pip install numpy
pip install timbral_models
```


# Dependencies
If the code is to be installed using a method other than pip, dependencies will need to be installed.  The timbral models rely on several other easily accessible python packages: `numpy`, `soundfile`, `librosa`, and `scipy`.  These are all easily installed using the `pip install` command.  e.g.
=======
# Dependencies
The timbral models rely on several other easily accessible python packages: `numpy`, `soundfile`, `librosa`, and `scipy`.  These are all easily installed using the `pip install` command.  e.g.
```
$ pip install numpy
$ pip install soundfile
$ pip install librosa
$ pip install scipy

```

# Using the models
The models are written as a package which can be imported into a Python script.  Within the package are six methods that predict the brightness, depth, hardness, metallic-nature, reverb, and roughness of an audio file.  These are named `timbral_xxx(fname)`, where `xxx` represents the timbral model.

To calculate the timbral attribute, give the method a string of the file name.  The method will then read in the audio file internally and return the timbral characteristic, as described below.

# Model output
The *hardness*, *depth*, and *brightness* models predict subjective ratings of their respective attributes.  Each model returns a float.  These models were trained on subjective ratings ranging from 0 to 100, but can extend beyond this range.

The *metalic-nature* model returns the probability of an audio file sounding metallic as a float ranging from 0.0 to 1.0.

The *roughness* model returns a float representing the roughness of the audio file.  The minimum roughness value is 0.0, but there is no upper limit on the maximum roughness value.

The *reverb* model returns an approximation of the RT60 of an audio file in milliseconds. If the algorithm cannot estimate an RT60, 0 is returned.  

# Model output
The *hardness*, *depth*, and *brightness* models predict subjective ratings of their respective attributes.  Each model returns a float.  These models were trained on subjective ratings ranging from 0 to 100, but can extend beyond this range.

The *metalic-nature* model returns the probability of an audio file sounding metallic as a float ranging from 0.0 to 1.0.

The *roughness* model returns a float representing the roughness of the audio file.  The minimum roughness value is 0.0, but there is no upper limit on the maximum roughness value.

The *reverb* model returns an approximation of the RT60 of an audio file in milliseconds. If the algorithm cannot estimate an RT60, 0 is returned.  

# Example usage

```
from timbral_models import timbral_brightness 

# generic file location
fname = '/Users/User/Music/AudioFileToTest.wav'

# calculate brightness
brightness = timbral_brightness(fname) 
```


