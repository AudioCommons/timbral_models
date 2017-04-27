# AudioCommons Timbral Models
This project contains Python scripts developed for extracting timbral attributes of audio files.

More detailed explanations of how the models function can be found in Deliverable D5.2: First prototype of timbral characterisation tools for semantically annotating non-musical content, available: http://www.audiocommons.org/materials/

# Dependencies
The timbral models rely on several other easily accessible python packages: `numpy`, `soundfile`, `librosa`, and `scipy`.  These are all easily installed using the `pip install` command.  e.g.
```
$ pip install numpy
$ pip install soundfile
$ pip install librosa
$ pip install scipy

```

# Using the models
Currently, the models are written in a format so they can be imported into a Python script.  
Each script may contain many methods, but the method which should be called is the `timbral_xxx(fname)` method.
To calculate the timbral attribute, give the method a string of the file name.  The method will then read in the audio file internally.

# Model output
The *hardness*, *depth*, and *brightness* models predict subjective ratings of their respective attributes.  Each model returns a float.  These models were trained on subjective ratings ranging from 0 to 100, but can extend beyond this range.

The *metalic-nature* model returns the probability of an audio file sounding metallic as a float ranging from 0.0 to 1.0.

The *roughness* model returns a float representing the roughness of the audio file.  The minimum roughness value is 0.0, but there is no upper limit on the maximum roughness value.

The *reverb* model returns an approximation of the RT60 of an audio file in milliseconds. If the algorithm cannot estimate an RT60, 0 is returned.  

# Example usage

```
import Timbral_Brightness as bright

# generic file location
fname = '/Users/User/Music/AudioFileToTest.wav'

# calculate brightness
brightness = bright.timbral_brightness(fname) 
```
