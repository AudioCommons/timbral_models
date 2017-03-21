# AudioCommons Timbral Models
This project contains Python scripts developed for extracting timbral attributes of audio files.

# Using the models
Currently, the models are written in a format so they can be imported into a Python script.  
Each script may contain many methods, but the method which should be called is the `timbral_xxx(fname)` method.
To calculate the timbral attribute, give the method a string of the file name.  The method will then read in the audio file internally.

# Example usage

```
import Timbral_Brightness as bright

# generic file location
fname = '/Users/User/Music/AudioFileToTest.wav'

# calculate brightness
brightness = bright.timbral_brightness(fname) 
```
