from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='timbral_models',
      version='0.4.0',
      description='Algorithms for predicting the timbral characteristics of audio files',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/AudioCommons/timbral_models',
      author='Andy Pearce',
      author_email='andy.pearce@surrey.ac.uk',
      license='Apache Software',
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Topic :: Multimedia :: Sound/Audio :: Analysis',
          'License :: OSI Approved :: Apache Software License',
          'Programming Language :: Python :: 2.7',
      ],
      packages=['timbral_models'],
      install_requires=[
          'numpy',
          'soundfile',
          'scipy',
          'librosa',
          'sklearn',
          'pyloudnorm',
          'six'
      ],
      zip_safe=False)
