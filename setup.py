from setuptools import setup

setup(name='timbral_models',
      version='0.1',
      description='Algorithms for predicting the timbral characteristics of audio files',
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
          'librosa'
      ],
      zip_safe=False)
