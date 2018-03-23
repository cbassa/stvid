# STVID

**stvid** will be a set of *python* programs to detect and identify satellite tracks on video observations of the night sky, and measure the satellite positions to determine and/or update their orbits.

This software will take over the satellite tracking functionality of my [sattools](https://github.com/cbassa/sattools). By porting the functionality to *python*, and using [astropy](https://github.com/astropy/astropy), and [opencv](https://opencv-python-tutroals.readthedocs.io/en/latest/), the software will be easier to install and operate.

This repository will primarily be used for development, and will rely, for the moment, on programs from the [sattools](https://github.com/cbassa/sattools) repository.

## Todo
Features to be implemented.

#### High priority
* ~~se sunset/sunrise times for starting/stopping data acquisition.~~
* Manual and automatic astrometric calibration.
* Recognize unidentified satellite/meteor tracks using [3D Hough transform](http://www.ipol.im/pub/art/2017/208/).

#### Medium priority
* Pause data acquisition of the current line-of-sight (alt/az) is in the Earth's shadow for a particular orbital altitude.
* Investigate sensitivity loss of `significance=(max-mean)/sigma` if the four frame images are stored as 8bit integers instead of floats.


#### Low priority
* Implement python based star finding (stick with *source extractor* for now).
* Migrate to [python based SGP4/SDP4 algorithms](https://github.com/brandon-rhodes/python-sgp4)
* Use masks to mask unilluminated CCD areas.
* Investigate automatic submission of IOD measurements to [SeeSat-L](http://www.satobs.org/seesat/).
* Migrate user settings to a configuration file.
