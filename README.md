# STVID

**stvid** is a set of *python* programs to detect and identify satellite tracks on video observations of the night sky, and measure the satellite positions to determine and/or update their orbits.

This software will take over the satellite tracking functionality of my [sattools](https://github.com/cbassa/sattools). By porting the functionality to *python*, and using [astropy](https://github.com/astropy/astropy), and [opencv](https://opencv-python-tutroals.readthedocs.io/en/latest/), the software is easier to install and operate.

This repository is currently used for development, and relies, for the moment, on programs from the [sattools](https://github.com/cbassa/sattools) repository.

## Installation
**stvid** handles requirements using pip. You can install requirements by running `pip install -r requirements.txt`. You should consider using a VirtualEnv to run stvid on a separate python virtual environment. You will also need to download the [Tycho 2 catalog](https://www.dropbox.com/s/isbq1wsmfr24zx5/tyc2.fits?dl=0).

## Configuration
* Copy the `configuration.ini-dist` file to `configuration.ini`
* Edit `configuration.ini` with your preferred settings

## Running
TBD

## Todo
Features to be implemented.

#### High priority
* ~~Use sunset/sunrise times for starting/stopping data acquisition.~~
* ~~Automatic astrometric calibration.~~
* Recognize unidentified satellite/meteor tracks using [3D Hough transform](http://www.ipol.im/pub/art/2017/208/).

#### Medium priority
* Pause data acquisition of the current line-of-sight (alt/az) is in the Earth's shadow for a particular orbital altitude.
* Investigate sensitivity loss of `significance=(max-mean)/sigma` if the four frame images are stored as 8bit integers instead of floats.


#### Low priority
* Implement python based star finding (stick with *source extractor* for now).
* Migrate to [python based SGP4/SDP4 algorithms](https://github.com/brandon-rhodes/python-sgp4)
* Use masks to mask unilluminated CCD areas.
* Investigate automatic submission of IOD measurements to [SeeSat-L](http://www.satobs.org/seesat/).
* ~~Migrate user settings to a configuration file.~~

## Run acquisition at startup

* Add user to video group (`sudo adduser <username> video`).
* Add video device to udev rules (add `SUBSYSTEM=="video1", GROUP="video", MODE="0660"` in `/etc/udev/rules.d/10-webcam.rules`).
* Create start up script in `/etc/init.d`. Call capture script as user with `su <username> -c "capture_1.sh"`.

## License
&copy; 2018-2019 Cees Bassa

Licensed under the [GPLv3](LICENSE).
