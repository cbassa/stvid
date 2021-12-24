# STVID

**stvid** is a set of *python* programs to detect and identify satellite tracks on observations of the night sky, and measure the satellite positions to determine and/or update their orbits.

## Background
This software will take over the satellite tracking functionality of C.Bassa's [sattools](https://github.com/cbassa/sattools). By porting the functionality to *python*, and using [astropy](https://github.com/astropy/astropy), and [opencv](https://opencv-python-tutroals.readthedocs.io/en/latest/), the software is indented to be easier to install and operate.

This repository is currently used for development, and relies, for the moment, on programs from the [sattools](https://github.com/cbassa/sattools) repository (namely `satid`).

## Installation
**stvid** handles requirements using pip. You can install requirements by running :

`pip install -r requirements.txt`

Consider using a VirtualEnv to run stvid on a separate python virtual environment.

### Additional requirements
* [Hough3dLines](https://gitlab.com/pierros/hough3d-code)
  * `hough3dlines` needs to be executable in your path
* [Sextractor](https://github.com/astromatic/sextractor)
  * `sextractor` needs to be executable in your path
* [Astromentry.net](https://github.com/dstndstn/astrometry.net)
  * `solve-field` needs to be executable in your path, with index files appropriate for your FOV
* [ZWO ASI SDK](https://astronomy-imaging-camera.com/software-drivers) if using a ZWO camera
* [Picamerax](https://github.com/labthings/picamerax) if using a Raspberry Pi camera. It can be installed by running:
`sudo pip install "picamerax[array]"`

## Installation on a Raspberry Pi for Acquistion.

**stvid** can be run on a Raspberry Pi 4B 8Gb with the Raspberry HQ Camera. Installation is not as easy as a on regular Linux and this installation step only supports the filming using acquire.py. The processing can be done on a dektop PC or laptop which has the complete installations of **sattools** and **stvid**. 

### Requirements

Make sure the camera is enabled and tested using raspistill. For more info see: 

[https://www.raspberrypi.com/documentation/accessories/camera.html#raspistill](https://www.raspberrypi.com/documentation/accessories/camera.html#raspistill)


Using the HQ camera under python is currently (December 2021) only supported on Buster and not on the newer Bullseye. Run:

`cat /etc/os-release`

Check that it says `buster`. If you are reading this in the future, you should check if your version of the OS supports the camera under python.

For installing we need the OS to be up to date:

```
sudo apt-get update --allow-releaseinfo-change
sudo apt-get upgrade
```

Also cmake and some build tools need to be installed:

```
sudo apt-get install cmake
sudo apt install build-essential
```

**stvid** runs in pyhton3, it is probably already installed, but if not:

`sudo apt install python3`

Next the installation of OpenCv is needed. This is based on [https://www.pyimagesearch.com/2018/09/19/pip-install-opencv/](https://www.pyimagesearch.com/2018/09/19/pip-install-opencv/)

Consider using a VirtualEnv to run stvid on a separate python virtual environment.

```
cd
sudo apt-get install libhdf5-dev libhdf5-serial-dev libhdf5-103
sudo apt-get install libqtgui4 libqtwebkit4 libqt4-test python3-pyqt5
sudo apt-get install libatlas-base-dev
sudo apt-get install libjasper-dev
wget https://bootstrap.pypa.io/get-pip.py
sudo python3 get-pip.py
sudo pip install opencv-contrib-python==4.1.0.25
```

This is a good moment to test OpenCV:

```
pi@allsky:~/stvid $ python3
Python 3.7.3 (default, Jan 22 2021, 20:04:44) 
[GCC 8.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import cv2
>>> cv2.__version__
'4.1.0'
>>> quit()
pi@allsky:~/stvid $ 
```

Install support for the Rapberry Pi cam in python:

`sudo pip install "picamerax[array]"`


The requirements are handled using pip. You can install requirements by running :

```
pip install -r requirements-pi.txt

```

### STVID

Stvid can be installed from this repository:

```
cd
git clone https://github.com/EelkeVisser/stvid.git
cd stvid
cp configuration.ini-dist configuration.ini
```

## Installation on a Raspberry Pi for Processing.
This chapter is TBD

## Configuration
* Edit `configuration.ini` with your preferred settings. Especially [Common], [Credentials] and [Camera]: camera_type = PI. You may also want to experiment with the analog_gain and digital_gain settings.
* Make sure ST_DATADIR is set as an env variable pointing to the /data folder of sattools

## Running

Run acquire.py to start filming:
```
cd
cd stvid
./acquire.py -c /home/pi/stvid/configuration.ini 
```

You may want to add `-l` for a live view window. And for testing during daylight use `-t 120`, for a test of 120 seconds. 

## Todo
Features to be implemented.

#### High priority
* ~~Use sunset/sunrise times for starting/stopping data acquisition.~~
* ~~Automatic astrometric calibration.~~
* ~~Recognize unidentified satellite/meteor tracks using [3D Hough transform](http://www.ipol.im/pub/art/2017/208/).~~

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
* Create start up script in `/etc/init.d`. Call capture script as user with `su <username> -c "acquire.py"`.

## License
&copy; 2018-2022 Cees Bassa

Licensed under the [GPLv3](LICENSE).
