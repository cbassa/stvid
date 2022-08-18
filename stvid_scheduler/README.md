# STVID Scheduler

The STVID scheduler is a Python daemon which has the following tasks:
* Read configuration from the client-ansible ("SatNOGS Cloud Config")
* Start/Stop the the aquisition using STVID `acquire.py`
* Run calibration (plate-solver), satellite detection and identification using STVID `process.py`
* Regularly download TLEs from
  - McCants "classfd" TLEs
  - McCants "integrated" TLEs
  - Space-Track
* Upload observation results to SatNOGS

## Configuration

The STVID scheduler is expecting a configuration as defined by specs/development.json (there will be releases of this specification in the future).
It will then use this configuration to write the configuration file expected by the various commands of STVID.

## Installation

The STVID scheduler can be run by loading the provided SystemdD service unit file.
