#!/bin/bash

# Sleep
sleep 10

# Does /dev/video1 exist
ls -l /dev/video1 >/tmp/camera_1.log 2>&1

# Run mplayer for 10 frames to setup easycap
mplayer -frames 10 tv:// -tv device=/dev/video1 >/tmp/camera_1.log 2>&1

# Start script
python /home/bassa/software/stvid/acquire.py 1 >/tmp/camera_1.log 2>&1 &
