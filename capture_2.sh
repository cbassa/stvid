#!/bin/bash

# Sleep
sleep 10

# Does /dev/video2 exist
ls -l /dev/video2 >/tmp/camera_2.log 2>&1

# Run mplayer for 10 frames to setup easycap
mplayer -frames 10 tv:// -tv device=/dev/video2 >/tmp/camera_2.log 2>&1

# Start script
python /home/bassa/software/stvid/acquire.py 2 >/tmp/camera_2.log 2>&1 &
