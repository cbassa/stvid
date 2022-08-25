#!/usr/bin/env python3
import sys
import os
import numpy as np
import cv2
import time
import ctypes
import multiprocessing
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.io import fits
import astropy.units as u
from stvid.utils import get_sunset_and_sunrise
import logging
import configparser
import argparse
import zwoasi as asi
from picamerax.array import PiRGBArray
from picamerax import PiCamera


# Capture images from pi
def capture_pi(image_queue, z1, t1, z2, t2, nx, ny, nz, tend, device_id, live, cfg):
    # Intialization
    first = True
    slow_CPU = False

    # Initialize cv2 device
    camera = PiCamera(sensor_mode=2)
    camera.resolution = (nx, ny)    
    # Turn off any thing automatic.
    camera.exposure_mode = 'off'        
    camera.awb_mode = 'off'
    # ISO needs to be 0 otherwise analog and digital gain won't work.
    camera.iso = 0
    # set the camea settings
    camera.framerate = cfg.getfloat(camera_type, 'framerate')
    camera.awb_gains = (cfg.getfloat(camera_type, 'awb_gain_red'), cfg.getfloat(camera_type, 'awb_gain_blue'))    
    camera.analog_gain = cfg.getfloat(camera_type, 'analog_gain')
    camera.digital_gain = cfg.getfloat(camera_type, 'digital_gain')
    camera.shutter_speed = cfg.getint(camera_type, 'exposure')

    rawCapture = PiRGBArray(camera, size=(nx, ny))
    # allow the camera to warmup
    time.sleep(0.1)

    try:
        # Loop until reaching end time
        while float(time.time()) < tend:
            # Wait for available capture buffer to become available
            if (image_queue.qsize() > 1):
                logger.warning("Acquiring data faster than your CPU can process")
                slow_CPU = True
            while (image_queue.qsize() > 1):
                time.sleep(0.1)
            if slow_CPU:
                lost_video = time.time() - t
                logger.info("Waited %.3fs for available capture buffer" % lost_video)
                slow_CPU = False

            # Get frames
            i = 0
            for frameA in camera.capture_continuous(rawCapture, format="bgr", use_video_port=True):
                            
                # Store start time
                t0 = float(time.time())                
                # grab the raw NumPy array representing the image, then initialize the timestamp                
                frame = frameA.array
                                    
                # Compute mid time
                t = (float(time.time())+t0)/2.0
                
                # Skip lost frames
                if frame is not None:
                    # Convert image to grayscale
                    z = np.asarray(cv2.cvtColor(
                        frame, cv2.COLOR_BGR2GRAY)).astype(np.uint8)
                    # optionally rotate the frame by 2 * 90 degrees.    
                    # z = np.rot90(z, 2)
                
                    # Display Frame
                    if live is True:                            
                        cv2.imshow("Capture", z)    
                        cv2.waitKey(1)
                    
                    # Store results
                    if first:
                        z1[i] = z
                        t1[i] = t
                    else:
                        z2[i] = z
                        t2[i] = t
                        
                # clear the stream in preparation for the next frame
                rawCapture.truncate(0)
                # count up to nz frames, then break out of the for loop.
                i += 1
                if i >= nz:
                    break
                
            if first: 
                buf = 1
            else:
                buf = 2
            image_queue.put(buf)
            logger.debug("Captured z%d" % buf)

            # Swap flag
            first = not first
        reason = "Session complete"
    except KeyboardInterrupt:
        print()
        reason = "Keyboard interrupt"
    except ValueError as e:
        logger.error("%s" % e)
        reason = "Wrong image dimensions? Fix nx, ny in config."
    finally:
        # End capture
        logger.info("Capture: %s - Exiting" % reason)
        camera.close()



# Capture images from cv2
def capture_cv2(image_queue, z1, t1, z2, t2, nx, ny, nz, tend, device_id, live):
    # Intialization
    first = True
    slow_CPU = False

    # Initialize cv2 device
    device = cv2.VideoCapture(device_id)

    # Set properties
    device.set(3, nx)
    device.set(4, ny)

    try:
        # Loop until reaching end time
        while float(time.time()) < tend:
            # Wait for available capture buffer to become available
            if (image_queue.qsize() > 1):
                logger.warning("Acquiring data faster than your CPU can process")
                slow_CPU = True
            while (image_queue.qsize() > 1):
                time.sleep(0.1)
            if slow_CPU:
                lost_video = time.time() - t
                logger.info("Waited %.3fs for available capture buffer" % lost_video)
                slow_CPU = False

            # Get frames
            for i in range(nz):
                # Store start time
                t0 = float(time.time())

                # Get frame
                res, frame = device.read()

                # Compute mid time
                t = (float(time.time())+t0)/2.0

                # Skip lost frames
                if res is True:
                    # Convert image to grayscale
                    z = np.asarray(cv2.cvtColor(
                        frame, cv2.COLOR_BGR2GRAY)).astype(np.uint8)

                    # Display Frame
                    if live is True:
                        cv2.imshow("Capture", z)
                        cv2.waitKey(1)

                    # Store results
                    if first:
                        z1[i] = z
                        t1[i] = t
                    else:
                        z2[i] = z
                        t2[i] = t

            if first: 
                buf = 1
            else:
                buf = 2
            image_queue.put(buf)
            logger.debug("Captured z%d" % buf)

            # Swap flag
            first = not first
        reason = "Session complete"
    except KeyboardInterrupt:
        print()
        reason = "Keyboard interrupt"
    except ValueError as e:
        logger.error("%s" % e)
        reason = "Wrong image dimensions? Fix nx, ny in config."
    finally:
        # End capture
        logger.info("Capture: %s - Exiting" % reason)
        device.release()


# Capture images
def capture_asi(image_queue, z1, t1, z2, t2, nx, ny, nz, tend, device_id, live, cfg):
    first    = True  # Array flag
    slow_CPU = False # Performance issue flag

    
    camera_type  = "ASI"
    gain         = cfg.getint(camera_type, 'gain')
    maxgain      = cfg.getint(camera_type, 'maxgain')
    autogain     = cfg.getboolean(camera_type, 'autogain')
    exposure     = cfg.getint(camera_type, 'exposure')
    binning      = cfg.getint(camera_type, 'bin')
    brightness   = cfg.getint(camera_type, 'brightness')
    bandwidth    = cfg.getint(camera_type, 'bandwidth')
    high_speed   = cfg.getint(camera_type, 'high_speed')
    hardware_bin = cfg.getint(camera_type, 'hardware_bin')
    sdk          = cfg.get(camera_type, 'sdk')
    try:
        software_bin = cfg.getint(camera_type, 'software_bin')
    except configparser.Error:
        software_bin = 0

    # Initialize device
    asi.init(sdk)

    num_cameras = asi.get_num_cameras()
    if num_cameras == 0:
        logger.error("No ZWOASI cameras found")
        raise ValueError
        sys.exit()

    cameras_found = asi.list_cameras()  # Models names of the connected cameras

    if num_cameras == 1:
        device_id = 0
        logger.info("Found one camera: %s" % cameras_found[0])
    else:
        logger.info("Found %d ZWOASI cameras" % num_cameras)
        for n in range(num_cameras):
            logger.info("    %d: %s" % (n, cameras_found[n]))
        logger.info("Using #%d: %s" % (device_id, cameras_found[device_id]))

    camera = asi.Camera(device_id)
    camera_info = camera.get_camera_property()
    logger.debug("ASI Camera info:")
    for (key, value) in camera_info.items():
        logger.debug("  %s : %s" % (key,value))

    camera.set_control_value(asi.ASI_BANDWIDTHOVERLOAD, bandwidth)
    camera.disable_dark_subtract()
    camera.set_control_value(asi.ASI_GAIN, gain, auto=autogain)
    camera.set_control_value(asi.ASI_EXPOSURE, exposure, auto=False)
    camera.set_control_value(asi.ASI_AUTO_MAX_GAIN, maxgain)
    camera.set_control_value(asi.ASI_AUTO_MAX_BRIGHTNESS, 20)
    camera.set_control_value(asi.ASI_WB_B, 99)
    camera.set_control_value(asi.ASI_WB_R, 75)
    camera.set_control_value(asi.ASI_GAMMA, 50)
    camera.set_control_value(asi.ASI_BRIGHTNESS, brightness)
    camera.set_control_value(asi.ASI_FLIP, 0)
    try:
        camera.set_control_value(asi.ASI_HIGH_SPEED_MODE, high_speed)
    except:
        pass
    try:
        camera.set_control_value(asi.ASI_HARDWARE_BIN, hardware_bin)
    except:
        pass
    camera.set_roi(bins=binning)
    camera.start_video_capture()
    camera.set_image_type(asi.ASI_IMG_RAW8)

    try:
        # Fix autogain
        if autogain:
            while True:
                # Get frame
                z = camera.capture_video_frame()

                # Break on no change in gain
                settings = camera.get_control_values()
                if gain == settings["Gain"]:
                    break
                gain = settings["Gain"]
                camera.set_control_value(asi.ASI_GAIN, gain, auto=autogain)

        # Loop until reaching end time
        while float(time.time()) < tend:
            # Wait for available capture buffer to become available
            if (image_queue.qsize() > 1):
                logger.warning("Acquiring data faster than your CPU can process")
                slow_CPU = True
            while (image_queue.qsize() > 1):
                time.sleep(0.1)
            if slow_CPU:
                lost_video = time.time() - t
                logger.info("Waited %.3fs for available capture buffer" % lost_video)
                slow_CPU = False

            # Get settings
            settings = camera.get_control_values()
            gain = settings["Gain"]
            temp = settings["Temperature"]/10.0
            logger.info("Capturing frame with gain %d, temperature %.1f" % (gain, temp))

            # Set gain
            if autogain:
                camera.set_control_value(asi.ASI_GAIN, gain, auto=autogain)

            # Get frames
            for i in range(nz):
                # Store start time
                t0 = float(time.time())

                # Get frame
                z = camera.capture_video_frame()

                # Apply software binning
                if software_bin > 1:
                    my, mx = z.shape
                    z = cv2.resize(z, (mx // software_bin, my // software_bin))
                
                # Compute mid time
                t = (float(time.time())+t0)/2.0

                # Display Frame
                if live is True:
                    cv2.imshow("Capture", z)
                    cv2.waitKey(1)

                # Store results
                if first:
                    z1[i] = z
                    t1[i] = t
                else:
                    z2[i] = z
                    t2[i] = t

            if first: 
                buf = 1
            else:
                buf = 2
            image_queue.put(buf)
            logger.debug("Captured buffer %d (%dx%dx%d)" % (buf, nx, ny, nz))

            # Swap flag
            first = not first
        reason = "Session complete"
    except KeyboardInterrupt:
        print()
        reason = "Keyboard interrupt"
    except ValueError as e:
        logger.error("%s" % e)
        reason = "Wrong image dimensions? Fix nx, ny in config."
    except MemoryError as e:
        logger.error("Capture: Memory error %s" % e)
    finally:
        # End capture
        logger.info("Capture: %s - Exiting" % reason)
        camera.stop_video_capture()
        camera.close()


def compress(image_queue, z1, t1, z2, t2, nx, ny, nz, tend, path, device_id, cfg):
    """ compress: Aggregate nframes of observations into a single FITS file, with statistics.

        ImageHDU[0]: mean pixel value nframes         (zmax)
        ImageHDU[1]: standard deviation of nframes    (zstd)
        ImageHDU[2]: maximum pixel value of nframes   (zmax)
        ImageHDU[3]: maximum pixel value frame number (znum)

    Also updates a [observations_path]/control/state.txt for interfacing with satttools/runsched and sattools/slewto
    """
    # Force a restart
    controlpath = os.path.join(path, "control")
    if not os.path.exists(controlpath):
        try:
            os.makedirs(controlpath)
        except PermissionError:
            logger.error("Can not create control path directory: %s" % controlpath)
            raise
    with open(os.path.join(controlpath, "state.txt"), "w") as fp:
        fp.write("restart\n")

    try:
        # Start processing
        while True:
            # Check mount state
            restart = False
            with open(os.path.join(controlpath, "state.txt"), "r") as fp:
                line = fp.readline().rstrip()
                if line == "restart":
                    restart = True

            # Restart
            if restart:
                # Log state
                with open(os.path.join(controlpath, "state.txt"), "w") as fp:
                    fp.write("observing\n")

                # Get obsid
                trestart = time.gmtime()
                obsid = "%s_%d/%s" % (time.strftime("%Y%m%d", trestart), device_id, time.strftime("%H%M%S", trestart))
                filepath = os.path.join(path, obsid)
                logger.info("Storing files in %s" % filepath)

                # Create output directory
                if not os.path.exists(filepath):
                    try:
                        os.makedirs(filepath)
                    except PermissionError:
                        logger.error("Can not create output directory: %s" % filepath)
                        raise

                # Get mount position
                with open(os.path.join(controlpath, "position.txt"), "r") as fp:
                    line = fp.readline()
                with open(os.path.join(filepath, "position.txt"), "w") as fp:
                    fp.write(line)

            # Wait for completed capture buffer to become available
            while (image_queue.qsize == 0):
                time.sleep(0.1)
                
            # Get next buffer # from the work queue
            proc_buffer = image_queue.get()
            logger.debug("Processing buffer %d" % proc_buffer)

            # Log start time
            tstart = time.time()

            # Process first buffer
            if proc_buffer == 1:
                t = t1                
                z = z1
            elif proc_buffer == 2:
                t = t2
                z = z2

            # Format time
            nfd = "%s.%03d" % (time.strftime("%Y-%m-%dT%T",
                            time.gmtime(t[0])), int((t[0] - np.floor(t[0])) * 1000))
            t0 = Time(nfd, format='isot')
            dt = t - t[0]

            # Cast to 32 bit float
            z = z.astype("float32")
            
            # Compute statistics
            zmax = np.max(z, axis=0)
            znum = np.argmax(z, axis=0)
            zs1 = np.sum(z, axis=0) - zmax
            zs2 = np.sum(z * z, axis=0) - zmax * zmax 
            zavg = zs1 / float(nz - 1)
            zstd = np.sqrt((zs2 - zs1 * zavg) / float(nz - 2))

            # Convert to float and flip
            zmax = np.flipud(zmax.astype("float32"))
            znum = np.flipud(znum.astype("float32"))
            zavg = np.flipud(zavg.astype("float32"))
            zstd = np.flipud(zstd.astype("float32"))

            # Generate fits
            fname = "%s.fits" % nfd

            # Format header
            hdr = fits.Header()
            hdr['DATE-OBS'] = "%s" % nfd
            hdr['MJD-OBS']  = t0.mjd
            hdr['EXPTIME']  = dt[-1]-dt[0]
            hdr['NFRAMES']  = nz
            hdr['CRPIX1']   = float(nx)/2.0
            hdr['CRPIX2']   = float(ny)/2.0
            hdr['CRVAL1']   = 0.0
            hdr['CRVAL2']   = 0.0
            hdr['CD1_1']    = 1.0/3600.0
            hdr['CD1_2']    = 0.0
            hdr['CD2_1']    = 0.0
            hdr['CD2_2']    = 1.0/3600.0
            hdr['CTYPE1']   = "RA---TAN"
            hdr['CTYPE2']   = "DEC--TAN"
            hdr['CUNIT1']   = "deg"
            hdr['CUNIT2']   = "deg"
            hdr['CRRES1']   = 0.0
            hdr['CRRES2']   = 0.0
            hdr['EQUINOX']  = 2000.0
            hdr['RADECSYS'] = "ICRS"
            hdr['COSPAR']   = cfg.getint('Common', 'observer_cospar')
            hdr['OBSERVER'] = cfg.get('Common', 'observer_name')
            hdr['SITELONG'] = cfg.getfloat('Common', 'observer_lon')
            hdr['SITELAT'] = cfg.getfloat('Common', 'observer_lat')
            hdr['ELEVATIO'] = cfg.getfloat('Common', 'observer_height')
            if cfg.getboolean('Astrometry', 'tracking_mount'):
                hdr['TRACKED'] = 1
            else:
                hdr['TRACKED'] = 0
            for i in range(nz):
                hdr['DT%04d' % i] = dt[i]
            for i in range(10):
                hdr['DUMY%03d' % i] = 0.0

            # Write fits file
            hdu = fits.PrimaryHDU(data=np.array([zavg, zstd, zmax, znum]),
                                header=hdr)
            hdu.writeto(os.path.join(filepath, fname))
            logger.info("Compressed %s in %.2f sec" % (fname, time.time() - tstart))

            # Exit on end of capture
            if t[-1] > tend:
                break
            logger.debug("Processed buffer %d" % proc_buffer)
            

    except KeyboardInterrupt:
        pass
    except MemoryError as e:
        logger.error("Compress: Memory error %s" % e)
    finally:
        # Exiting
        logger.info("Exiting compress")


def main():
    # Read commandline options
    conf_parser = argparse.ArgumentParser(description='Capture and compress' +
                                                      ' live video frames.')
    conf_parser.add_argument('-c', '--conf_file',
                             help="Specify configuration file. If no file" +
                             " is specified 'configuration.ini' is used.",
                             metavar="FILE")
    conf_parser.add_argument('-t', '--test', 
                             nargs='?',
                             action='store', 
                             default=False,
                             help='Testing mode - Start capturing immediately for (optional) seconds',
                             metavar="s")
    conf_parser.add_argument('-l', '--live', action='store_true',
                             help='Display live image while capturing')

    args = conf_parser.parse_args()

    # Process commandline options and parse configuration
    cfg = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
    
    conf_file = args.conf_file if args.conf_file else "configuration.ini"
    result = cfg.read([conf_file])

    if not result:
        print("Could not read config file: %s\nExiting..." % conf_file)
        sys.exit()

    # Setup logging
    logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] " +
                                     "[%(levelname)-5.5s]  %(message)s")
    logger = logging.getLogger()

    # Generate directory
    path = os.path.abspath(cfg.get('Common', 'observations_path'))
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except PermissionError:
            logger.error("Can not create observations_path: %s" % path)
            sys.exit()

    fileHandler = logging.FileHandler(os.path.join(path, "acquire.log"))
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler(sys.stdout)
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)
    logger.setLevel(logging.DEBUG)

    logger.info("Using config: %s" % conf_file)

    # Testing mode
    if args.test is None:
        test_duration = 31
        testing = True
    elif args.test is not False:
        test_duration = int(args.test)
        testing = True
    else:
        testing = False
    logger.info("Test mode: %s" % testing)
    if (testing):
        logger.info("Test duration: %ds" % test_duration)

    # Live mode
    live = True if args.live else False
    logger.info("Live mode: %s" % live)

    # Get camera type
    camera_type = cfg.get('Camera', 'camera_type')

    # Get device id
    device_id = cfg.getint(camera_type, 'device_id')

    # Current time
    tnow = Time.now()

    # Set location
    loc = EarthLocation(lat=cfg.getfloat('Common', 'observer_lat')*u.deg,
                        lon=cfg.getfloat('Common', 'observer_lon')*u.deg,
                        height=cfg.getfloat('Common', 'observer_height')*u.m)

    if not testing:
        # Reference altitudes
        refalt_set  = cfg.getfloat('Control', 'alt_sunset')*u.deg
        refalt_rise = cfg.getfloat('Control', 'alt_sunrise')*u.deg

        # FIXME: The following will fail without internet access
        #        due to failure to download finals2000A.all
        # Get sunrise and sunset times
        state, tset, trise = get_sunset_and_sunrise(tnow, loc, refalt_set, refalt_rise)

        # Start/end logic
        if state == "sun never rises":
            logger.info("The sun never rises. Exiting program.")
            sys.exit()
        elif state == "sun never sets":
            logger.info("The sun never sets.")
            tend = tnow+24*u.h
        elif (trise < tset):
            logger.info("The sun is below the horizon.")
            tend = trise
        elif (trise >= tset):
            dt = np.floor((tset-tnow).to(u.s).value)
            logger.info("The sun is above the horizon. Sunset at %s."
                        % tset.isot)
            logger.info("Waiting %.0f seconds." % dt)
            tend = trise
            try:
                time.sleep(dt)
            except KeyboardInterrupt:
                sys.exit()
    else:
        tend = tnow + test_duration*u.s

    logger.info("Starting data acquisition")
    logger.info("Acquisition will end after "+tend.isot)

    # Get settings
    nx = cfg.getint(camera_type, 'nx')
    ny = cfg.getint(camera_type, 'ny')
    nz = cfg.getint(camera_type, 'nframes')

    # Initialize arrays
    z1base = multiprocessing.Array(ctypes.c_uint8, nx*ny*nz)
    z1 = np.ctypeslib.as_array(z1base.get_obj()).reshape(nz, ny, nx)
    t1base = multiprocessing.Array(ctypes.c_double, nz)
    t1 = np.ctypeslib.as_array(t1base.get_obj())
    z2base = multiprocessing.Array(ctypes.c_uint8, nx*ny*nz)
    z2 = np.ctypeslib.as_array(z2base.get_obj()).reshape(nz, ny, nx)
    t2base = multiprocessing.Array(ctypes.c_double, nz)
    t2 = np.ctypeslib.as_array(t2base.get_obj())

    image_queue = multiprocessing.Queue()

    # Set processes
    pcompress = multiprocessing.Process(target=compress,
                                        args=(image_queue, z1, t1, z2, t2, nx, ny,
                                              nz, tend.unix, path, device_id, cfg))
    if camera_type == "PI":
        pcapture = multiprocessing.Process(target=capture_pi,
                                           args=(image_queue, z1, t1, z2, t2,
                                                 nx, ny, nz, tend.unix, device_id, live, cfg))
    elif camera_type == "CV2":
        pcapture = multiprocessing.Process(target=capture_cv2,
                                           args=(image_queue, z1, t1, z2, t2,
                                                 nx, ny, nz, tend.unix, device_id, live))
    elif camera_type == "ASI":
        pcapture = multiprocessing.Process(target=capture_asi,
                                           args=(image_queue, z1, t1, z2, t2,
                                                 nx, ny, nz, tend.unix, device_id, live, cfg))

    # Start
    pcapture.start()
    pcompress.start()

    # End
    try:
        pcapture.join()
        pcompress.join()
    except (KeyboardInterrupt, ValueError):
        time.sleep(0.1) # Allow a little time for a graceful exit
    except MemoryError as e:
        logger.error("Memory error %s" % e)
    finally:
        pcapture.terminate()
        pcompress.terminate()

    # Release device
    if live is True:
        cv2.destroyAllWindows()


if __name__ == '__main__':
    main()
