"""Module for shutter control."""
import gpiozero
import logging

LOGGER = logging.getLogger(__name__)

class Shutter:
    """
    Class to control the satnogs-eye shutter motor.

    Can be used as context manager in a with-statement,
    to ensure the shutter is always closed when acquistion
    is stopped.
    """
    def __init__(self, pin):
        self.motor = gpiozero.LED(pin=pin)
    def open(self):
        self.motor.on()
        LOGGER.info("Shutter opened.")
    def close(self):
        self.motor.off()
        LOGGER.info("Shutter closed.")
