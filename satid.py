#!/usr/bin/env python
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv



line1 = ('1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753')
line2 = ('2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667')

satellite = twoline2rv(line1, line2, wgs72)
position, velocity = satellite.propagate(2000, 6, 29, 12, 50, 19.733571)

print(position,velocity)
