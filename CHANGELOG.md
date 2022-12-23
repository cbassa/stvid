# Changelog

## dev-20221122
### Breaking Changes / Porting Guide
- Improved configuration for orbital elements
  - New configuration section `Elements`
  - Configuration field `Common.tle_path` moved to `Elements.tlepath`
- Removed old `process.py` script; Renamed `process_new.py` to `process.py` (dec640a)
- New dependency: pyyaml
- New data product: `{timestamp}_data.yaml`

## dev-20221021
- First tagged version of the 'dev'-branch
- Removed sattools dependency; (Addded much more lightweight satpredict dependency)
- Renamed configuration fields section Common to `Observer`, specifically:
  `Common.observer_lat` to `Observer.latitude`
  `Common.observer_lon` to `Observer.longitude`
  `Common.observer_height` to `Observer.height`
  `Common.cospar` to `Observer.cospar`
  `Common.name` to `Observer.name`

# Changes between `process_new.py` and `process.py`
- No `results_path` anymore

