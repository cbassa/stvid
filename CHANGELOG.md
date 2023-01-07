# Changelog

## dev-20230104
### Breaking Changes
- Improved calibration logic (Fixes #90):
  process.py will wait for new files to be processed, even if there has been no successful plate solving.
  (see below for new CLI option `--batch` to restore old behavior)

- Changed timestamp format in output filenames,
  use hyphens instead of colons as time seperator (Commit 11da524);
  old: `2023-01-05T00:58:00.988.fits`, new: `2023-01-05T00-58-00.988.fits`

### Other Fixes & Improvements
- Added new CLI option to process.py: `--batch`:
  "Batch process observations, exit when done" (Commit fb064b5)
- Added new CLI option to process.py: `--reprocess` (Commit fd1bb88):
  "Remove processed files and start from scratch".
  This will remove files matching the following patterns:
  - `test.fits`
  - `*.png`,
  - `*.cat`
  - `*.cal`
  - `*.csv`
  - `*.dat`
- Added row/column masks to process.py (Commit 3267c3f)
  New configuration fields:
  - `LineDetection.rows_to_mask`
  - `LineDetection.columns_to_mask`
- Added support for splitting configuration into multiple config files (Commit 585b43e)
- Moved from YAML to JSON format for the process data otput, renamed data product
  from `{timestamp}_data.yaml` to `{timestamp}_data.json` (Commit 4c49c48)
- Added more fields to the process data output (Fixes #92):
  - start
  - exptime
  - ra, dec
  - sx, sy
  - wx, wy
- Speed up argmax by reconfiguring z array to allow looping over axis=2 (Closes #88)
- Fixed bug in processing, overhauled identification strategy (Commit 00fe0f4)
- Added correct COSPAR for unknowns which was hard-coded to `` before (Commit 4f4eca8)
- Improved process.py: Write FITS to temp first, then rename (Commit 6f8545a)
- Added multithreaded processing (Commit 4254025); New configuration field:
  `LineDetection.cpu_count`

## dev-20221122
### Breaking Changes / Porting Guide
- Improved configuration for orbital elements
  - New configuration section `Elements`
  - Configuration field `Common.tle_path` moved to `Elements.tlepath`
- Removed old `process.py` script; Renamed `process_new.py` to `process.py` (Commit dec640a)
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

