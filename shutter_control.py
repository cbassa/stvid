"""Module for CLI methods."""

import argparse
import readchar

from stvid.shutter import Shutter
from stvid.config import add_argument_conf_file, load_config


def run_shutter_tui(cfg):
    """Entry point for CLI."""

    print("Eye Control")
    print("Press 'o' to open, 'c' to close, or 'q' to quit.")

    shutter = Shutter(cfg.getint("Shutter", "pin"))

    while True:
        print("Waiting for command (o/c/q)...", end='', flush=True)
        key = readchar.readchar().lower()
        print(key)  # Print the pressed key

        if key == "o":
            print("Opening...")
            shutter.open()
        elif key == "c":
            print("Closing...")
            shutter.close()
        elif key == "q":
            print("Exiting... Shutter will be closed.")
            break
        else:
            print("Invalid input. Please use 'o', 'c', or 'q'.")

    # Ensure the motor is off when exiting
    shutter.close()


if __name__ == "__main__":
    conf_parser = argparse.ArgumentParser(description="Control the satnogs-eye shutter.")
    conf_parser = add_argument_conf_file(conf_parser)

    args = conf_parser.parse_args()
    cfg = load_config(args.conf_files)

    if not cfg.has_section("Shutter"):
        print("ERROR: Configuration section 'Shutter' not found.")
        return

    run_shutter_tui(cfg)
