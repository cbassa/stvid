#!/usr/bin/env python3
import argparse
import readchar
import configparser

from stvid.shutter import Shutter

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


def main():
    conf_parser = argparse.ArgumentParser(description="Control the satnogs-eye shutter.")
    conf_parser.add_argument("-c", "--conf_file",
                             help="Specify configuration file(s). If no file" +
                             " is specified 'configuration.ini' is used.",
                             action="append",
                             nargs="?",
                             metavar="FILE")
    args = conf_parser.parse_args()
    cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))

    if not cfg.has_section("Shutter"):
        print("ERROR: Configuration section 'Shutter' not found.")
        return

    run_shutter_tui(cfg)

    
if __name__ == "__main__":
    main()
