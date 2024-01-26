import sys
import configparser


def load_config(args):
    """
    Load STVID configuration from one or more files, as configured by the command-line arguments.

    Arguments
    args (argparse.Namespace):
        The result of argparse.parse_args
    """
    cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    conf_file = args.conf_file if args.conf_file else "configuration.ini"
    result = cfg.read(conf_file)

    if not result:
        print("Could not read config file: %s\nExiting..." % conf_file)
        sys.exit(1)
    return cfg
