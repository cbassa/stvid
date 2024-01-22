import sys
import configparser


def load_config(conf_files):
    """
    Load STVID configuration from one or more files.

    Arguments
    conf_files (list[str]): List of configuration files
    """
    cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    result = cfg.read(conf_files)

    if not result:
        print(f"Could not read config files: {conf_files}\nExiting...")
        sys.exit(1)
    return cfg


def add_argument_conf_file(parser):
    parser.add_argument("-c", "--conf_file",
                            help="Specify configuration file(s). If no file" +
                            " is specified 'configuration.ini' is used.",
                            action="append",
                            nargs="?",
                            metavar="FILE",
                            dest="conf_files",
                            default=["configuration.ini"])
    return parser
