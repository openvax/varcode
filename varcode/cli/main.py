"""Dispatch optional subcommands while preserving the annotation CLI."""

import sys


def main(args_list=None):
    args = list(sys.argv[1:] if args_list is None else args_list)
    if args and args[0] == "check-samples":
        from .sample_checks import main as check_main
        return check_main(args[1:])
    from .effects_script import main as effects_main, arg_parser
    arg_parser.epilog = "Sample integrity: varcode check-samples --help"
    return effects_main(args)


if __name__ == "__main__":
    sys.exit(main())
