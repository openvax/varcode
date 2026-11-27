# PYTHON_ARGCOMPLETE_OK
"""The ``varcode`` command: one entry point with subcommands.

Subcommands:

``effects``        annotate variants with predicted coding effects
``genes``          annotate variants with overlapping genes

For backward compatibility, ``varcode --vcf x.vcf ...`` with no subcommand is
the same as ``varcode effects --vcf x.vcf ...``.
"""

import os
import sys

import argcomplete

from ..version import __version__
from .help import HelpArgumentParser

SUBCOMMANDS = ("effects", "genes")


def build_parser(*, attach_subcommands=False):
    """Build the top-level ``varcode`` argument parser without heavy imports."""
    parser = HelpArgumentParser(
        prog="varcode",
        description="Variant annotation in Python.",
        epilog="Run `varcode <subcommand> --help` for the options of each subcommand.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="varcode %s" % __version__,
    )
    subparsers = parser.add_subparsers(
        title="subcommands",
        dest="subcommand",
        metavar="{effects,genes}",
    )
    effects_parser = subparsers.add_parser(
        "effects",
        help="annotate variants with overlapping genes and predicted coding effects (default when no subcommand is given)",
    )
    genes_parser = subparsers.add_parser(
        "genes",
        help="annotate variants with overlapping gene names",
    )
    if attach_subcommands:
        from .effects_script import attach_effects_parser
        from .genes_script import attach_genes_parser

        for subparser in (effects_parser, genes_parser):
            subparser.add_argument(
                "--version",
                action="version",
                version="varcode %s" % __version__,
            )
        attach_effects_parser(effects_parser)
        attach_genes_parser(genes_parser)
    return parser


def run_subcommand(name, args):
    """
    Run one ``varcode`` subcommand, importing its module only when needed.

    Parameters
    ----------
    name : str
        One of ``effects`` or ``genes``.
    args : list of str
        Arguments for the subcommand, excluding its name.

    Returns
    -------
    int or None
        Exit status of the subcommand; ``None`` means success.

    Raises
    ------
    ValueError
        If ``name`` is not a subcommand.
    """
    if name == "effects":
        from .effects_script import main as subcommand_main
    elif name == "genes":
        from .genes_script import main as subcommand_main
    else:
        raise ValueError(
            "Unknown varcode subcommand %r; expected one of %s"
            % (name, ", ".join(SUBCOMMANDS)))
    return subcommand_main(args)


def main(argv=None):
    """
    Run the ``varcode`` command.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments excluding the program name (default: sys.argv[1:]).

    Returns
    -------
    int or None
        Exit status; ``None`` means success.
    """
    if argv is None and "_ARGCOMPLETE" in os.environ:
        parser = build_parser(attach_subcommands=True)
        argcomplete.autocomplete(parser)
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] in SUBCOMMANDS:
        return run_subcommand(args[0], args[1:])
    if args and args[0] in ("-h", "--help"):
        build_parser().print_help()
        return 0
    if args and args[0] == "--version":
        print("varcode %s" % __version__)
        return 0
    # Legacy form: no subcommand means `varcode effects`.
    return run_subcommand("effects", args)


if __name__ == "__main__":
    sys.exit(main())
