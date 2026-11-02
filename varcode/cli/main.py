"""The ``varcode`` command: one entry point with subcommands.

Subcommands:

``effects``        annotate variants with predicted coding effects
``genes``          annotate variants with overlapping genes
``check-samples``  screen VCFs for sample mix-ups

For backward compatibility, ``varcode --vcf x.vcf ...`` with no subcommand is
the same as ``varcode effects --vcf x.vcf ...``, and the separate
``varcode-genes`` script is an alias for ``varcode genes``.
"""

import sys

from .. import __version__

USAGE = """\
usage: varcode [-h] [--version] {effects,genes,check-samples} ...

Variant annotation in Python.

subcommands:
  effects         annotate variants with overlapping genes and predicted coding
                  effects (default when no subcommand is given)
  genes           annotate variants with overlapping gene names
                  (alias: the varcode-genes command)
  check-samples   screen VCFs for sample mix-ups

Run `varcode <subcommand> --help` for the options of each subcommand.
"""


def run_subcommand(name, args):
    """
    Run one ``varcode`` subcommand, importing its module only when needed.

    Parameters
    ----------
    name : str
        One of ``effects``, ``genes`` or ``check-samples``.
    args : list of str
        Arguments for the subcommand, excluding its name.

    Returns
    -------
    int or None
        Exit status of the subcommand; ``None`` means success.

    Raises
    ------
    KeyError
        If ``name`` is not a subcommand.
    """
    if name == "effects":
        from .effects_script import main as subcommand_main
    elif name == "genes":
        from .genes_script import main as subcommand_main
    elif name == "check-samples":
        from .sample_checks import main as subcommand_main
    else:
        raise KeyError(name)
    return subcommand_main(args)


SUBCOMMANDS = ("effects", "genes", "check-samples")


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
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] in SUBCOMMANDS:
        return run_subcommand(args[0], args[1:])
    if args and args[0] in ("-h", "--help"):
        sys.stdout.write(USAGE)
        return 0
    if args and args[0] == "--version":
        print("varcode %s" % __version__)
        return 0
    # Legacy form: no subcommand means `varcode effects`.
    return run_subcommand("effects", args)


if __name__ == "__main__":
    sys.exit(main())
