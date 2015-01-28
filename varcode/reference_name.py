

def infer_reference_name(path):
    # NCBI builds and hg releases aren't identical
    # but the differences are all on chrM and unplaced contigs
    candidates = {
        'NCBI36' : ['hg18', 'b36', 'B36', 'GRCh36', 'NCBI36'],
        'GRCh37' : ['hg19', 'b37', 'B37', 'GRCh37', 'NCBI37'],
        'GRCh38' : ['b38', 'B38', 'GRCh38', 'NCBI38'],
    }

    for name in sorted(candidates.keys(), reverse=True):
        aliases = candidates[name]
        for alias in aliases:
            if alias in path:
                return name

    raise ValueError(
        "Failed to infer human genome assembly name for %s" % path)

def ensembl_release_number_for_reference_name(name):
    if name == "NCBI36":
        return 54
    elif name == "GRCh37":
        return 75
    else:
        assert name == "GRCh38", "Unrecognized reference name: %s" % name
        return 77