

def infer_reference_name(path):
    # NCBI builds and hg releases aren't identical
    # but the differences are all on chrM and unplaced contigs
    candidates = {
        'GRCh36' : ['hg18', 'b36', 'B36', 'GRCh36'],
        'GRCh37' : ['hg19', 'b37', 'B37', 'GRCh37'],
        'GRCh38' : ['b38', 'B38', 'GRCh38'],
    }

    for name in sorted(candidates.keys(), reverse=True):
        aliases = candidates[name]
        for alias in aliases:
            if alias in path:
                return name

    raise ValueError(
        "Failed to infer human genome assembly name for %s" % path)

