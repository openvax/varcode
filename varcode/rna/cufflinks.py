# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import

import logging

import pandas as pd
import numpy as np

# default column names from cufflinks tracking files
# for gene and isoform expression levels
STATUS_COLUMN = "FPKM_status"
ID_COLUMN = "tracking_id"
FPKM_COLUMN = "FPKM"
LOCUS_COLUMN = "locus"
GENE_NAMES_COLUMN = "gene_short_name"

def load_cufflinks_dataframe(
        filename,
        id_column=ID_COLUMN,
        fpkm_column=FPKM_COLUMN,
        status_column=STATUS_COLUMN,
        locus_column=LOCUS_COLUMN,
        gene_names_column=GENE_NAMES_COLUMN,
        drop_failed=True,
        drop_lowdata=False,
        drop_hidata=True,
        replace_hidata_fpkm_value=None,
        drop_nonchromosomal_loci=True,
        drop_novel=False):
    """
    Loads a Cufflinks tracking file, which contains expression levels
    (in FPKM: Fragments Per Kilobase of transcript per Million fragments)
    for transcript isoforms or whole genes. These transcripts/genes may be
    previously known (in which case they have an Ensembl ID) or a novel
    assembly from the RNA-Seq data (in which case their IDs look like "CUFF.1")

    Parameters
    ----------

    filename : str
        Filename of tracking file e.g. "genes.tracking_fpkm"

    id_column : str, optional

    fpkm_column : str, optional

    status_column : str, optional
        Name of column which indicates the FPKM estimate status. The column
        name is typically "FPKM_status". Possible contained within this column
        will be OK, FAIL, LOWDATA, HIDATA.

    locus_column : str, optional

    gene_names_column : str, optional

    drop_failed : bool, optional
        Drop rows whose FPKM status is "FAIL" (default=True)

    drop_lowdata : bool, optional
        Drop rows whose FPKM status is "LOWDATA", meaning that Cufflinks thought
        there were too few reads to accurately estimate the FPKM (default=False)

    drop_hidata : bool, optional
        Drop rows whose FPKM status is "HIDATA", meaning that too many
        fragments aligned to a feature for Cufflinks to process. Dropping
        the most expressed genes seems like a stupid idea so: default=False

    replace_hidata_fpkm_value : float, optional
        If drop_hidata=False, the HIDATA entries will still have an FPKM=0.0,
        this argument lets you replace the FPKM with some known constant.

    drop_nonchromosomal_loci : bool, optional
        Drop rows whose location isn't on a canonical chromosome
        i.e. doesn't start with "chr" (default=True)

    drop_novel : bool, optional
        Drop genes or isoforms that aren't found in Ensembl (default = False)

    Returns DataFrame with columns:
        id : str
        novel : bool
        fpkm : float
        chr : str
        start : int
        end : int
        gene_names : str list
    """
    df = pd.read_csv(filename, sep="\s+")

    for flag, status_value in [
            (drop_failed, "FAIL"),
            (drop_lowdata, "LOWDATA"),
            (drop_hidata, "HIDATA")]:
        mask = df[status_column] == status_value
        mask_count = mask.sum()
        total_count = len(df)
        if flag and mask_count > 0:
            verb_str = "Dropping"
            df = df[~mask]
        else:
            verb_str = "Keeping"
        logging.info(
            "%s %d/%d entries from %s with status=%s",
            verb_str,
            mask_count,
            total_count,
            filename,
            status_value)

    if drop_nonchromosomal_loci:
        loci = df[locus_column]
        chromosomal_loci = loci.str.startswith("chr")
        n_dropped = (~chromosomal_loci).sum()
        if n_dropped > 0:
            logging.info("Dropping %d/%d non-chromosomal loci from %s" % (
                n_dropped, len(df), filename))
            df = df[chromosomal_loci]

    if replace_hidata_fpkm_value:
        hidata_mask = df[status_column] == "HIDATA"
        n_hidata = hidata_mask.sum()
        logging.info(
            "Setting FPKM=%s for %d/%d entries with status=HIDATA",
            replace_hidata_fpkm_value,
            n_hidata,
            len(df))
        df[fpkm_column][hidata_mask] = replace_hidata_fpkm_value

    if len(df) == 0:
        raise ValueError("Empty FPKM tracking file: %s" % filename)

    ids = df[id_column]
    known = ids.str.startswith("ENS")

    if known.sum() == 0:
        raise ValueError("No Ensembl IDs found in %s" % filename)

    if drop_novel:
        n_dropped = (~known).sum()
        if n_dropped > 0:
            logging.info("Dropping %d/%d novel entries from %s",
                n_dropped, len(df), filename)
            df = df[known]
            known = np.ones(len(df), dtype='bool')

    loci = df[locus_column]

    # capture all characters after 'chr' but before ':'
    chromosomes = loci.str.extract("chr([^:]*):.*")
    # capture all characters after e.g. 'chr1:', which look like '132-394'
    ranges = loci.str.extract("chr[^:]*:(.*)")
    # capture all numbers before the dash
    starts = ranges.str.extract("(\d*)-\d*").astype(int)
    # capture all numbers after the dash
    ends = ranges.str.extract("\d*-(\d*)")

    # gene names are given either as "-" or a comma separated list
    # e.g. "BRAF1,PFAM2"
    gene_names_strings = df[gene_names_column]
    gene_names_strings[gene_names_strings == "-"] = ""
    # split each entry into a list of zero or more strings
    gene_names_lists = gene_names_strings.str.split(",")

    return pd.DataFrame({
        "id": df[id_column],
        "novel": ~known,
        "fpkm": df[fpkm_column],
        "chr": chromosomes,
        "start": starts,
        "end": ends,
        "gene_names": gene_names_lists
    })

def load_cufflinks_dict(*args, **kwargs):
    """
    Returns dictionary mapping feature identifier (either transcript or gene ID)
    to a DataFrame row with fields:
        id : str
        novel : bool
        fpkm : float
        chr : str
        start : int
        end : int
        gene_names : str list
    """
    return {
        row.id: row
        for (_, row)
        in load_cufflinks_dataframe(*args, **kwargs).iterrows()
    }

def load_cufflinks_fpkm_dict(*args, **kwargs):
    """
    Returns dictionary mapping feature identifier (either transcript or gene ID)
    to FPKM expression value.
    """
    return {
        row.id: row.fpkm
        for (_, row)
        in load_cufflinks_dataframe(*args, **kwargs).iterrows()
    }
