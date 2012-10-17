# BioRanges - Ranges for biological data in Python

BioRanges is a library for storing Range data in Python, especially
designed for ranges on genomic sequences. In the future, it will have
two modules: lightweight ranges for storing a few ranges on lots of
sequences (i.e. when BLASTing hundreds of thousands of contigs), and a
module with an interval tree backend for processing up to hundreds of
thousands of ranges on a few sequences (i.e. when storing genes). Both
would have idential interfaces so they could be swapped easily via one
`import` line, *a la*:

    from BioRanges.lightweight import SeqRanges
    from BioRanges.intervaltree import SeqRanges

## Maturity

This is an immature module, so use with caution. Currently the
lightweight module is being developed as a prototype for accessor
methods and operations. 

If you need something more mature, use the libraries below. I am only
writing this because I can customize the object interface and I can
design it with BLAST HSP ranges in mind (via lightweight module),
rather than huge sets of genes, exons, etc.

 - Bioconductor's [GenomicRanges](http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
 - Brent Pedersen's [Quicksect](https://github.com/brentp/quicksect)
 - James Taylor's [bx-python](https://bitbucket.org/james_taylor/bx-python/wiki/Home)

## Example

Like Bioconductor's excellent `IRanges` and `GenomicRanges` packages,
BioRanges has abstractions for a generic range and a genomic
(sequence) range. Respectively, these are `Range` and
`SeqRange`. Unlike Bioconductor, there are separate classes for a
single range and a collection of ranges (this makes sense, since R is
a vector-based langauge and Python is not).

`Range` has intuitive functionality:

    >>> from BioRanges.lightweight import Range, Ranges, SeqRange, SeqRanges
    >>> a = Range(100, 140)
    >>> b = Range(104, 105)
    >>> a in b
    True
    >>> c = Range(200, 200) # e.g. a SNP
    >>> c in b
    False
    >>> 200 in c
    True
    >>> a.overlaps(c)
    False
    >>> a.overlaps(b)
    True

`Ranges` collections behave like lists:

    >>> x = Ranges()
    >>> x.append(a)
    >>> x.append(b)
    >>> x.append(c)
    >>> x
    Ranges with 3 ranges
    start end width name
      100 140    40 None
      104 105     1 None
      200 200     0 None

`SeqRange` requires strand and sequence name, and optionally sequence
lengths and data.

    >>> sa = SeqRange(a, "chr1", "+", data={"gene_name":"fake-gene-1a"})
    >>> sa
    SeqRange on 'chr1', strand '+' at [100, 140], 1 data keys

Unlike GenomicRanges, data is shapeless and potentially ragged (for
better or worse). It can be accessed or set from the `SeqRange` object
like an element from a dictionary.

    >>> sa['gene_name'] = "other-fake-gene-1b"

Collections of `SeqRange` objects behave like lists, and be created
via lists of elements:

    >>> genes = SeqRanges(ranges, ["chr1"]*4, ["+"]*4, data_list=[{"gene_id":x} for x in range(4)])
    >>> genes
    SeqRanges with 4 ranges
    seqnames   ranges strand
       chr1 [0, 100]      +
       chr1 [1, 101]      +
       chr1 [2, 102]      +
       chr1 [3, 103]      +

Because `SeqRange` is potentially ragged, it isn't printed in the
`repr` method. However, some keys can optionally be provided:

    >>> print genes.show(["gene_id"])
    SeqRanges with 4 ranges
    seqnames   ranges strand | gene_id
        chr1 [0, 100]      + |       0
        chr1 [1, 101]      + |       1
        chr1 [2, 102]      + |       2
        chr1 [3, 103]      + |       3

Ordered lists of this data can be accessed by attribute (if it's
related to a SeqRange attribute), or by the `getdata()` method if it's
`SeqRange` data:

    >>> genes.width
    [100, 100, 100, 100]
    >>> genes.start
    [0, 1, 2, 3]
    >>> genes.getdata("gene_id")
    [0, 1, 2, 3]

## Development

Please help! Email me at vsbuffaloAAAA@ucdavis.edu (sans poly-A tail)
if you wish to join, or just clone and send a pull request.

## Todo

 - When interval trees are implemented, perhaps make GenericRange,
   GenericSeqRange, GenericRanges, GenericSeqRanges, etc to allow both
   lightweight and the interval tree to inherit from.

 - Unit tests 

 - All interval tree backend code, efficiency testing.
