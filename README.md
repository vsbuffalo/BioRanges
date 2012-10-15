# BioRanges - Ranges for biological data in Python

This module is experimental. I do not recommend you use this for
anything at this stage.

## Design 

In some cases, we don't want heavyweight ranges, and in other cases we
do. Consider two cases: a few HSPs per 100,000 contigs, and 100,000
ranges on 23 chromosomes. In the first case, constructing an interval
tree per contig would be overkill. Any efficiency gained in some
operations (coverage, overlap calculation) in the former case would be
washed out by setting up the interval tree data structure. In the
second case, interval trees would be absolutely necessary.

For this reason, I think a good ranges module would have two
submodules: one for lightweight ranges and one for ranges with an
interval tree backend. Currently I have started working on
`BioRanges.lightweight`.

The attributes and methods of all classes must be generic: the
backends *must* be transparent. Minor differences in private
attrbiutes are acceptable of course, to allow for different data
structures. We may wish to implement interval trees in Cython too. 

All interfaces are modelled after Bioconductor's `GenomicRanges` and
`IRanges`. I highly recommend these over this module if your end goal
is analysis. These are much more mature packages, and R with
Bioconductor is a better environment for analysis in my opinion. But
for processing lots of data, Python can be a more comfortable
environment.

There are some key differences between BioRanges and GenomicRanges:
data is stored in a dictionary, so it won't have the same structure as
a GenomicRanges `GRanges` `elementMetaData` `DataFrame`. But this
allows us to store BioPython's HSPs and other objects more
easily. BioRanges won't have the same expressivity in terms of
interval operations.

## Development

Please help! Email me at vsbuffaloAAAA@ucdavis.edu (sans poly-A tail)
if you wish to join, or just clone and send a pull request.

## Todo

 - Unit tests
 - All interval tree backend code, efficiency testing.
