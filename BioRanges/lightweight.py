"""
Lightweight Ranges for Biological Data.

These range classes are not meant for more than a few ranges per
sequence, nor for large overlap calculations. For that, use
BioRanges.Ranges. These will work fine for BLAST hits.

Much of the credit for interfaces goes to Bioconductor's GenomicRanges
and IRanges. These are *better* for analysis than these Python
implementation, which is designed more for processing scripts than
analysis.
"""

# Design Notes:
#
# We could also have a GenericRange and a GenericRangeCollections
# classes to prevent some code duplication. When interval trees are
# implemented as the back processing end for these classes' methods,
# we will likely go this approach.

STRAND_OPTIONS = ("+", "-", "*")
NUM_RANGES_DISPLAY = 10

import pdb
from collections import Counter, OrderedDict

def verify_arg_length(msg, args):
    """
    Check whether the lists of arguments supplied are the same (if
    they are not still None). If not, raise ValueError with message
    msg. Return the length of arguments.
    """
    arg_lens = set([len(x) for x in args if x is not None])
    if len(arg_lens) > 1:
        raise ValueError(msg)
    return list(arg_lens)[0]

class Range(object):
    """
    A basic range/interval class.
    """

    def __init__(self, start=None, end=None, width=None, name=None):
        """
        Constructor methods for creating a new range.
        """
        # check consistency of arguments
        if Counter((start, end, width))[None] > 2:
            raise ValueError("too few arguments for Range(): "
                             "need two of [start, end, width]")
        if start > end or (width is not None and width < 0):
            raise ValueError("negative range widths not allowed "
                             "(end > start and width >= 0)")

        # infer missing values (as Bioconductor's IRanges does)
        if start is None:
            start = end - width
        if end is None:
            end = start + width
        if width is None:
            width = end - start

        self.start = start
        self.end = end
        self.width = width

    def __repr__(self):
        return "Range"

    def overlaps(self, other):
        """
        Return a boolean indicating whether two ranges overlap.
        """

        if other.__class__.__name__ != "Range":
            raise ValueError("overlaps() method requires another Range object")

        return other.start <= self.end and self.start <= other.end

class Ranges(object):
    """
    Container class for Range objects.
    """

    def __init__(self, starts=None, ends=None, widths=None, names=None):
        """
        Create multiple Range objects.
        """
        # check whether the lists of arguments supplied are the same
        # (if they are not still None).
        args = [starts, ends, widths, names]
        arg_len = verify_arg_length("list of starts, ends, widths, and "
                                     "names must be of the same length", args)

        self._ranges = list()
        for i in range(arg_len):
            self._ranges.append(Range(starts[i], ends[i], widths[i], names[i]))
        
    def __len__(self):
        """
        Return number of ranges in this collection.
        """
        return len(self._ranges)

    def __delitem__(self, i):
        """
        Remove item from Ranges collection.
        """
        del(self._ranges[i])
        
    def __setitem__(self, i, range):
        """
        Set item in Ranges collection.
        """
        self._ranges[i] = range
    
    def __getitem__(self, i):
        """
        Get item from Ranges collection.
        """
        return self._ranges[i]

    @property
    def start(self):
        """
        Get list of all start positions.
        """
        return [r.start for r in self._ranges]

    @property
    def end(self):
        """
        Get list of all end positions.
        """
        return [r.end for r in self._ranges]

    @property
    def width(self):
        """
        Get list of all widths of ranges.
        """
        return [r.width for r in self._ranges]

    def overlaps(self):
        """
        Placeholder for overlaps, telling users to use non-lightweight
        version.
        """
        raise ValueError("lightweight Ranges objects do not the "
                         "support overlap() method")


class SeqRange(object):
    """
    A range on a sequence (chromosome, contig, etc).
    """

    def __init__(self, range, seqname, strand, data=None):
        """
        Constructor method for SequenceRange objects.
        """
        self.range = range
        self.seqname = seqname

        if strand not in STRAND_OPTIONS:
            raise ValueError("strand must be either: %s" % ', '.join(STRAND_OPTIONS))
        self.strand = strand
        self.data = data

    def __repr__(self):
        return "SeqRange(%d, %d, %s)" % (range.start, range.end, strand)

    def overlaps(self, other):
        """
        Return a boolean indicating whterh two ranges overlap. Since
        these are SeqRanges, we have to consider strand and
        seqname. Following GRanges, we will require the are the same;
        to test overlaps ignoring strand, either a different method
        will be added, or strands should be changed to "*".
        """
        if self.seqname != other.seqname or self.strand != other.strand:
            return False
        return self.range.overlaps(other)

class SeqRanges(object):
    """
    A container class for a set of ranges on a sequence (chromosome,
    contig, etc).
    """

    def __init__(self, ranges, strands, seqnames, datas=None):
        """
        Constructor method for SeqRange objects.
        """

        # Data structure notes:
        #
        # We use a dictionary with seqnames as the key, with the
        # strand being a defaultdict containing a list of
        # SeqRanges. We could have another layer of keys corresponding
        # to strand, but both seqnames and strand are already stored
        # in SeqRange objects, so the redundancy is just to achieve
        # O(1) lookup time. The non-lightweight implementation will do
        # this with interval trees and handle these issues throught
        # that.
        
        args = [a for a in [ranges, strands, seqnames, datas] if a is not None]
        arg_len = verify_arg_length("list of ranges, strands, seqnames, and "
                                     "datas must be of the same length", args)

        self._ranges = list()
        for i in range(arg_len):
            rng = ranges[i]
            if datas is not None:
                self._ranges.append(SeqRange(rng, strands[i], seqnames[i], datas[i]))
            else:
                self._ranges.append(SeqRange(rng, strands[i], seqnames[i]))

    def __repr__(self):
        """
        Representation of SeqRanges collection using a few sample
        rows.
        """
        lines = ["SeqRanges with %d ranges" % len(self)]
        header = ["seqnames", "ranges", "strand"]
        rows = [header]
        ncols = range(len(header))
        max_col_width = [len(c) for c in header]
        for i, seqrange in enumerate(self._ranges[:NUM_RANGES_DISPLAY]):
            rng = seqrange.range
            this_row = [seqrange.seqname,
                        "[%d, %d]" % (rng.start, rng.end),
                        str(seqrange.strand)]
            max_col_width = [max((len(this_row[j]), max_col_width[j])) for j in ncols]
            rows.append(this_row)

        # now, add appropriate formating and spacing
        for row in rows:
            tmp_line = ""
            for i, col in enumerate(row):
                if i > 0:
                    tmp_line += " "
                tmp_line += " "*(max_col_width[i] - len(col)) + col
            lines.append(tmp_line)

        return "\n".join(lines)

    def __len__(self):
        """
        Return the number of ranges in this object.
        """
        return len(self._ranges)

    def __setitem__(self, i, seqrange):
        """
        Set item in SeqRanges collection; these are done by index
        only.
        """
        if seqrange.__class__.__name__ != "SeqRange":
            raise ValueError("assignment can only handle SeqRange objects")
        self._ranges[i] = seqrange

    def __getitem__(self, i):
        """
        Get a SeqRange from a SeqRanges collection.
        """
        return self._ranges[i]

    def overlaps(self):
        raise ValueError("lightweight Ranges objects do not the "
                         "support overlap() method")
