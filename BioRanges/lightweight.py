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
from collections import Counter

def verify_arg_length(msg, args):
    """
    Check whether the lists of arguments supplied are the same (if
    they are not still None). If not, raise ValueError with message
    msg. Return the length of arguments.
    """
    arg_lens = set([len(x) for x in args if x is not None])
    if len(arg_lens) > 1:
        raise ValueError(msg)
    if len(arg_lens) == 0:
        return 0
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
        if (start > end or (width is not None and width < 0) or
            start < 0 or end < 0):
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
        self.name = name

    def __repr__(self):
        if self.name is not None:
            return "Range '%s' over [%d, %d]" % (self.name, self.start, self.end)
        return "Range over [%d, %d]" % (self.start, self.end)

    def overlaps(self, other):
        """
        Return a boolean indicating whether two ranges overlap.
        """

        if not isinstance(other, "Range"):
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
            widths_i = widths[i] if widths is not None else None
            names_i = names[i] if names is not None else None
            self._ranges.append(Range(starts[i], ends[i], widths_i, names_i))

    def __repr__(self):
        """
        Representation of Ranges collection using a few sample rows.
        """
        lines = ["Ranges with %d ranges" % len(self)]
        header = ["start", "end", "width", "name"]
        rows = [header]
        ncols = range(len(header))
        max_col_width = [len(c) for c in header]
        for i, rng in enumerate(self._ranges):
            this_row = map(str, [rng.start, rng.end, rng.width, rng.name])
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

    def append(self, other):
        """
        Add on a new Range object to this Ranges object.
        """
        if isinstance(other, Range):
            self._ranges.append(other)
        elif isinstance(other, Ranges):
            self._ranges.extend(other._ranges)
        elif isinstance(other, list):
            class_ok = [isinstance(x, Range) for x in other]
            if not all(class_ok):
                raise ValueError("append() method can only handle lists "
                                 "where each element is a Range")
            for rng in other:
                self._ranges.append(rng)
        else:
            raise ValueError("append() method can only objects of class "
                             "list, Range, and Ranges")


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

    def __init__(self, range, seqname, strand, data=dict()):
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
        repr_str = "SeqRange on '%s', strand '%s' at [%d, %d], %d data keys"
        return repr_str % (self.seqname, self.strand,
                           self.range.start, self.range.end, len(self.data))

    def overlaps(self, other):
        """
        Return a boolean indicating whether two ranges overlap. Since
        these are SeqRanges, we have to consider strand and
        seqname. Following GRanges, we will require the are the same;
        to test overlaps ignoring strand, either a different method
        will be added, or strands should be changed to "*".
        """
        if self.seqname != other.seqname or self.strand != other.strand:
            return False
        return self.range.overlaps(other)

    @property
    def start(self):
        """
        Return start position (accessor for range.start).
        """
        return self.range.start

    @property
    def end(self):
        """
        Return end position (accessor for range.end).
        """
        return self.range.end

    @property
    def width(self):
        """
        Return width (accessor for range.width)
        """
        return self.range.width

    def sliceseq(self, seq):
        """
        Given a sequence, return the sequence in the region.
        """
        return seq[self.range.start:self.range.end]

    def maskseq(self, seq, mask_char="X"):
        """
        Mask a sequence based on SeqRange, uses standard BioRanges
        0-based indexing.
        """
        start = self.range.start
        end = self.range.end
        mask_len = end - start + 1
        masked = seq[:start] + mask_char*mask_len + seq[end+1:]
        return masked

class SeqRanges(object):
    """
    A container class for a set of ranges on a sequence (chromosome,
    contig, etc).
    """

    def __init__(self, ranges=None, seqnames=None, strands=None, data_list=None,
                 seqlengths=dict()):
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
        args = [ranges, seqnames, strands, data_list]
        not_none_args = [a for a in args if a is not None]
        arg_len = verify_arg_length("list of ranges, seqnames, strands, and "
                                     "data_list must be of the same length", not_none_args)

        self._ranges = list()
        for i in range(arg_len):
            rng = ranges[i]
            if data_list is not None:
                self._ranges.append(SeqRange(rng, seqnames[i], strands[i], data=data_list[i]))
            else:
                self._ranges.append(SeqRange(rng, seqnames[i], strands[i]))

        self.seqlengths = seqlengths        

    def __repr__(self):
        return self.show()
        
    def show(self, keys=list()):
        """
        Representation of SeqRanges collection using a few sample
        rows, with possible specified data keys also displayed
        """
        if not isinstance(keys, list) or not isinstance(keys, tuple):
            raise ValueError("'keys' argument must be list or tuple")
        lines = ["SeqRanges with %d ranges" % len(self)]
        header = ["seqnames", "ranges", "strand"]
        ncols_seqranges = len(header)
        header.extend(keys)
        rows = [header]
        ncols = range(len(header))
        max_col_width = [len(c) for c in header]
        for i, seqrange in enumerate(self._ranges):
            rng = seqrange.range
            this_row = [seqrange.seqname,
                        "[%d, %d]" % (rng.start, rng.end),
                        str(seqrange.strand)]
            for key in keys:
                this_row.append(str(seqrange.data.get(key, "")))
            max_col_width = [max((len(this_row[j]), max_col_width[j])) for j in ncols]
            rows.append(this_row)

        # now, add appropriate formating and spacing
        for row in rows:
            tmp_line = ""
            for i, col in enumerate(row):
                if i > 0:
                    tmp_line += " "
                if i == ncols_seqranges:
                    tmp_line += "| "
                tmp_line += " "*(max_col_width[i] - len(col)) + col
            lines.append(tmp_line)

        print "\n".join(lines)

    def append(self, other):
        """
        Add on a new SeqRange object to this SeqRanges object.
        """
        if isinstance(other, SeqRange):
            self._ranges.append(other)
        elif isinstance(other, SeqRanges):
            self._ranges.extend(other._ranges)

        elif isinstance(other, list):
            class_ok = [isinstance(x, SeqRange) for x in other]
            if not all(class_ok):
                raise ValueError("append() method can only handle lists "
                                 "where each element is a SeqRange")
            for seqrange in other:
                self._ranges.append(seqrange)
        else:
            raise ValueError("append() method can only objects of class "
                             "list, SeqRange, and SeqRanges")

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
        if not isinstance(seqrange, SeqRange):
            raise ValueError("assignment can only handle SeqRange objects")
        self._ranges[i] = seqrange

    def __getitem__(self, i):
        """
        Get a SeqRange from a SeqRanges collection.
        """
        return self._ranges[i]

    def __delitem__(self, i):
        """
        Delete a SeqRange from a SeqRanges collection.
        """
        del(self._ranges[i])

    @property
    def start(self):
        """
        Get list of all start positions.
        """
        return [r.range.start for r in self._ranges]

    @property
    def end(self):
        """
        Get list of all end positions.
        """
        return [r.range.end for r in self._ranges]

    @property
    def width(self):
        """
        Get list of all widths of ranges.
        """
        return [r.range.width for r in self._ranges]

    @property
    def strand(self):
        """
        Get list of all strand.
        """
        return [r.range.strand for r in self._ranges]

    def overlaps(self):
        raise ValueError("lightweight Ranges objects do not the "
                         "support overlap() method")
