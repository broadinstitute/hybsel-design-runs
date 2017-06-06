"""Utilities for working with sequence i/o.
"""

from collections import OrderedDict
import re
import textwrap

import numpy as np

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def read_fasta(fn, data_type='str', replace_degenerate=False):
    """Read a FASTA file.

    Args:
        fn: path to FASTA file to read
        data_type: determines whether to store a sequence as
            a native Python string ('str') or as a numpy array
            ('np')
        replace_degenerate: when True, replace the degenerate
            bases ('Y','R','W','S','M','K') with 'N'

    Returns:
        dict mapping the name of each sequence to the sequence
        itself. The mapping is ordered by the order in which
        the sequence is encountered in the FASTA file; this
        helps in particular with replicating past results,
        where the input order could affect the output.
    """
    degenerate_pattern = re.compile('[YRWSMK]')

    m = OrderedDict()
    with open(fn) as f:
        curr_seq_name = ""
        for line in f:
            line = line.rstrip()
            if len(line) == 0:
                continue
            if curr_seq_name == "":
                # Must encounter a new sequence
                assert line.startswith('>')
            if len(line) == 0:
                # Reset the sequence being read on an empty line
                curr_seq_name = ""
            elif line.startswith('>'):
                curr_seq_name = line[1:]
                m[curr_seq_name] = ''
            else:
                # Append the sequence
                if replace_degenerate:
                    line = degenerate_pattern.sub('N', line)
                m[curr_seq_name] += line

    if data_type == 'str':
        # Already stored sequence as string
        m_converted = m
    elif data_type == 'np':
        m_converted = OrderedDict()
        for seq_name, seq in m.items():
            m_converted[seq_name] = np.fromstring(seq, dtype='S1')
    else:
        raise ValueError("Unknown data_type " + data_type)

    return m_converted


def iterate_fasta(fn, data_type='str', replace_degenerate=False):
    """Scan through a FASTA file and yield each sequence.

    This is a generator that scans through a given FASTA file and,
    upon completing the read of a sequence, yields that sequence.

    Args:
        fn: path to FASTA file to read
        data_type: determines whether to store a sequence as
            a native Python string ('str') or as a numpy array
            ('np')
        replace_degenerate: when True, replace the degenerate
            bases ('Y','R','W','S','M','K') with 'N'

    Yields:
        each sequence in the FASTA file
    """
    degenerate_pattern = re.compile('[YRWSMK]')

    def format_seq(seq):
        if data_type == 'str':
            # Already stored as str
            return seq
        elif data_type == 'np':
            return np.fromstring(seq, dtype='S1')
        else:
            raise ValueError("Unknown data_type " + data_tyoe)

    with open(fn) as f:
        curr_seq = ''
        for line in f:
            line = line.rstrip()
            if len(line) == 0:
                continue
            if line.startswith('>'):
                # Yield the current sequence (if there is one) and reset the
                # sequence being read
                if len(curr_seq) > 0:
                    yield format_seq(curr_seq)
                curr_seq = ''
            else:
                # Append the sequence
                if replace_degenerate:
                    line = degenerate_pattern.sub('N', line)
                curr_seq += line
        if len(curr_seq) > 0:
            yield format_seq(curr_seq)


def write_fasta(seqs, out_fn, chars_per_line=70):
    """Write sequences to a FASTA file.

    Args:
        seqs: dict (or OrderedDict) mapping the name of each
            sequence to the sequence itself
        out_fn: path to FASTA file to write
        chars_per_line: the number of characters to put on
            each line when writing a sequence
    """
    with open(out_fn, 'w') as f:
        for seq_name, seq in seqs.items():
            if isinstance(seq, np.ndarray):
                # Convert seq to a string
                seq = ''.join(seq)
            f.write('>' + seq_name + '\n')
            seq_wrapped = textwrap.wrap(seq, chars_per_line)
            for seq_line in seq_wrapped:
                f.write(seq_line + '\n')
            f.write('\n')
