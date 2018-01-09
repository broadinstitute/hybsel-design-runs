"""Download FASTAs of genomes for each dataset.
"""

import argparse
from collections import defaultdict
import hashlib
import os
import re
import shutil
import textwrap
import time

from Bio import Entrez

__author__ = 'Hayden Metsky <hayden@mit.edu>'


Entrez.email = "hayden@mit.edu"

ADDITIONAL_HUMAN_HOST = [
    "Hepatitis A virus",
    "Middle East respiratory syndrome coronavirus",
    "Human herpesvirus 6A",
    "Human herpesvirus 6B",
    "Eastern equine encephalitis virus",
    "Venezuelan equine encephalitis virus",
    "Western equine encephalitis virus",
    ## Sindbis virus ##
    "Sindbis virus",
    "Babanki virus",
    "Sindbis-like virus",
    "Kyzylagach virus",
    "Ockelbo virus",
    "Sindbis-like virus YN87448",
    ####
    "Semliki Forest virus",
    "Semliki forest virus",
    ## Human mastadenovirus E ##
    "Human adenovirus 4",
    "Human adenovirus E",
    "Simian adenovirus 22",
    "Simian adenovirus 23",
    "Simian adenovirus 24",
    "Simian adenovirus 25",
    "Simian adenovirus 25.2",
    "Simian adenovirus 26",
    "Simian adenovirus 30",
    "Simian adenovirus 36",
    "Simian adenovirus 37.1",
    "Simian adenovirus 37.2",
    "Simian adenovirus 38",
    "Simian adenovirus 39",
    ####
    "Human adenovirus 41",
    "Human adenovirus F",
    "Human adenovirus 52",
    "Simian adenovirus 1",
    ## Parechoviruses ##
    "Human parechovirus",
    "Human parechovirus 1",
    "Human parechovirus 2",
    "Human parechovirus 3",
    "Human parechovirus 4",
    "Human parechovirus 5",
    "Human parechovirus 6",
    "Human parechovirus 7",
    "Human parechovirus 8",
    "Human parechovirus type 1 PicoBank/HPeV1/a",
    "Ljungan virus",
    "Ljungan virus strain 145SL",
    "Ljungan virus strain 174F",
    ####
    "Monkeypox virus",
    "Cowpox virus",
    "Saffold virus",
    "Human coronavirus OC43",
    "Human enteric coronavirus 4408",
    "Borna disease virus 1",
    "Borna disease virus 2"
]

DATASET_PYTHON_TEMPLATE_UNSEGMENTED = "dataset_unsegmented.template.py"
DATASET_PYTHON_TEMPLATE_SEGMENTED = "dataset_segmented.template.py"

class Dataset:

    def __init__(self, name, description, is_dna,
                 is_not_believed_to_be_human_pathogen,
                 is_erv, compare_to_explicit_tax_name,
                 subset_of):
        self.name = name
        self.description = description
        self.is_dna = is_dna
        self.is_not_believed_to_be_human_pathogen = \
            is_not_believed_to_be_human_pathogen
        self.is_erv = is_erv
        self.compare_to_explicit_tax_name = compare_to_explicit_tax_name
        self.subset_of = subset_of
        self.description_pattern = re.compile('^' + description + '$')

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    @staticmethod
    def from_line(line):
        ls = line.split('\t')
        name = ls[0]
        description = ls[1]
        if len(ls) > 2:
            options = ls[2].split(',')
        else:
            options = []
        is_dna = 'dna' in options
        is_not_believed_to_be_human_pathogen = \
            'not_believed_to_be_human_pathogen' in options
        is_erv = 'erv' in options
        compare_to_explicit_tax_name = \
            'compare_to_explicit_tax_name' in options

        subset_of = None
        for option in options:
            if option.startswith('subset:'):
                if subset_of is not None:
                    raise ValueError("More than one subset was given")
                subset_of = option[len('subset:'):]

        available_options = ['dna', 'not_believed_to_be_human_pathogen',
                             'erv', 'compare_to_explicit_tax_name']
        for option in options:
            if (not option.startswith('subset:') and
                    option not in available_options):
                raise ValueError("Unknown option %s" % option)

        return Dataset(name, description, is_dna,
                       is_not_believed_to_be_human_pathogen,
                       is_erv, compare_to_explicit_tax_name,
                       subset_of)


def read_dataset_list(fn):
    datasets = []
    with open(fn) as f:
        for line in f:
            datasets += [Dataset.from_line(line.rstrip())]
    return datasets


class SequenceFromAccessionList:

    def __init__(self, representative, name, host,
                 lineage, taxonomy_name, segment):
        self.representative = representative
        self.name = name
        self.host = host
        self.lineage = lineage
        self.taxonomy_name = taxonomy_name
        self.segment = segment
        self.human_is_host = 'human' in host
        self.is_segmented = segment != 'segment'

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    @staticmethod
    def from_line(line):
        ls = line.split('\t')

        segment = ls[5]
        removals = ['RNA', 'DNA']
        for removal in removals:
            if (segment.startswith('segment ' + removal + ' ') and
                    len(segment) > len('segment ' + removal + ' ')):
                # change 'segment removal X' to 'segment X'
                segment = 'segment ' + segment[len('segment ' + removal + ' '):]

        return SequenceFromAccessionList(ls[0], ls[1], ls[2], ls[3], ls[4],
                                         segment)


def read_genome_accession_list(fn):
    sequences = []
    with open(fn) as f:
        for line in f:
            if line.startswith('#'):
                continue
            sequences += [SequenceFromAccessionList.from_line(line.rstrip())]
    return sequences


def uniqueify_genome_accession_list(sequences):
    # some sequences have the same name (neighbor) but different representatives
    # (i.e., different reference sequences) and therefore appear more than once
    # in the list; the sequences of these genomes are all we care about, so
    # "unique-ify" them so that each sequence name appears just once
    # (arbitrarily select the sequence, so in effect arbitrarily pick the
    # representative)
    return list(set(sequences))


def filter_sequences_with_nonhuman_host(sequences):
    # return only those sequences where humans are a host
    # sometimes (e.g., for 'Hepatitis A') the accession list incorrectly
    # omits human as a host, so also check against the manual
    # ADDITIONAL_HUMAN_HOST list
    return [s for s in sequences if (s.human_is_host or
                                     s.taxonomy_name in ADDITIONAL_HUMAN_HOST)]


def verify_dataset_list(datasets):
    # each dataset name should appear only once
    assert len(set([d.name for d in datasets])) == len(datasets)


def verify_sequence_names_are_unique(sequences):
    assert len(set([s.name for s in sequences])) == len(sequences)


def find_datasets_for_name(datasets, lineage_most_specific, tax_name):
    # search through dataset descriptions, which are based on taxonomy
    # names, to find the one(s) matching lineage_most_specific
    # (unless dataset is explicitly supposed to look at tax_name)
    # search linearly and find all in case there's more than one (we'll
    # want to check later on that there is exactly one)
    matches = []
    for dataset in datasets:
        if dataset.compare_to_explicit_tax_name:
            compare_to = tax_name
        else:
            compare_to = lineage_most_specific
        if dataset.description_pattern.match(compare_to):
            matches += [dataset]
    return matches


def pair_each_sequence_with_dataset(sequences, datasets, datasets_to_skip,
                                    allow_multiple_dataset_matches):
    matches = {}
    for sequence in sequences:
        # there is sequence.taxonomy_name, but the dataset descriptions
        # are based on the most specific label in the lineage
        # (sequence.lineage), which is often (but not always) the same as
        # sequence.taxonomy_name
        sequence_lineage_most_specific = sequence.lineage.split(',')[-1]
        matching_datasets = find_datasets_for_name(datasets,
                                sequence_lineage_most_specific,
                                sequence.taxonomy_name)
        if len(matching_datasets) < 1:
            raise ValueError("No matching datasets for %s" % sequence.lineage)

        matching_datasets = [d for d in matching_datasets if d.name not in \
                             datasets_to_skip]
        if len(matching_datasets) == 0:
            # skip this dataset
            continue
        if not allow_multiple_dataset_matches and len(matching_datasets) > 1:
            raise ValueError("More than one matching dataset for %s" %
                             sequence.lineage)
        matches[sequence] = matching_datasets
    return matches


def read_dataset_skip_list(fn):
    datasets_to_skip = set()
    with open(fn) as f:
        for line in f:
            datasets_to_skip.add(line.rstrip())
    return datasets_to_skip


def read_extra_sequences_paths(fn):
    extra_sequences_paths = {}
    with open(fn) as f:
        for line in f:
            ls = line.rstrip().split('\t')
            extra_sequences_paths[ls[0]] = ls[1]
    return extra_sequences_paths


def map_dataset_to_sequences(dataset_for_sequence):
    # using map of sequence->[dataset], return map of dataset->[sequence]
    sequences_for_dataset = defaultdict(list)
    for sequence, datasets in dataset_for_sequence.items():
        for dataset in datasets:
            sequences_for_dataset[dataset].append(sequence)
    return dict(sequences_for_dataset)


def download_dataset(dataset, sequences, extra_sequences_path, out_dir):
    print("Starting download for", dataset.name)

    num_sequences_segmented = sum([s.is_segmented for s in sequences])
    if num_sequences_segmented > 0 and num_sequences_segmented != len(sequences):
        # either none or all sequences should be labeled as segmented
        raise ValueError("There are %d sequences and %d are marked as segmented"
                         % (len(sequences), num_sequences_segmented))

    is_segmented = num_sequences_segmented > 0
    segments = set(s.segment for s in sequences)
    # Note that len(segments) may equal 1 even if is_segmented is True if
    # sequences for dataset are labeled with segments (e.g., assigned
    # 'segment X') but only 1 segment's sequence is available among the
    # sequences (i.e., all sequences are assigned 'segment X')

    if is_segmented:
        gb = download_raw_from_genbank(sequences, results_type='gb')
        strains = parse_strain_from_gb_results(gb)
        num_found = sum(True for s in sequences if strains[s.name] != None)
        sequences_for_strain = breakup_sequences_by_strain(sequences, strains,
                                                           segments)

        print(("%s (%d sequences) is segmented with %d sequenced segments; found "
               "strains for %d of these sequences (%d of the sequences "
               "could not be placed in a grouping)") % (dataset.name,
              len(sequences), len(segments), num_found,
              len(sequences_for_strain[None])))

        write_dir = os.path.join(out_dir, 'data', dataset.name)
        if not os.path.exists(write_dir):
            os.makedirs(write_dir)

        header_names = set()

        num_genomes = 0
        # make a fasta for each strain (or isolate)
        for strain, strain_sequences in sequences_for_strain.items():
            if strain != None:
                num_genomes += 1
                added = make_fasta_for_genomes(strain_sequences, write_dir)
                header_names.update(added)

        # make a separate fasta for each sequence that could not be grouped
        # with a strain
        for sequence in sequences_for_strain[None]:
            num_genomes += 1
            added = make_fasta_for_genomes([sequence], write_dir)
            header_names.update(added)

        num_sequences = len(sequences)
    else:
        print("%s (%d sequences) is not segmented" % (dataset.name,
              len(sequences)))

        write_dir = os.path.join(out_dir, 'data')
        if not os.path.exists(write_dir):
            os.makedirs(write_dir)

        # make a fasta for this dataset
        header_names = make_fasta_for_genomes(sequences, write_dir, dataset.name)
        num_genomes = len(sequences)
        num_sequences = len(sequences)

    if extra_sequences_path:
        sequences_added, genomes_added = merge_with_extra_sequences(dataset,
                                            sequences, extra_sequences_path,
                                            is_segmented, write_dir)
        num_sequences += sequences_added
        num_genomes += genomes_added

    make_dataset_python_file(dataset, num_genomes, num_sequences, segments,
                             out_dir)
        

def breakup_sequences_by_strain(sequences, strains, segments):
    # strains map sequence_name->strain; make and return a map of
    # strain->[sequences]
    # in the returned dict, None->[sequences with no identified strain]

    sequence_for_name = {s.name: s for s in sequences}
    sequences_for_strain = defaultdict(list)
    for sequence_name, strain in strains.items():
        sequence = sequence_for_name[sequence_name]
        sequences_for_strain[strain].append(sequence)

    # ensure each strain has at most 1 of each segment
    strain_additions_to_none = set()
    seq_additions_to_none = set()
    for strain in sequences_for_strain.keys():
        if strain == None:
            continue
        for segment in segments:
            num_occurrences_of_segment = sum(True for g in
                sequences_for_strain[strain] if g.segment == segment)
            if num_occurrences_of_segment > 1:
                # segment appears more than once in this strain, so
                # invalidate it by adding all of these sequences to the
                # 'None' strain
                strain_additions_to_none.add(strain)
                for s in sequences_for_strain[strain]:
                    seq_additions_to_none.add(s)
    for strain in strain_additions_to_none:
        del sequences_for_strain[strain]
    for s in seq_additions_to_none:
        sequences_for_strain[None].append(s)

    sequences_for_strain = dict(sequences_for_strain)
    if None not in sequences_for_strain:
        sequences_for_strain[None] = []

    return sequences_for_strain


def merge_with_extra_sequences(dataset, sequences, extra_sequences_path,
                               is_segmented, write_dir):
    sequences_added = 0
    genomes_added = 0

    acc_nums_present = set([s.name for s in sequences])

    if is_segmented:
        # since this dataset is segmented, the extra sequences should be
        # in individual genome files and the given path should be a directory
        assert os.path.isdir(extra_sequences_path)

        for fn in os.listdir(extra_sequences_path):
            # check whether none or all of the sequences in fn are 'extra'
            # (i.e., are needed to be copied)
            acc_nums_needed = set()
            num_sequences = 0
            with open(os.path.join(extra_sequences_path, fn)) as f:
                for line in f:
                    if line.startswith('>'):
                        header = line.rstrip()[1:]
                        acc_num = extract_accession_num_from_header(header)
                        if (acc_num not in acc_nums_present and
                                acc_num.split('.')[0] not in acc_nums_present):
                            # acc_num nor the num without the verison number
                            # is already present, so we need it
                            acc_nums_needed.add(acc_num)
                        num_sequences += 1

            num_needed = len(acc_nums_needed)
            if num_needed == 0:
                continue
            if num_needed != num_sequences:
                raise ValueError(("In %s, some but not all sequences are "
                                  "needed; unsure what to do") % 
                                 os.path.join(extra_sequences_path, fn))

            # find the segment needed based on a header
            def extract_segment(header):
                segment_match = re.search(
                    'segment (.+?)(?: |,)|\|Segment:(.+?)\|', header)
                if segment_match is None:
                    raise ValueError("In %s, could not determine a segment" %
                                     os.path.join(extra_sequences_path, fn))

                if segment_match.group(1):
                    return segment_match.group(1)
                else:
                    return segment_match.group(2)

            # copy over fn while changing each header to include the segment
            # at the end (i.e., 'header [segment X]') as a suffix
            # (but not if it already ends in the suffix)
            print("Copying genome to", os.path.join(write_dir, fn))
            with open(os.path.join(write_dir, fn), 'w') as fw:
                with open(os.path.join(extra_sequences_path, fn)) as fr:
                    for line in fr:
                        if line.startswith('>'):
                            header = line.rstrip()[1:]
                            acc_num = extract_accession_num_from_header(header)
                            acc_nums_present.add(acc_num)
                            segment = extract_segment(header)
                            suffix = ' [segment ' + segment + ']'
                            if header.endswith(suffix):
                                updated_header = header
                            else:
                                updated_header = header + suffix
                            fw.write('>' + updated_header + '\n')
                        else:
                            fw.write(line)

            sequences_added += num_needed
            genomes_added += 1
    else:
        # since this dataset is not segmented, the extra sequences should
        # be in a single fasta file
        assert os.path.isfile(extra_sequences_path)

        out_path = os.path.join(write_dir, dataset.name + '.fasta')

        # append sequences to out_path
        with open(out_path, 'a') as fw:
            fw.write('\n')
            with open(extra_sequences_path) as fr:
                currently_appending = False
                for line in fr:
                    if line.startswith('>'):
                        header = line.rstrip()[1:]
                        acc_num = extract_accession_num_from_header(header)
                        if (acc_num not in acc_nums_present and
                                acc_num.split('.')[0] not in acc_nums_present):
                            # append the current sequence
                            currently_appending = True
                            sequences_added += 1
                            genomes_added += 1
                            acc_nums_present.add(acc_num)
                        else:
                            currently_appending = False
                    if currently_appending:
                        fw.write(line)

    return sequences_added, genomes_added


def make_fasta_for_genomes(sequences, write_dir, name=None,
                           max_tries=5,
                           base_delay=5):
    # GenBank sporadically appears to give back a FASTA with missing
    # sequences (i.e., does not include all that were requested);
    # if this is the case, a ValueError is thrown and retry a few times
    # before crashing the entire program
    try_num = 1
    while try_num <= max_tries:
        try:
            return _make_fasta_for_genomes(sequences, write_dir, name=name)
        except ValueError as e:
            if try_num == max_tries:
                # used up all tries
                raise e
            time.sleep(2**(try_num - 1) * base_delay)
            try_num += 1


def _make_fasta_for_genomes(sequences, write_dir, name=None):
    if name is None:
        # make a name from the hash of the sequence names
        name = hashlib.sha224(''.join([s.name for s in sequences]).encode()).\
                    hexdigest()[-8:]

    out_path = os.path.join(write_dir, name + '.fasta')
    fasta_txt = download_raw_from_genbank(sequences, results_type='fasta')

    # Map header->sequence (sequence_for_header)
    # Rather than doing this while reading/copying the fasta file, this
    # allows us to ensure every sequence has an associated header in the
    # return from GenBank
    header_names = set(line[1:] for line in fasta_txt.split('\n')
                       if line.startswith('>'))
    sequence_for_header = {}
    for s in sequences:
        found = None
        for header in header_names:
            if s.name in header:
                if found is not None:
                    raise ValueError(("Found more than one possible header "
                                      "for sequence %s") % s.name)
                found = header
        if found is None:
            raise ValueError("Could not find header for sequence %s" % s.name)
        sequence_for_header[found] = s

    with open(out_path, 'w') as f:
        for line in fasta_txt.split('\n'):
            if line.startswith('>'):
                header = line.rstrip()[1:]
                if header not in sequence_for_header:
                    raise ValueError("Unknown sequence for header %s" % header)
                sequence = sequence_for_header[header]
                if sequence.is_segmented:
                    header = header + ' [' + sequence.segment + ']'
                f.write('>' + header)
            else:
                f.write(line)
            f.write('\n')

    return header_names


def download_raw_from_genbank(sequences,
                              results_type='fasta',
                              batch_size=50,
                              max_tries=5,
                              base_delay=5):
    # Entrez gives sporadic exceptions (usually RuntimeErrors);
    # retry the call a few times if this happens before crashing
    # the entire program
    try_num = 1
    while try_num <= max_tries:
        try:
            return _download_raw_from_genbank(sequences,
                                              results_type=results_type,
                                              batch_size=batch_size)
        except Exception as e:
            if try_num == max_tries:
                # used up all tries
                raise e
            time.sleep(2**(try_num - 1) * base_delay)
            try_num += 1


def _download_raw_from_genbank(sequences,
                               results_type='fasta',
                               batch_size=50):
    # download raw data from GenBank of type 'gb' or 'fasta', as specified
    # by results_type

    accession_names = [s.name for s in sequences]

    # first fetch GI numbers in batches of size batch_size_large
    gi = []
    for i in range(0, len(accession_names), batch_size):
        gi_batch_query = ' '.join(accession_names[i:(i + batch_size)])
        gi_batch = Entrez.read(Entrez.esearch(db='nuccore',
                                              term=gi_batch_query,
                                              retmax=10**6))['IdList']
        gi.extend(gi_batch)

    raw_results = ''

    # now query GenBank using the GI numbers, fetching results in batches
    # of size batch_size
    query = ','.join(gi)
    reader = Entrez.read(Entrez.epost(db='nuccore', id=query))
    for i in range(0, len(gi), batch_size):
        results = Entrez.efetch(db='nuccore',
                                rettype=results_type,
                                retstart=i,
                                retmax=batch_size,
                                webenv=reader['WebEnv'],
                                query_key=reader['QueryKey'])
        raw_results += results.read()

    return raw_results


def parse_strain_from_gb_results(gb_results):
    # gb_results is output of download_raw_from_genbank with results_type='gb'
    # returns map of accession_name->strain

    # each result is separated by the line '//', so split on this
    gb_results_split = gb_results.split('//\n')

    accession_pattern = '^ACCESSION\s+(\w+)( |$)'
    strain_pattern = '/strain="(.+?)"'
    isolate_pattern = '/isolate="(.+?)"'

    strains = {}

    for result in gb_results_split:
        if len(result) == 0 or result.isspace():
            continue

        # find accession number
        accession_match = re.search(accession_pattern, result, re.MULTILINE)
        if not accession_match:
            print(accession_pattern, result)
            raise ValueError("Unknown accession number in result")
        accession = accession_match.group(1)

        # find strain
        strain_match = re.search(strain_pattern, result, re.MULTILINE)
        if strain_match:
            strain = strain_match.group(1)
        else:
            # 'strain' is not present, so look for isolate
            isolate_match = re.search(isolate_pattern, result, re.MULTILINE)
            if isolate_match:
                strain = isolate_match.group(1)
            else:
                # no strain or isolate
                strain = None
        strains[accession] = strain

    return strains


def make_dataset_python_file(dataset, num_genomes, num_sequences,
                             segments, out_dir):
    # sort segments alphanumerically (numerically if they're numbers)
    segments = sorted(segments,
                      key=lambda x: (int(x) if x.isdigit() else float('inf'), x))

    is_segmented = len(segments) > 1
    if is_segmented:
        template_file = DATASET_PYTHON_TEMPLATE_SEGMENTED
    else:
        template_file = DATASET_PYTHON_TEMPLATE_UNSEGMENTED

    fillins = {}
    fillins['VIRUS_REGEX'] = dataset.description
    fillins['NUM_GENOMES'] = str(num_genomes)
    fillins['DATASET_NAME'] = dataset.name

    if dataset.subset_of is None:
        fillins['SUBSET_NOTE'] = ''
    else:
        fillins['SUBSET_NOTE'] = ('\n' + "Note that the sequences in this "
                                  "dataset are a subset of those in the '%s' "
                                  "dataset." % dataset.subset_of + '\n')

    if is_segmented:
        segments_brief = [s.replace('segment ', '') for s in segments]
        fillins['NUM_SEGMENTS'] = str(len(segments))
        fillins['NUM_SEQUENCES'] = str(num_sequences)
        fillins['LIST_OF_SEGMENTS_PYTHON_FORM'] = \
            '[' + ', '.join(["'" + s + "'" for s in segments_brief]) + ']'
        fillins['LIST_OF_SEGMENTS_REGEX_FORM'] = '|'.join(segments_brief)

    out_file = os.path.join(out_dir, dataset.name + '.py')

    # Use a custom instance of TextWrapper so that we can set
    # replace_whitespace to False and therefore allow '\n' in the
    # values of fillins
    tw = textwrap.TextWrapper(replace_whitespace=False)

    with open(out_file, 'w') as fw:
        with open(template_file) as fr:
            num_comment_delims = 0
            for line in fr:
                line = line.rstrip()
                num_comment_delims += line.count('"""')
                line_filledin = line
                for k, v in fillins.items():
                    line_filledin = line_filledin.replace('[[' + k + ']]', v)
                if num_comment_delims % 2 == 1:
                    # in multiline comment
                    has_final_newline = line_filledin.endswith('\n')
                    line_filledin = tw.fill(line_filledin)
                    if has_final_newline:
                        # fill() always strips final newlines, so if it had
                        # one add it back
                        line_filledin += '\n'
                fw.write(line_filledin + '\n')


def read_extra_sequences_headers(extra_sequences_path, is_segmented):
    headers = []
    if is_segmented:
        # since this dataset is segmented, the extra sequences should be
        # in individual genome files and the given path should be a directory
        assert os.path.isdir(extra_sequences_path)

        for fn in os.listdir(extra_sequences_path):
            with open(os.path.join(extra_sequences_path, fn)) as f:
                for line in f:
                    if line.startswith('>'):
                        header = line.rstrip()[1:]
                        headers += [header]
    else:
        # since this dataset is not segmented, the extra sequences should
        # be in a single fasta file
        assert os.path.isfile(extra_sequences_path)

        with open(extra_sequences_path) as fr:
            for line in fr:
                if line.startswith('>'):
                    header = line.rstrip()[1:]
                    headers += [header]
    return headers

def extract_accession_num_from_header(header):
    acc_match = re.search(
        '\|(?:gb|emb|dbj|ref)\|(.+?)\||gb:(.+?)\|', header)
    if acc_match is None:
        raise ValueError("In %s, could not determine accession" %
                         header)

    if acc_match.group(1):
        return acc_match.group(1)
    else:
        return acc_match.group(2)

def write_accession_nums(dataset, sequences, extra_sequences_path, out_dir):
    segments = set(s.segment for s in sequences)
    is_segmented = len(segments) > 1

    accession_nums = set([s.name for s in sequences])

    if extra_sequences_path:
        headers = read_extra_sequences_headers(extra_sequences_path,
                                               is_segmented)
        for header in headers:
            acc_num = extract_accession_num_from_header(header)
            if acc_num in accession_nums:
                # skip because we already have it
                continue
            if acc_num.split('.')[0] in accession_nums:
                # skip because we already have the prefix of the accession
                # number (i.e., the number without the version); e.g.,
                # acc_num is 'KJ123.1' but we have 'KJ123'
                continue
            accession_nums.add(acc_num)

    accession_nums = sorted(list(accession_nums))
    out_file = os.path.join(out_dir, dataset.name + '.accession_nums')
    with open(out_file, 'w') as fw:
        for an in accession_nums:
            fw.write(str(an) + '\n')


def main(args):
    datasets = read_dataset_list(args.dataset_list)
    sequences = read_genome_accession_list(args.genome_accession_list)
    sequences = filter_sequences_with_nonhuman_host(sequences)
    sequences = uniqueify_genome_accession_list(sequences)

    if args.datasets_to_skip:
        datasets_to_skip = read_dataset_skip_list(args.datasets_to_skip)
    else:
        datasets_to_skip = set()

    if args.extra_sequences:
        extra_sequences_paths = read_extra_sequences_paths(
            args.extra_sequences)
    else:
        extra_sequences_paths = {}

    verify_dataset_list(datasets)
    verify_sequence_names_are_unique(sequences)

    dataset_for_sequence = pair_each_sequence_with_dataset(sequences, datasets,
        datasets_to_skip, args.allow_multiple_dataset_matches)
    sequences_for_dataset = map_dataset_to_sequences(dataset_for_sequence)

    for dataset in datasets:
        if dataset not in sequences_for_dataset:
            if dataset.name not in datasets_to_skip:
                raise ValueError("No sequences for dataset: %s" % dataset.name)

    if args.skip_download_and_print_sequences:
        for dataset, sequences in sequences_for_dataset.items():
            for s in sequences:
                print('\t'.join([dataset.name, s.representative, s.name, s.lineage]))
        return 0

    for dataset, sequences in sequences_for_dataset.items():
        if dataset.name in extra_sequences_paths:
            extra_sequences_path = extra_sequences_paths[dataset.name]
        else:
            extra_sequences_path = None

        if args.skip_download_and_write_accession_nums:
            write_accession_nums(dataset, sequences, extra_sequences_path,
                                 args.skip_download_and_write_accession_nums)
        else:
            download_dataset(dataset, sequences, extra_sequences_path,
                             args.out_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-dl', '--dataset_list', required=True,
        help="File with list of datasets")
    parser.add_argument('-gl', '--genome_accession_list', required=True,
        help="File with accession list of all viral genomes")
    parser.add_argument('-ds', '--datasets_to_skip',
        help="File with list of datasets to not download")
    parser.add_argument('-es', '--extra_sequences',
        help=("File with list of datasets and corresponding paths to extra "
              "sequences to merge with those downloaded from GenBank"))
    parser.add_argument('--skip_download_and_write_accession_nums',
        help=("When set, do not perform the download and instead write "
              "a list of accession nums for each dataset in the specified "
              "directory"))
    parser.add_argument('--skip_download_and_print_sequences',
        dest="skip_download_and_print_sequences",
        action="store_true",
        help=("When set, do not perform the download and instead print "
              "a list of sequences (accession and lineage) for each "
              "dataset"))
    parser.add_argument('--allow_multiple_dataset_matches',
        dest="allow_multiple_dataset_matches",
        action="store_true",
        help=("When set, do not fail when a sequence matches more than one "
              "dataset"))
    parser.add_argument('-o', '--out_dir', required=True,
        help="Directory in which to place output data")

    args = parser.parse_args()  

    main(args)
