# ----------Reading k-mers--------------------------
# import the required modules
from Bio import SeqIO  # to handle sequences
from collections import defaultdict  # for the kmer collections

# Step 1: Write code to read in the DNA sequence from the above FASTA files,
# and count the occurences of 14-mers in each of the four files
# (storing these in a dictionary or similar).

# example - making kmers from a sequence
seq = "TTGAAAGAAAAACAATTTTG"
k_dict = defaultdict(int)  # empty dict
for idx in range(0, len(seq) - 13):  # last kmer starts at the 14th nt
    kmer = seq[idx:(idx + 14)]
    if kmer in k_dict:  # increase count of kmer exists in the dict
        k_dict[kmer] = k_dict.get(kmer, 0) + 1
    else:
        k_dict[kmer] = 1  # add kmer as key to the dict
print(len(k_dict))


# function to count the kmers from a fasta file
def count_kmers_from_fasta(fasta_file, k):
    # read the input file
    fasta_entry = SeqIO.parse(open(fasta_file), 'fasta')
    # make empty dict
    kmers_dict = defaultdict(int)
    # read each contig separately
    for contig in fasta_entry:
        sequence = str(contig.seq)
        # obtain kmers across the contig seq
        for idx in range(0, len(sequence) - (k - 1)):  # last kmer is the last k nt in the sequence
            kmer = sequence[idx:(idx + k)]
            # store kmers to a dict if doesn't already exist, but increase count if it already exists
            if kmer in kmers_dict:
                kmers_dict[kmer] = kmers_dict.get(kmer, 0) + 1
            else:
                kmers_dict[kmer] = 1
    # count the occurences of kmers the file
    return f'{fasta_file}: {len(kmers_dict)} kmers'


# count the occurences of 14-mers in each of the four files
files = ("s_pneumoniae_genomes/14412_3#82.contigs_velvet.fa",
         "s_pneumoniae_genomes/14412_3#84.contigs_velvet.fa",
         "s_pneumoniae_genomes/R6.fa", "s_pneumoniae_genomes/TIGR4.fa")
for file in files:
    print(count_kmers_from_fasta(file, 14))

# --------------- Simple Jaccard distances --------------------
# calculate the Jaccard distance between two of the input samples
# example of two kmer lists

A = ["TTGAAAGAAAAACA", "TGAAAGAAAAACAA",
     "GAAAGAAAAACAAT", "AAAGAAAAACAATT",
     "AAGAAAAACAATTT", "AGAAAAACAATTTT",
     "GAAAAACAATTTTG"]
B = ["TTGAAAGAAAAACA", "TGAAAGAAAAACAA",
     "GAAAGAAAAACAAT", "AAAGAAAAACAATT",
     "GAAAAACAATTTTG", "CTCGATCCATGTAT",
     "TCGATCCATGTATG"]
# could compare the lists but let's first make dicts as in the first example
# first make empty dicts, then add the kmers
dictA = defaultdict(int)
for kmer in A:
    # check if kmer is already in the dict
    if kmer in dictA:
        dictA[kmer] = dictA.get(kmer, 0)+1 # just increase the value count if already exists
    else:
        dictA[kmer] = 1 # add it if it doesn't

dictB = defaultdict(int)
for kmer in B:
    # check if kmer is already in the dict
    if kmer in dictB:
        dictB[kmer] = dictB.get(kmer, 0)+1 # just increase the value count if already exists
    else:
        dictB[kmer] = 1 # add it if it doesn't

# count shared and unique kmers between samples
AnB, AuB = 0, 0 # initiate count as 0 for both intersection and union
for k in dictA: # compare all A against B
    if k in dictB: # if shared between the two
        AnB += 1 # increase intersection count
        AuB += 1 # increase the union as well
        del dictB[k] # remove the shared kmer from the dict
    else:
        AuB += 1 # it is unique to A, so add it the union
# the remains in B are unique to B, add them to the union as well
AuB += len(dictB)
print(f'AnB: {AnB}, AuB: {AuB}')

# calculate the jaccard index
j_index = AnB/AuB
print(f'jaccard index : {j_index}')

# calculate the jaccard distance
j_dist = 1-j_index
print(f'jaccard distance: {j_dist}')


# Step:2 Using your dictionaries of 14-mer counts, write code which
# calculates the Jaccard distance between two of the input samples.

# on the sample files
# i like to use functions
def fasta_to_kmers(fasta, k_value):
    fasta_entry = SeqIO.parse(open(fasta), 'fasta')
    kmers_dict = defaultdict(int)
    # treat each contig as separate
    for contig in fasta_entry:
        sequence = str(contig.seq)
        # obtain kmers across the contig
        for idx in range(0, len(sequence) - (k_value - 1)):  # last kmer is the last k nt in the sequence
            kmer = sequence[idx:(idx + k_value)]
            # store kmers to a dict if doesn't already exist, but increase count if it already exists
            if kmer in kmers_dict:
                kmers_dict[kmer] = kmers_dict.get(kmer, 0) + 1
            else:
                kmers_dict[kmer] = 1
    return kmers_dict


def cal_jaccard_distance(fasta1, fasta2, k):
    # read the input files and obtain kmers
    kmers_dictA = fasta_to_kmers(fasta1, k)
    kmers_dictB = fasta_to_kmers(fasta2, k)

    # count shared and unique kmers between samples
    AnB, AuB = 0, 0
    for kmer in kmers_dictA:
        if kmer in kmers_dictB:  # shared kmers between the two samples
            AnB += 1
            AuB += 1
            del kmers_dictB[kmer]  # remove the shared kmer from the B dict
        else:  # unique kmers from A
            AuB += 1
    # add the unique kmers in B to the union
    AuB += len(kmers_dictB)

    # calculate the jaccard index
    j_index = AnB / AuB
    print(f'jaccard index : {j_index}')

    # calculate the jaccard distance
    j_dist = 1 - j_index
    #     print(f'jaccard distance: {j_dist}')
    return f'Jaccard Distance : {j_dist}'


cal_jaccard_distance("s_pneumoniae_genomes/14412_3#82.contigs_velvet.fa",
                     "s_pneumoniae_genomes/14412_3#84.contigs_velvet.fa", 14)
cal_jaccard_distance("s_pneumoniae_genomes/R6.fa", "s_pneumoniae_genomes/TIGR4.fa", 14)

# -------------------- MinHash Jaccard distances ------------------------
# Step 3: Find an implementation of a hash function and import it into your code.
# Calculate the hash of some example 14-mers, and confirm that the same 14-mer input maps to the same integer output.
import xxhash

# some 14-mers for this example
kmers = ["TTGAAAGAAAAACA", "TGAAAGAAAAACAA", "GAAAGAAAAACAAT",
         "AAAGAAAAACAATT", "AAGAAAAACAATTT", "AGAAAAACAATTTT",
         "GAAAAACAATTTTG"]
# calculate a hash of the first k-mer, this comes back as a hex string.
print(xxhash.xxh32(kmers[0]).hexdigest())

# Run on all the k-mers
hash1 = []
for kmer in kmers:
    hash1.append(xxhash.xxh32(kmer).hexdigest())
print(hash1)

# Confirm the same hashes are generated on repeated runs (so they are comparable between samples)
hash2 = []
for kmer in kmers:
    hash2.append(xxhash.xxh32(kmer).hexdigest())
print(hash2 == hash1)

# Step 4: Make a sketch of each of the input sequences by following the steps above:
# mapping any present 14-mers (ignoring their actual count as before)
# to integers using a hash function, then taking the lowest 1000 values and storing these in a list/array.
# Save your sketches for each input sequence: either as a text file, or some other representation.

# get dictionary of kmers
dict82 = fasta_to_kmers("s_pneumoniae_genomes/14412_3#82.contigs_velvet.fa", 14)
dict84 = fasta_to_kmers("s_pneumoniae_genomes/14412_3#84.contigs_velvet.fa", 14)


# make hashes for the kmers
def dict_hashes(the_dict):
    result = []
    for key in the_dict:
        result.append(xxhash.xxh32(key).hexdigest())
    return result


# loop through the files
files = ("s_pneumoniae_genomes/14412_3#82.contigs_velvet.fa", "s_pneumoniae_genomes/14412_3#84.contigs_velvet.fa")
for file in files:
    dict_file = fasta_to_kmers(file, 14) # use previous function to get kmers from the fasta
    hashes = dict_hashes(dict_file) # get hashes for the kmers using the function
    sorted_hashes = sorted(hashes) # sort the hashes
    s = 1000
    smallest = sorted_hashes[0:1000] # get the smallest 1000 hashes
    # save the smallest hashes to a text file
    text = open(f'{file[0:-3]}_1000smallest_hashes', "w") # name text file accordingly
    for elt in smallest:
        text.write(elt + "\n")
    text.close()

# Step 5: For each pair of inputs, load their sketches, and
# calculate their Jaccard distance using the sizes of the intersection and the union of the saved hash values.
# load their sketches
# calculate their Jaccard distance using the sizes of the intersection and the union of the saved hash values


# function to calculate haccard distance from a sketch of hashes
def cal_hashes_jaccard_distance(sketch1, sketch2):
    # read the input files for the hashes
    sketchA = open(sketch1, "r")
    sketchB = open(sketch2, "r")
    Ahashes = sketchA.readlines()
    Bhashes = sketchB.readlines()

    # count shared and unique hashes between the two inputs
    AnB, AuB = 0, 0
    for hsh in Ahashes:
        if hsh in Bhashes:  # shared hashes between the two inputs
            AnB += 1
            AuB += 1
            Bhashes.remove(hsh)  # remove the shared hash from B
        else:  # unique hashes from A
            AuB += 1
    # add the unique (those that remained) hashes in B to the union
    AuB += len(Bhashes)

    # calculate the jaccard index
    j_index = AnB / AuB
    print(f'jaccard index : {j_index}')

    # calculate the jaccard distance
    j_dist = 1 - j_index
    #     print(f'jaccard distance: {j_dist}')
    return j_dist


cal_hashes_jaccard_distance("s_pneumoniae_genomes/14412_3#82.contigs_velvet_1000smallest_hashes",
                            "s_pneumoniae_genomes/14412_3#84.contigs_velvet_1000smallest_hashes")
