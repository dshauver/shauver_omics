#!python3

def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return (occurrences)

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def naive_with_rc(p, t):
    occurrences = []
    occurrences_rc = []
    occurrences = naive(p, t)
    rc_p = reverseComplement(p)
    if not p == rc_p:
        occurrences_rc = naive(rc_p, t)
    if not occurrences_rc:
        return occurrences
    return occurrences + occurrences_rc

def naive_2mm(p, t, err_count):
    occurences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        error_counter = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                error_counter += 1
                if error_counter > err_count:
                    match = False
                    break        
        if match:
            occurences.append(i)  # all chars matched; record
    return occurences

def phred33ToQ(qual):
    return ord(qual) - 33

def createHist(qualities):
    # Create a histogram of quality scores
    hist = [0]*50
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist

def find_lowqual_reads(quals):
    seq_length = len(quals[0])
    num_seqs = len(quals)

    # Initialize a list to store the sum of quality scores for each position
    sum_scores = [0] * seq_length

    # For each position
    for pos in range(seq_length):
        total_score = 0
        for seq in quals:
            char = seq[pos]
            score = phred33ToQ(char)
            total_score += score
        # Compute the average score for this position
        sum_scores[pos] = total_score / num_seqs
    
    # Find the position with the minimum average score
    min_score = min(sum_scores)
    min_position = sum_scores.index(min_score)

    return(min_score, min_position)

def has_stop_codon(dna, frame=0) :
    "This function checks if given dna sequence has in frame stop codons"
    stop_codon_found = False
    stop_codons = ['tga','tag','taa']
    for i in range(frame, len(dna), 3) :
        codon = dna[i:i+3].lower()
        if codon in stop_codons :
            stop_codon_found = True
            break
    return(stop_codon_found)

def reversecomplement(seq):
    "Return the reverse complement of the dna string."

    seq = reverse_string(seq)

    seq = complement(seq)

    return(seq)

def reverse_string(seq):
    return seq[::-1]

def complement(seq):
    "Return the complementary sequence string"

    basecomplement = {
        'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
        'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'
    }

    letters = list(seq)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)

def gc(seq):
    seq = seq.lower()

    no_c = seq.count('c')

    #count the # of g's
    no_g = seq.count('g')

    #get the lenght of the dna sequence
    seq_length = len(seq)

    #calculate the gc ratio
    gc_percent = (no_c + no_g) * 100.0 / seq_length

    return(gc_percent)