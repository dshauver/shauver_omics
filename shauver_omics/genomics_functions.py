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

def naive_Nmm(p, t, err_count):
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
    read_length = len(quals[0])
    num_reads = len(quals)

    # Initialize a list to store the sum of quality scores for each position
    sum_scores = [0] * read_length

    # For each position
    for pos in range(read_length):
        total_score = 0
        for read in quals:
            char = read[pos]
            score = phred33ToQ(char)
            total_score += score
        # Compute the average score for this position
        sum_scores[pos] = total_score / num_reads
    
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

def queryIndex_Nmm(p, t, index, err_count):
    k = index.k
    offsets = []
    for i in index.query(p):
        match = True
        error_counter = 0
        if p[k:] != t[i+k:i+len(p)]:  # verify that rest of P matches
            error_counter += 1
            if error_counter > err_count:
              match = False
              break
        if match:
          offsets.append(i)
    return offsets

def approximate_match(p, t, n):
    segment_length = int(round(len(p) / (n+1))) # find length of segments
    all_matches = set() # use a set for indices to avoid duplicates
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
        matches = boyer_moore(p[start:end], p_bm, t)
        # Extend matching segments to see if whole p matches
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches)

def approximate_match_index(p, t, n, k=8):
    """ Find all approximate matches of pattern p in text t with up to n mismatches
        using Index class for k-mer matching. Also counts total index hits.
        
        Args:
            p: pattern string
            t: text string
            n: maximum number of mismatches allowed
            k: k-mer length for Index (default 8)
            
        Returns:
            Tuple: (List of starting positions where p approximately matches t, total index hits)
    """
    segment_length = int(round(len(p) / (n + 1)))  # find length of segments
    all_matches = set()  # use a set for indices to avoid duplicates
    index = Index(t, k)  # create index for text with k-mer length
    total_hits = 0  # track total index hits
    
    for i in range(n + 1):
        start = i * segment_length
        end = min((i + 1) * segment_length, len(p))
        # Get matches for the current segment using Index
        segment = p[start:end]
        if len(segment) < k:
            continue  # skip segments shorter than k
        matches = index.query(segment[:k])
        total_hits += len(matches)  # count hits for this segment
        
        # Extend matching segments to verify if whole pattern matches
        for m in matches:
            # Adjust match position to account for segment start
            match_start = m - start
            if match_start < 0 or match_start + len(p) > len(t):
                continue
                
            # Count mismatches across the entire pattern
            mismatches = 0
            for j in range(len(p)):
                if p[j] != t[match_start + j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(match_start)
    
    return sorted(list(all_matches)), total_hits

class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def approximate_match_SubSeqIndex(p, t, n, k=8, ival=3):
    """ Find all approximate matches of pattern p in text t with up to n mismatches
        using SubseqIndex class for k-mer matching. Also counts total index hits.
        
        Args:
            p: pattern string to search for
            t: text string to search within
            n: maximum number of mismatches allowed in a match
            k: length of k-mers used in SubseqIndex (default 8)
            ival: interval between sampled characters in subsequences (default 3)
            
        Returns:
            Tuple: (List of starting positions where p approximately matches t, total index hits)
    """
    # Create a SubseqIndex object for the text using the specified k-mer length and interval
    index = SubseqIndex(t, k, ival)
    
    # Calculate the span: the total length of characters covered by a k-mer subsequence
    # Formula: 1 (first character) + ival * (k - 1) (intervals between k characters)
    span = 1 + ival * (k - 1)
    
    # Use a set to collect unique starting positions of approximate matches
    all_matches = set()
    
    # Initialize a counter to track the total number of index hits across all queries
    total_hits = 0
    
    # Iterate over different offsets (0 to ival-1) to cover all possible alignments
    for o in range(ival):
        # Skip this offset if it would make the subsequence longer than the pattern
        if o + span > len(p):
            break
        
        # Query the index with a subsequence of p starting at offset o
        matches = index.query(p[o:])
        
        # Add the number of matches found at this offset to the total hit count
        total_hits += len(matches)
        
        # Process each match position returned by the index query
        for m in matches:
            # Adjust the match position to account for the offset in the pattern
            match_start = m - o
            
            # Ensure the adjusted position is valid within the text bounds
            if match_start >= 0 and match_start + len(p) <= len(t):
                # Initialize a counter for mismatches between pattern and text
                mismatches = 0
                
                # Compare each character of the pattern with the corresponding text character
                for j in range(len(p)):
                    if p[j] != t[match_start + j]:
                        mismatches += 1
                        # Exit early if mismatches exceed the allowed limit
                        if mismatches > n:
                            break
                
                # If mismatches are within the allowed limit, record this position
                if mismatches <= n:
                    all_matches.add(match_start)
    
    # Convert the set of matches to a sorted list and return it along with total hits
    return sorted(list(all_matches)), total_hits