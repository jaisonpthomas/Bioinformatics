def space_sep_list(lst):
    return ' '.join([str(item) for item in lst])


def pattern_count(text, pattern):
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if text[i: i+len(pattern)] == pattern:
            count += 1
    return count


def frequent_words(text, k):
    from collections import Counter
    word_counts = Counter()
    for i in range(len(text)-k+1):
        word_counts[text[i: i+k]] += 1
    return word_counts


def reverse_comp(seq):
    comps = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
    output = []
    for ch in seq:
        output.append(comps[ch])
    output.reverse()
    return ''.join(output)


def pat_occurs(genome, pat):
    positions = []
    for i in range(len(genome)-len(pat)+1):
        if genome[i: i+len(pat)] == pat:
            positions.append(i)
    return positions


def get_clumps(seq, k, L, t):
    clumps = {}

    slice_kmer_positions = {}
    # Pass thru first L-snippet to create cache table
    initial_kmer = seq[0: k]
    for pos in range(L-k+1):
        curr_kmer = seq[pos: pos+k]
        if curr_kmer not in slice_kmer_positions:
            slice_kmer_positions[curr_kmer] = set()
        slice_kmer_positions[curr_kmer].add(pos)

    for kmer, pos_set in slice_kmer_positions.items():
        if len(pos_set) >= t:
            clumps[kmer] = pos_set.copy()

    # Grabs L-nuc seqs
    for pos in range(1, len(seq)-L+1):
        slice_kmer_positions[initial_kmer].remove(
            min(slice_kmer_positions[initial_kmer]))
        initial_kmer = seq[pos: pos+k]
        last_kmer = seq[pos+L-k: pos+L]

        if last_kmer not in slice_kmer_positions:
            slice_kmer_positions[last_kmer] = set()
        slice_kmer_positions[last_kmer].add(pos+L-k)

        if len(slice_kmer_positions[last_kmer]) >= t:
            if last_kmer not in clumps:
                clumps[last_kmer] = set()
            clumps[last_kmer] = clumps[last_kmer] | slice_kmer_positions[last_kmer]

    return clumps


def pattern_to_number(seq):
    nuc_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    num_list = []
    for nuc in seq:
        num_list.append(nuc_map[nuc])
    num_list.reverse()

    output = 0
    exp = 0
    for num in num_list:
        output += num * 4**exp
        exp += 1

    return output


def number_to_pattern(num, k):
    rem_list = []
    nuc_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    for _ in range(k):
        num, rem = divmod(num, 4)
        rem_list.append(rem)
    rem_list.reverse()
    for idx, val in enumerate(rem_list):
        rem_list[idx] = nuc_map[val]
    return ''.join(rem_list)


def skew(genome):
    output = [0]
    for i, nuc in enumerate(genome):
        if nuc == 'C':
            output.append(output[i]-1)
        elif nuc == 'G':
            output.append(output[i]+1)
        else:
            output.append(output[i])
    return output


def find_min(skew):
    return [idx for idx, val in enumerate(skew) if val == min(skew)]


def hamming_distance(seq_a, seq_b):
    hd = 0
    for i in range(len(seq_a)):
        if seq_a[i] != seq_b[i]:
            hd += 1
    return hd


def pat_match_hamming(pattern, genome, hd):
    output = []
    k = len(pattern)
    for i in range(len(genome)-k+1):
        gen_subseq = genome[i: i+k]
        if hamming_distance(pattern, gen_subseq) <= hd:
            output.append(i)
    return output


def freq_words_with_mismatch(genome, k, d, rc):
    def freq_words(genome, k, d):

        def add_neighbor(neighbor, next_cache, output):
            if neighbor not in next_cache:
                next_cache.add(neighbor)
                if neighbor.lower() not in output:
                    output[neighbor.lower()] = 0
                output[neighbor.lower()] += 1

        def build_cache(prefix, suffix, edits):
            cache = set()

            def recursive_builder(prefix, suffix, edits):
                if suffix and edits:
                    recursive_builder(prefix+suffix[0], suffix[1:], edits)
                    for nuc in "ACGT":
                        if suffix[0] != nuc.lower():
                            recursive_builder(prefix+nuc, suffix[1:], edits-1)
                else:
                    cache.add(prefix+suffix)
            recursive_builder(prefix, suffix, edits)
            return cache

        genome = genome.lower()
        first_kmer = genome[: k]
        curr_cache = build_cache('', first_kmer, d)
        output = {item.lower(): 1 for item in curr_cache}

        for i in range(len(genome)-k):
            next_cache = set()
            for key in curr_cache:
                next_key_base = key[1:]
                num_edits = len([ch for ch in next_key_base if ch.isupper()])
                next_letter = genome[i+k]
                if num_edits == d:
                    neighbor = next_key_base+next_letter
                    add_neighbor(neighbor, next_cache, output)
                elif num_edits < d:
                    for nuc in "AGCT":
                        if next_letter == nuc.lower():
                            neighbor = next_key_base+next_letter
                        else:
                            neighbor = next_key_base+nuc
                        add_neighbor(neighbor, next_cache, output)
            curr_cache = next_cache
        return output

    output = freq_words(genome, k, d)
    if rc:
        rc_genome = reverse_comp(genome)
        rc_pass = freq_words(rc_genome, k, d)
        output = {k: output.get(k, 0) + rc_pass.get(k, 0)
                  for k in set(output) | set(rc_pass)}

    max_occurs = max(output.values())
    return [(k, v) for k, v in output.items() if v == max_occurs]
