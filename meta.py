import os, re, csv
from typing import Dict, List, Tuple

def IOHandler(dir_path: str, file_path: str) -> Tuple[Dict[int, str], Dict[int, str]]:
    genome_sequences = {}

    for file_name in os.listdir(dir_path):
        with open(os.path.join(dir_path, file_name), 'r') as file:
            genome_number = None
            sequence = []
            for line in file:
                if line.startswith('>'):
                    match = re.search(r'Genome_Number_(\d+)', line)
                    if match:
                        genome_number = int(match.group(1))
                else:
                    sequence.append(line.strip())
            if genome_number is not None:
                genome_sequences[genome_number] = ''.join(sequence)

    read_sequences = {}
    current_read_number = None
    current_sequence = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                if current_read_number is not None:
                    read_sequences[current_read_number] = ''.join(current_sequence)
                current_read_number = re.search(r'read_(\d+)', line).group(1)
                current_sequence.clear()
            else:
                current_sequence.append(line.strip())

    if current_read_number is not None:
        read_sequences[current_read_number] = ''.join(current_sequence)

    return genome_sequences, read_sequences

def reader(index: Dict[int, str], r) -> Dict[int, int]:
    matches = {}
    if len(r) > 45:
        r = r[:45]

    def helper(key: str):
        if key in index:
            for i in index[key]:
                if i not in matches:
                    matches[i] = 0
                matches[i] += 1

    helper(r[:15])
    helper(r[15:30])
    helper(r[30:])

    return matches
    
def finder(index: Dict[int, str], read: str) -> int:
    matches = reader(index, read)
    if len(matches) == 1:
        for key in matches.keys():
            return key
    elif len(matches) > 1:
        bestmatch = -1
        maxcount = 0
        for genome in matches:
            if matches[genome] > maxcount:
                maxcount = matches[genome]
                bestmatch = genome
        return bestmatch 
    return -1

if __name__ == "__main__":
    reference_genomes_directory = "./1d-references"
    reads_file = "./project1d_reads.fasta"

    refs, reads = IOHandler(reference_genomes_directory, reads_file)
    numGenomes = len(refs)

    k = 15
    fullMatches = {}
    for r in refs:
        for i in range(len(refs[r]) - k):
            kmer = refs[r][i:i + k]
            if kmer not in fullMatches:
                fullMatches[kmer] = []
            fullMatches[kmer].append(r)

    threshold = 0.007
    counts = [0] * numGenomes
    for key, value in reads.items():
        matches = reader(fullMatches, value)
        for m in matches:
            counts[m] += 1
    abovethreshold = []
    total = sum(counts)
    for i in range(len(counts)):
        if counts[i] > threshold * total:
            abovethreshold.append(i)
    partialMatches = {}
    for genome in abovethreshold:
        for i in fullMatches:
            if genome in fullMatches[i]:
                if i not in partialMatches:
                    partialMatches[i] = []
                partialMatches[i].append(genome)

    res = {}
    for key, value in reads.items():
        matchingGenome = finder(partialMatches, value)
        res[key] = matchingGenome
    
    with open("predictions.csv", "w", newline = "") as f:
        writer = csv.writer(f)
        for r in res:
            writer.writerow([f">read_{r} Genome_Number_{res[r]}"])