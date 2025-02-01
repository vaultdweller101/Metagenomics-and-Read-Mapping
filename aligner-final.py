from sys import argv
from Bio import SeqIO
import csv
from typing import List, Tuple, Dict
# susu

def stringDistance(a: str, b: str) -> int:
    res = 0
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]:
            res += 1
    return res

def mapRead(indexDict: Dict[str, int], read: str, referenceGenome: str) -> List[str]:

    def findSub(read: str, referenceGenome: str, index: int) -> List[str]:
        res = []
        for i in range(min(len(read), len(referenceGenome))):
            if read[i] != referenceGenome[i]:
                res.append(">S" + str(index + i) + " " + referenceGenome[i] + " " + read[i])
        return res
    
    def findInsertDeletion(read: str, referenceGenome: str, index: int) -> List[str]:
        res = []
        i = 0
        while i < len(read) and read[i] == referenceGenome[i]:
            i += 1
        i -= 1
        og = stringDistance(read, referenceGenome)
        ifDel = stringDistance(read[:i] + referenceGenome[i] + read[i:]  , referenceGenome)
        ifInsert = stringDistance(read[:i] + read[i + 1:], referenceGenome)
        if ifDel < 2 and ifDel < og and ifDel < ifInsert:
            res.append(">D" + str(index + i) + " " + referenceGenome[i])
        elif ifInsert < 2 and ifInsert < og and ifInsert < ifDel:
            res.append(">I" + str(index + i) + " " + read[i])
        return res

    res = []
    if len(read) > 48:
        read = read[:48]
    if len(read) < 48:
        return res
    matches_set = set()
    threshold = len(referenceGenome) - 48

    def helper(key: str, delta: int):
        if key in indexDict:
            for i in indexDict[key]:
                if i - delta not in matches_set and 0 <= i - delta < threshold:
                    matches_set.add(i - delta)

    helper(read[:16], 0)
    helper(read[16:32], 16)
    helper(read[32:], 32)

    for m in matches_set:
        refmatch = referenceGenome[m : m + 48]
        if stringDistance(read, refmatch) <= 3:
            res.extend(findSub(read, refmatch, m))
        else:
            res.extend(findInsertDeletion(read, refmatch, m))
    return res

def isError(freq, subLimit: int, indelLimit: int) -> List[str]:
    res = []
    for m in freq:
        if (m[1] == "S" and freq[m] >= subLimit) \
            or (m[1] in "ID" and freq[m] >= indelLimit):
            res.append(m)
    return res

def IOHandler(pathRef: str, pathReads: str) -> Tuple[str, List[str]]:
    used = {}
    for i in SeqIO.parse(pathRef, "fasta"):
        used[i.id] = i.seq
    genomes = str(used['genome'])
    
    with open(pathReads) as f:
        reads = [read.strip() for read in f.readlines() if read[0] != '>']
    return (genomes, reads)

if __name__ == "__main__":
    refAddress = '.\\project1b-s_reference_genome.fasta'
    readsAddress = '.\\project1b-s_with_error_paired_reads.fasta'
    ref, reads = IOHandler(refAddress, readsAddress)
    used = {}
    k = 16
    for i in range(len(ref) - k):
        kmer = ref[i:i + k]
        if kmer not in used:
            used[kmer] = []
        used[kmer].append(i)

    freq = {}

    for r in reads:
        for m in mapRead(used, r, ref):
            if m not in freq:
                freq[m] = 0
            freq[m] += 1

    trueMutations = isError(freq, 4, 2)

    with open("predictions.csv", "w", newline = "") as f:
        writer = csv.writer(f)
        for m in trueMutations:
            writer.writerow([m])