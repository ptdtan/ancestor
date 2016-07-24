MIN_GAP_SIZE = 11
#from ultilities import *
from ultilities import *
from collections import defaultdict

class Contig:
    def __init__(self, name=None, seq=None, start = 0, end = 0 ,
                 link=MIN_GAP_SIZE, sign="+", region=(0,0)):
        self.name = name
        self.seq = seq
        self.link = link
        self.start = start
        self.region = region
        if seq:
            self.end = len(seq)
        else:
            self.end = end
        if sign == "+":
            self.sign = 1
        else:
            self.sign = -1

    def __hash__(self):
        return hash(self.name)

class Scaffold:
    def __init__(self, name, contigs):
        self.contigs = contigs
        self.name = name
        self.hash_cnts = {cnt:self.contigs.index(cnt) for cnt in contigs}

    def __hash__(self):
        return hash(self.name)

    def _join(self, cnt, est_gap):
        self.contigs.append(cnt)
        self.contigs[-2].link = est_gap
        return

    def _break(self, cnt, slide=1):
        """
        break to the rights
        """
        cnt_idx = self.hash_cnts[cnt]
        cnts1 = self.contigs[:cnt_idx+slide]
        cnts2 = self.contigs[cnt_idx+slide:]
        name1 = "%s.broken.%s.1" %(self.name, cnt.name)
        name2 = "%s.broken.%s.2" %(self.name, cnt.name)
        return Scaffold(name1, cnts1), Scaffold(name2, cnts2)

    def _add_seq(self, seqDict):
        for i in range(len(self.contigs)):
            start, end = self.contigs[i].region
            if start + end == 0:
                self.contigs[i].seq = seqDict[self.contigs[i].name][:]
            else:
                self.contigs[i].seq = seqDict[self.contigs[i].name][start:end]

class Assembly:
    def __init__(self, name, scaffolds = []):
        self.scaffolds = scaffolds
        self.name = name
        self.scf_hash = {scf:self.scaffolds.index(scf) for scf in scaffolds}
        self.seqs = None
    @staticmethod
    def with_links(name, links):
        return Assembly(name, scaffolds = parse_links(links))

    def _update(self, old=[], new=[]):
        map(lambda x: self.scaffolds.pop(x), [self.scf_hash[scf] for scf in old])
        map(lambda x: self.scaffolds.append(x), [self.scf_hash[scf] for scf in new])
        self.scf_hash = {scf:self.scaffolds.index(scf) for scf in scaffolds}
        return

    def _getSeq(self, fileSeq):
        self.seqs = fasta_parser(fileSeq)
        pass

def parse_links(links):
    """Parser for scaffolds_links file
    @param ifile: _scaffolds.links file
    @return
    """
    ofile = open(links)
    scaffolds_links = ofile.read().strip().split('\n\n')
    assembly = list()
    for block in scaffolds_links:
    	arr = block.split('\n')
        contigs = []
    	for cnt_raw in arr[1:]:
            try:
                name, start, end , gap, _ = cnt_raw.strip().split("\t")
            except ValueError:
                name, start, end , gap = cnt_raw.strip().split("\t")
            if "[" in name:
                raw_region = name[name.index("[")+1:name.index("]")]
                region = tuple(map(int, raw_region.split(":")))
            contigs.append(Contig(name=name[1:], start = int(start), end = int(end),
                                  link = int(gap), sign=name[0]))
    assembly.append(Scaffold(name = arr[0], contigs = contigs))
    return assembly
