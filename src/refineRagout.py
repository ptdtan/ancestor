"""
framework for refining Ragout's scaffolds by using read-pair mapping information.
@ptdtan
"""

MIN_GAP_SIZE = 11

import ultilities as ul
from collections import defaultdict
from Bio import SeqIO
from copy import copy

class Contig:
    def __init__(self, uname=None, seq="",
                 link=MIN_GAP_SIZE, sign="+", region=(0,0), supported = ""):
        self.uname = uname
        if "[" in uname:
            self.name = uname[:uname.index("[")]
        else:
            self.name = uname
        self.seq = seq
        self.link = link
        self.start = region[0]
        self.end = region[1]
        self.region = region
        if sign == "+":
            self.sign = 1
        else:
            self.sign = -1
        self.supporting_genomes = supported

class Scaffold:
    def __init__(self, name, contigs):
        self.contigs = contigs
        self.name = name
        self.hash_cnts = {cnt.uname:self.contigs.index(cnt) for cnt in contigs}
        self.seq = ''

    def __hash__(self):
        return hash(self.name)

    def _join(self, cnt, est_gap=MIN_GAP_SIZE):
        self.contigs.append(copy(cnt))
        self.contigs[-2].link = est_gap
        return

    def _break(self, cntName, slide=1, trimmed=False):
        """
        break to the rights: slide +1
        break to the left: slide 0
        """
        try:
            cnt_idx = self.hash_cnts[cntName]
            if trimmed:
                cnts1 = self.contigs[:cnt_idx]
                cnts2 = self.contigs[cnt_idx+slide:]
            else:
                cnts1 = self.contigs[:cnt_idx+slide]
                cnts2 = self.contigs[cnt_idx+slide:]
            cnts1[-1].link = 0
            name1 = "%s.broken.%s.1" %(self.name, cntName)
            name2 = "%s.broken.%s.2" %(self.name, cntName)
            return Scaffold(name1, cnts1), Scaffold(name2, cnts2)
        except KeyError:
            print "Can't break, contig not found"
        return

    def _add_seq(self, seqDict):
        for i in range(len(self.contigs)):
            start, end = self.contigs[i].region
            if start + end == 0:
                self.contigs[i].seq = seqDict[self.contigs[i].name][:]
                self.contigs[i].start = 0
                self.contigs[i].end = len(self.contigs[i].seq)
            else:
                self.contigs[i].seq = seqDict[self.contigs[i].name][start:end]
                self.contigs[i].start = start
                self.contigs[i].end = end

    def _to_sequence(self):
        sequence = []
        for cnt in self.contigs:
            if cnt.sign == 1:
                sequence.append(''.join([str(cnt.seq),'N'*cnt.link]))
            else:
                sequence.append(''.join([ul.reverse_complement(str(cnt.seq)), 'N'*cnt.link]))
        self.seq = ''.join(sequence)
        pass

class Assembly:
    """
    basic Usage:
    ass = Assembly.with_links("ghaGan1", "/path/to/ghaGan1_scaffolds.links") # tab-delimiter
    ass._chain_explotion(adjacency) #adjacency : returned from parse_adjacency("/path/to/adjacency/file")
    #until N50 < 100M
    ass._getSeq("/path/to/ghaGan1.fa")
    ass._addSeq()
    for scf in  ass.scaffolds:
        scf._to_sequence()
    ass._output_generate(filename="/path/to/output/ghaGan1.release.fasta")
    """
    def __init__(self, name, scaffolds = []):
        self.scaffolds = scaffolds
        self.name = name
        self.scf_hash = {scf.name:self.scaffolds.index(scf) for scf in self.scaffolds}
        self.cnts_hash = self._hash_cnts()
        self.seqs = dict()

    @staticmethod
    def with_links(name, links):
        return Assembly(name, scaffolds = parse_links(links))

    def _hash_cnts(self):
        ret = dict()
        for scf in self.scaffolds:
            for cnt in scf.contigs:
                ret[cnt.uname] = scf.name
        return ret

    def _split(self, cntname):
        """
        split scaffold to 2 part
        1: start to cnt
        2: cnt+1 to end
        then update and re-generate hash_cnts
        """
        scfidx = self.scf_hash[self.cnts_hash[cntname]]
        scf = self.scaffolds[scfidx]
        new_scfs = scf._break(cntname)
        self._update(old=[scf], new=new_scfs)
        self.cnts_hash = self._hash_cnts()
        pass

    def _merge(self, cnt1name, cnt2name):
        """
        merge two contigs cnt1 and cnt2:
        - break scaffolds of cnt1 and cnt2 to 2 fragments, return 4 fragments
        - join frag1 of scf1 with cnt2
        - update scaffolds, and re-generate hash_cnts
        """
        scf1idx = self.scf_hash[self.cnts_hash[cnt1name]]
        scf2idx = self.scf_hash[self.cnts_hash[cnt2name]]
        scf1, scf2 = self.scaffolds[scf1idx], self.scaffolds[scf2idx]
        try:
            cnt1, cnt2 = scf1.contigs[scf1.hash_cnts[cnt1name]] , scf2.contigs[scf2.hash_cnts[cnt2name]]
        except KeyError:
            print "Can't merge"
            return
        pair1, pair2 = scf1._break(cnt1name), scf2._break(cnt2name, trimmed=True)
        if (pair1 and pair2) and (cnt1.sign==1 and cnt2.sign==1):
            pair1[0]._join(cnt2)
        self._update(old=[scf1, scf2], new=pair1+pair2)
        self.cnts_hash = self._hash_cnts()
        pass

    def _update(self, old=[], new=[]):
        map(lambda x: self.scaffolds.pop(x), [self.scf_hash[scf.name] for scf in old])
        self.scaffolds+=new
        self.scf_hash = {scf.name:self.scaffolds.index(scf) for scf in self.scaffolds}
        return

    def _getSeq(self, fileSeq):
        for seq in SeqIO.parse(fileSeq, format="fasta"):
            self.seqs[seq.id] = seq.seq
        pass

    def _addSeq(self):
        for scf in self.scaffolds:
            scf._add_seq(self.seqs)
        pass

    def _validate(self):
        scf_lens = list()
        for scf in self.scaffolds:
            scf_len = sum([cnt.end - cnt.start + cnt.link for cnt in scf.contigs])
            scf_lens.append(scf_len)
        n50 = ul._calc_n50(scf_lens, sum(scf_lens))
        return n50

    def _chain_explotion(self, adj):
        for key, value in adj.items():
            if self._check_consistent(key, value):
                self._split(key)
            if self._validate() < 100000000:
                break
        pass

    def _check_consistent(self, cnt1, cnt2):
        """
        Only allow 2 contigs in regular orientation, 1 is forward and the other is reverse_complement
        and 2 contigs came from 2 different scaffolds
        """
        try:
            scfidx1 = self.cnts_hash[cnt1]
            scfidx2 = self.cnts_hash[cnt2]
        except KeyError:
            return False
        if scfidx1 != scfidx2:
            return True # tricky
        scf = self.scaffolds[self.scf_hash[scfidx1]]
        cntidx1, cntidx2 = scf.hash_cnts[cnt1], scf.hash_cnts[cnt1]
        if cntidx2 - cntidx1 == 1:
            return False

    def _output_generate(self, filename):
        assembly_fasta = dict()
        for scf in self.scaffolds:
            assembly_fasta[scf.name] = scf.seq
        ul.write_fasta_dict(fasta_dict=assembly_fasta, filename=filename)
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
                _ = ""
            if "[" in name:
                raw_region = name[name.index("[")+1:name.index("]")]
                region = tuple(map(int, raw_region.split(":")))
                contigs.append(Contig(uname=name[1:], link = int(gap),
                                      sign=name[0], region=region, supported = _ ))
            else:
                contigs.append(Contig(uname=name[1:], region = (0,int(end)),
                              link = int(gap), sign=name[0], supported = _))
        assembly.append(Scaffold(name = arr[0], contigs = contigs))
    return assembly

def parse_adjacency(filename):
    adj = {}
    for line in open(filename):
        left, right = line.strip().split("\t")
        adj[left] = right
    return adj
