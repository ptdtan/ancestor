MIN_GAP_SIZE = 11

import ultilities as ul
from collections import defaultdict
from Bio import SeqIO
from copy import copy

class Contig:
    def __init__(self, uname=None, seq=None,
                 link=MIN_GAP_SIZE, sign="+", region=(0,0)):
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
                sequence.append(''.join([cnt.seq,'N'*cnt.link]))
            else:
                sequence.append(''.join([ul.reverse_complement(cnt.seq), 'N'*cnt.link]))
        self.seq = ''.join(sequence)
        pass
class Assembly:
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
        scfidx = self.scf_hash[self.cnts_hash[cntname]]
        scf = self.scaffolds[scfidx]
        new_scfs = scf._break(cntname)
        self._update(old=[scf], new=new_scfs)
        self.cnts_hash = self._hash_cnts()
        pass

    def _merge(self, cnt1name, cnt2name):
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
            scf_len = sum([len(cnt.seq)+cnt.link for cnt in scf.contigs])
            scf_lens.append(scf_len)
        n50 = ul._calc_n50(scf_lens, sum(scf_lens))
        print "N50: %d" %(n50)
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
                contigs.append(Contig(uname=name[1:], link = int(gap), sign=name[0], region=region))
            else:
                contigs.append(Contig(uname=name[1:], region = (0,end-start+1),
                              link = int(gap), sign=name[0]))
        assembly.append(Scaffold(name = arr[0], contigs = contigs))
    return assembly
