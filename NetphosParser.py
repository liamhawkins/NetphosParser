from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pandas as pd
import fnmatch

class NetphosParser:
    def __init__(self):
        self.records = []
        
    def __len__(self):
        return len(self.records)
    
    def parse(self, file):
        '''Parse gff results from NetPhos 3.1'''
        with open(file) as f:
            content = f.readlines()
            content = [x.rstrip() for x in content]
        record_starts = []
        for index,line in enumerate(content):
            if line.startswith('##Type Protein'):
                record_starts.append(index)

        for i in range(len(record_starts)):
            if i+1 < len(record_starts):
                self.records.append(self.Record(content[record_starts[i]:record_starts[i+1]]))
            else:
                self.records.append(self.Record(content[record_starts[i]:]))
                
    def filter_(self, df):
        '''Remove hits from each record that do not correspond to peptides in df'''
        def nth_repl(s, sub, repl, nth):
            # From Padraic Cunningham https://stackoverflow.com/a/35092436/7292937
            find = s.find(sub)
            # if find is not p1 we have found at least one match for the substring
            i = find != -1
            # loop util we find the nth or we find no match
            while find != -1 and i != nth:
                # find + 1 means we start at the last match start index + 1
                find = s.find(sub, find + 1)
                i += 1
            # if i  is equal to nth we found nth matches so replace
            if i == nth:
                return s[:find]+repl+s[find + len(sub):]
            return s
        
        def matches_peptide(hit, peptides, orig_peptides):
            pad_peptides = []
            hit.matched_peptides = []
            for index,peptide in enumerate(peptides):
                unmod_peptide = peptide
                pept_residue_index = peptide.index('#')-1
                residue = peptide[pept_residue_index]
                if residue != hit.residue:
                    continue
                peptide = peptide.replace("#","")
                diff = hit.residue_index-pept_residue_index
                if diff > 0:
                    pad_peptide = '?'*diff+peptide[:len(hit.context)-diff]
                    pad_peptide = pad_peptide+'?'*(len(hit.context)-len(pad_peptide))
                elif diff < 0:
                    pad_peptide = peptide[-diff:len(hit.context)-diff]
                    pad_peptide = pad_peptide+'?'*(len(hit.context)-len(pad_peptide))
                else:
                    pad_peptide = peptide[:len(hit.context)]
                    pad_peptide = pad_peptide+'?'*(len(hit.context)-len(pad_peptide))
                
                if fnmatch.fnmatch(pad_peptide, hit.context):
                    hit.matched_peptides.append(unmod_peptide)
                    hit.orig_peptide = orig_peptides[peptides.index(unmod_peptide)]
                
                pad_peptides.append(pad_peptide)
            
            return len(fnmatch.filter(pad_peptides, hit.context)) > 0
            
        
        # Split each peptide with >1 phos-site into multiple peptides each containing 1
        new_rows = []
        for index, row in df.iterrows():
            if row['peptide'].count('#') > 1:
                for i in range(row['peptide'].count('#')):
                    new_row = row.copy(deep=True)
                    new_row['peptide'] = nth_repl(row['peptide'], '#', '', i+1)
                    new_rows.append(new_row)
                df = df.drop([index])
        df = df.append(new_rows)
                
        # Delete records that aren't in df
        del_records = []
        for index, record in enumerate(self.records):
            if not record.id in df['uniprot'].values:
                del_records.append(index)
        for i in sorted(del_records, reverse=True):
            del self.records[i]
        
        del_records = []
        for index, record in enumerate(self.records):
            matching_peptides = df[df['uniprot'] == record.id]['peptide'].values.tolist()
            orig_peptides = df[df['uniprot'] == record.id]['orig_peptide'].values.tolist()
            del_hits = []
            
            for h_index, hit in enumerate(record.hits):
                if not matches_peptide(hit, matching_peptides, orig_peptides):
                    del_hits.append(h_index)
            for i in sorted(del_hits, reverse=True):
                del record.hits[i]
            if len(record.hits) == 0:
                del_records.append(index)
        for i in sorted(del_records, reverse=True):
            del self.records[i]
            
    def to_df(self):
        df = pd.DataFrame(columns=['uniprot', 'orig_peptide', 'residue', 'residue_index', 'context', 'kinase', 'matching_peptides', 'score'])
        i=0
        for record in self.records:
            for hit in record.hits:
                df.loc[i] = [record.id, hit.orig_peptide, hit.residue, hit.site_index, hit.context, hit.kinase, '; '.join(hit.matched_peptides), hit.score]
                i+=1
        return df
                
                
    class Record:
        def __init__(self, inp):
            assert inp[0].startswith('##Type Protein')
            assert inp[1].startswith('##Protein')
            self.id = inp[0].split(' ')[2]
            self.hits = []
            self.get_sequence(inp)
            self.get_hits(inp)
        
        def __len__(self):
            return len(self.hits)
        
        def get_sequence(self, inp):
            start = 2
            end = inp.index('##end-Protein')
            seq = ''.join([x.replace('#', '') for x in inp[start:end]])
            self.sequence = SeqRecord(Seq(seq,IUPAC.protein),id=self.id)
        
        def get_hits(self, inp):
            start = inp.index('##end-Protein') + 3
            hit_list = inp[start:]
            for hit in hit_list:
                self.hits.append(self.Hit(self, hit))
            
        class Hit:
            def __init__(self, parent, inp):
                self.parent = parent
                inp_list = inp.split()
                self.seqname = inp_list[0]
                self.source = inp_list[1]
                self.kinase = inp_list[2]
                self.start = int(inp_list[3])
                self.score = float(inp_list[5])
                self.above_threshold = inp_list[-1] == 'YES'
                self.residue = str(self.parent.sequence.seq[self.start-1])
                self.get_context()
                
            def __str__(self):
                return str(self.__dict__)
            
            def get_context(self):
                self.num_flank = 7 # 7+1+7=15aa context
                self.site_index = self.start - 1
                if self.site_index - self.num_flank < 0:
                    self.start_index = 0
                    self.residue_index = self.num_flank + (self.site_index - self.num_flank)
                else:
                    self.start_index = self.site_index - self.num_flank
                    self.residue_index = self.num_flank

                if self.site_index + self.num_flank > len(self.parent.sequence.seq):
                    self.stop_index = len(self.parent.sequence.seq)
                else:
                    self.stop_index = self.site_index + self.num_flank
                self.context = str(self.parent.sequence.seq[self.start_index:self.stop_index+1])