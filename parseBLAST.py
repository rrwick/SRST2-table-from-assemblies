#!/usr/bin/env python

"""
This is a module providing classes to process BLAST outputs. It works for the script detector.py.

Python versions 2.7 and 3 compatible

Copyright (C) 2017 Yu Wan <wanyuac@gmail.com, https://github.com/wanyuac>
Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
First and the latest edition: 6-9 Sep 2017
"""

from collections import namedtuple


"""
The following named tuple class stores information extracted from each row of the BLAST outputs.
The allele name "allele" is extracted from the query name and is modified according to the logical variable unique_allele_symbols in detector.py.
"""
Hit = namedtuple("Hit", ["query", "cluster", "allele", "query_length", "coverage", "contig", "hit_length", "start", "end", \
                         "identity", "bit_score", "perfect_match", "hit_seq"])

""" Store each haplotype (represented by a hit with a unique consensus sequence) and its copy number in each assembly """
Hap = namedtuple("Hap", ["rep_hit", "copy_num"])  # rep_hit: a representative hit for this haplotype; copy_num: copy number of this haplotype in an assembly


class Cluster:
    """ Manages BLAST hits under the same sequence cluster (also known as a gene). """
    
    def __init__(self, name, hits = []):
        """ Initiate private properties """
        self.__name = name  # cluster (gene) name
        self.__hits = hits  # a list of Hit objects. Hits may be an empty list when a cluster object is initiated.
        self.__size = len(hits)  # number of hits within this cluster
        return
    
    
    @property
    def name(self):
        return self.__name  # Users are not allowed to change the name attribute
    
    
    @property
    def size(self):  # get the number of hits
        return self.__size
        
        
    def add_hit(self, hit):
        self.__hits.append(hit)  # add a hit in the hits list
        self.__size = len(self.__hits)  # refresh the hit count
        return
        
        
    def filter_hits(self, min_coverage = 90.0, max_divergence = 10.0):
        if self.__size > 0:  # No action is taken for an empty cluster.
            min_identity = 100.0 - max_divergence
            hits_kept = [h for h in self.__hits if h.identity >= min_identity and h.coverage >= min_coverage]  # use a list comprehension to simplify the expression
            self.__hits = hits_kept  # Notice hits_kept may be an empty list.
            self.__size = len(self.__hits)  # refresh the hit count
        return
        
        
    def find_best_hit(self, hits = None):
        """ Determine the best hit within this cluster """
        if hits is None:
            hits = self.__hits
        if len(hits) > 0:
            best_hit = hits[0]  # arbitarily treat the first hit as the best hit of the current cluster                
            for hit in hits:
                is_hit_perfect = (hit.identity == 100.0 and hit.coverage == 100.0)  # for the current hit
                is_current_best_perfect = (best_hit.identity == 100.0 and best_hit.coverage == 100.0)  # for the current best hit, which may not be the current hit
                
                # If this hit is perfect and the current best isn"t perfect, then this is the new best hit.
                if is_hit_perfect and (not is_current_best_perfect):
                    best_hit = hit

                # If this hit is perfect and the current best is also perfect, then this is the new best hit only if it is longer.
                elif is_hit_perfect and is_current_best_perfect and (hit.hit_length > best_hit.hit_length):
                    best_hit = hit

                # If neither this hit nor the current best are perfect, this is the new best only if it has a higher bit score.
                elif (not is_hit_perfect) and (not is_current_best_perfect) and (hit.bit_score > best_hit.bit_score):
                    best_hit = hit
        else:  # when some clusters lose all of their hits after filtering for the query coverage and nucleotide identity
            best_hit = None
        return best_hit  # a single Hit object
    
    
    def find_all_copies(self, max_overlapping_nt = 0):
        """
        Search for all physically separated alleles of the same cluster and return the result as a dictionary {contig : hits}.
        max_overlapping_nt: the maximal number of nucleotides allowed to share between to hits. Hits overlapping by more than
        this threshold will be considered as hits to the same allele.
        """
        if self.__size > 0:
            hit = self.__hits[0]  # Arbitrarily take the first hit to initiate a list of Hit objects.
            copies = {hit.contig : [hit]}  # self.__hits may be filtered for the coverage and nucleotide identity.
            if self.__size > 1:  # when there are other hits within this cluster
                for hit in self.__hits[1 : ]:  # go through the rest of hits
                    overlaps = self.__find_overlaps(hit, copies, max_overlapping_nt)
                    if overlaps is None:  # when there is no overlap at all
                        if hit.contig in list(copies.keys()):
                            copies[hit.contig].append(hit)
                        else:
                            copies[hit.contig] = [hit]
                    else:  # take the better/best hit when overlaps happen.
                        separate_hits = copies[overlaps.contig]  # Hits in copies[contig] must be physically separate, but they show overlaps with the query hit.
                        overlapping_hits = [separate_hits[i] for i in overlaps.hit_indices]  # a list comprehension extracting items from list A given a list of indices B
                        overlapping_hits.append(hit)  # Now the list of overlapping hits is complete.
                        best_hit = self.find_best_hit(overlapping_hits)
                        new_separate_hits = [hit for index, hit in enumerate(separate_hits) if index not in overlaps.hit_indices]  # https://stackoverflow.com/questions/497426/deleting-multiple-elements-from-a-list
                        new_separate_hits.append(best_hit)  # best_hit may replace two or more existing hits in copies[contig], although it is unlikely to happen.
                        copies[overlaps.contig] = new_separate_hits
        else:  # for an empty cluster
            copies = {}
        return copies
    
    
    def __find_overlaps(self, query, hits_dict, max_overlapping_nt = 0):
        """
        A private function determines whether a query hit (query) overlaps any one in a dictionary of hits (hits_dict).
        It returns a named tuple for the contig where the query overlaps with existing hits and the indices of
        these hits. This function works for the self.find_all_copies method. hits_dict: a dictionary of hits, namely,
        {contig : hit}. hit_index: the index of a hit on a given contig.
        """
        Overlaps = namedtuple("Overlaps", ["contig", "hit_indices"])
        overlaps = None  # the default return value when there are no overlapping hits at all
        for contig, hits in hits_dict.items():
            if query.contig == contig:  # when both hits come from the same contig
                i = 0  # The index starts from zero.
                indices = []  # where hits overlap
                for hit in hits:
                    coords = [query.start, query.end, hit.start, hit.end]
                    span = max(coords) - min(coords) + 1  # the minimal length (number of bases) covering both hits
                    overlapped_len = query.hit_length + hit.hit_length - span  # becomes negative when both hits are overlapping
                    if overlapped_len > max_overlapping_nt:
                        indices.append(i)  # find an overlap
                    i += 1  # move to the next item
                if len(indices) > 0:  # when the query overlaps hits on this contig
                    overlaps = Overlaps(contig = contig, hit_indices = indices)
                    break  # Break the for loop because a query can only show overlaps to hits on a single contig.
        return overlaps  # either an Overlaps object or a None value
    

class Assembly:
    """ Manages BLAST outputs for a given FASTA file """
    
    def __init__(self, name, clusters = {}, excl_empty = True):
        """
        Initiate private properties
        excl.empty: exclude empty clusters (ie. cluster.__size = 0)        
        """
        self.__name = name  # assembly name
        self.__clusters = clusters  # a dictionary of the Cluster objects
        if excl_empty:
            self.__clusters = {cluster_name : cluster for cluster_name, cluster in self.__clusters.items() if cluster.size > 0}
        self.__size = len(self.__clusters)  # number of valid clusters
        return
    
    
    @property
    def name(self):
        return self.__name
    
    
    @property
    def cluster_num(self):
        return self.__size
    
    
    @property
    def cluster_names(self):
        return list(self.__clusters.keys())
        
        
    def add_cluster(self, cluster):
        self.__clusters.update(cluster)  # In fact, this method is able to append multiple clusters at once.
        self.__size = len(self.__clusters)
        return
    
    
    def find_best_hits(self):
        best_hits = {}
        if self.__size > 0:
            for cluster_name, cluster in self.__clusters.items():
                bh = cluster.find_best_hit()  # bh: the best hit
                if bh is not None:  # bh may equal None when excl_empty = False.
                    best_hits[cluster_name] = bh
        return best_hits
    
    
    def find_all_copies(self, max_overlapping_nt = 0):
        """ Search for all putative allele calls for each cluster in the current assembly """
        copies = {}
        if self.__size > 0:
            for cluster_name, cluster in self.__clusters.items():
                hits = []
                copies_cluster = cluster.find_all_copies(max_overlapping_nt)  # return a dictionary {contig : [hit1, hit2, ...]}
                for contig in list(copies_cluster.keys()):
                    hits += copies_cluster[contig]
                copies[cluster_name] = hits  # We do not care contig information here.
        return copies
