#!/bin/bash

# use samtools view and target list bed file to generate sam file for reads overlapping targets
python errorpatterns.py sample.sam sample.fasta

# filter errorpatterns output using germline mpileup (around 15 sites to filter) with grep -Pv
grep -Pv '^chr3\t10188236|^chr4\t1807894|^chr4\t55141055|^chr5\t112175770|^chr5\t170837514|^chr5\t170837515|^chr7\t55249063|^chr9\t5073768|^chr9\t5073769|^chr9\t5073770|^chr9\t5073771|^chr9\t5073772|^chr9\t133748264|^chr10\t123274765|^chr17\t7579373' sample.sam.errorpatterns > sample.sam.errorpatternsNOGERM

# generate table of covariate cell counts as python pickle
python buildcovariatetable.py errorpatternsNOGERM

# analyze pileups
python analyzepileup.py sample.bam sample.fasta sample.sam.pik3ca339432.sam sample.sam.pik3ca339432.bed sample.sam.errorpatternsNOGERM.covariatetable

# not yet sure what to do with the output from above, since we only have Dirichlet MLE code in matlab
