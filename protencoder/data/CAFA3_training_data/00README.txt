CAFA3 training data

The files in this package can be used to train a protein function predictor. 

Content: 
1. uniprot_sprot_exp.fasta
2. uniprot_sprot_exp.txt

The .fasta file was downloaded from the UniProt database in September 2016; only Swiss-Prot sequences are included (in the fasta format). The file contains protein sequences of all experimentally annotated proteins. The experimental annotations (class labels) for these proteins are available in the .txt file in the following format:

<accession number> <go term> <namespace>

where <namespace> indicates one of the three GO ontologies:

F: molecular function
P: biological process
C: cellular component

=====
=====

Useful links:

1. All unlabeled protein sequences from the UniProt database can be downloaded from http://www.uniprot.org/. The Swiss-Prot release can be obtained from the UniProt web site.


2. Other sequences and experimental annotations (class labels) can be downloaded from the UniProt-GOA database, available at http://www.ebi.ac.uk/GOA. Merging UniProt-GOA annotations with Swiss-Prot may require conversion of sequence IDs and accession numbers. The Swiss-Port database is generally synchronized with UniProt-GOA.


3. Another set of sequences and annotations can be retrieved directly from the Gene Ontology web site, available at http://geneontology.org/. This web site is maintained by the Gene Ontology Consortium. The current versions of the Gene Ontology can be downloaded from http://geneontology.org/ontology/go.obo whereas the consortium also provides their own release of functionally annotated proteins through their ftp site. Using this resource may require additional ID conversion between these databases. Furthermore, identical proteins might be annotated differently. This requires further reconciliation. 


4. Other resources include the EcoCyc database (http://ecocyc.org/), FlyBase (http://flybase.org/), the Saccharomyces Genome Database (http://www.yeastgenome.org/), the Candida Genome Database (http://www.candidagenome.org/), etc.


5. CAFA1 and CAFA2 publications are available from

http://www.ncbi.nlm.nih.gov/pubmed/23353650

http://www.ncbi.nlm.nih.gov/pubmed/27604469

A gentle introduction to protein function prediction is available from

http://www.cs.indiana.edu/~predrag/papers/radivojac_cafa_2013.pdf


6. For questions and comments regarding this dataset, contact Yuxiang Jiang and Predrag Radivojac at {yuxjiang, predrag}@indiana.edu

======

The CAFA Organizers
September 24, 2016
Bloomington, Indiana 47405
United States of America
