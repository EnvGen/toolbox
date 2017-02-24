# Assign taxonomy from blast results


## How it works
Each blast hit for a gene gets a weight:

    fHit = float(alignment_length_in_query)/query_length
    fHit *= float(percent_identity)/100.0
    fHit = min(1.0,fHit)


Which is used to put a weight each taxa at level ```depth``` get a weight

    weight_for_taxa = (fHit - MIN_IDENTITY_TAXA[depth])/(1.0 - MIN_IDENTITY_TAXA[depth])


and MIN_IDENTITY_TAXA[depth] is the minimum level of identity we accept as a proof for a gene to be placed on that level. Only weights higer than 0.0 are considered.

## Example with one single hit
As an example if a query gene has a single hit for a gene in the database belonging to the genus Escherichia with a value for fHit=0.75:

    weight_for_Escherichia = (0.75 - MIN_IDENTITY_TAXA[genus])/(1.0 - MIN_IDENTITY_TAXA[genus])
    weight_for_Escherichia = (0.75 - 0.90])/(1.0 - 0.90) < 0.0

but
    weight_for_Enterobacteriales = (0.75 - MIN_IDENTITY_TAXA[order])/(1.0 - MIN_IDENTITY_TAXA[order])
    weight_for_Enterobacteriales = (0.75 - 0.70])/(1.0 - 0.70) = 0.1667

Since there are no other hits for this query gene, the normalized weight for Enterobacteriales will be 1.0 and the gene will be assigned to the order Enterobacteriales. However no assignment will be made for any more specific level.


## Example with three hits
Let's extend the previous example with two more hits in the database. One hit to a Citrobacter amalonaticus (Bacteria, Proteobacteria, Gammaproteobacteria, Enterobacteriales, Enterobacteriaceae) with fHit=0.79 and one hit to Pseudomonas aeruginosa ((Bacteria, Proteobacteria, Gammaproteobacteria, Pseudomonadales, Pseudomonadaceae) with fHit=0.
