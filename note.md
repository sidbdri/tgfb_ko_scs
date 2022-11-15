- The data contains cells from mouse whole brain sample(mix cell type).
- 4WT + 4 KO, induced KO of an exon(ENSE00003688775) of the TGFB2 gene. The KO is NOT induced in astrocyte cells.
- After cellranger run, we notice that there is not much difference between the genotype
- We example the KO reads
  - We found that there are only a few number of reads falling into these reagion. The gene is not express much.
  - With the low number, and the astrocyte expression, it is hard to tell if KOs are KO, but the WTs do seem to be WT
  - Jing wanted to check where these reads are from
    - I found out the cells expressing those reads in KOs and they look like they are from  astrocyte
  - Jing went back to qrtpcr to comfirm genotype, and give me a updated set of CMO/genotype mapping, I rerun cellranger aggr. 
  - Jing inspect the results and think it is worse than before. 
  - Jing decided to go back for more qrtpcr check for the genotype.
- Jing explain that
  - the cell type population is not as expect. She is expecting way less microglia and way more neurons.
  - something might have go wrong during the lib prep/sequencing
  - the cell compisition might be another reason why the genotype difference is subtle

Jing decided, for now, we are going to stop looking at this data due to the cell population is very wrong.

