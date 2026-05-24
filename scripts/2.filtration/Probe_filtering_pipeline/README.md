# Probe filtering pipeline

This folder contains the probe filtering workflow used to remove candidate padlock probes with potentially undesirable sequence properties, including strong secondary structures and strong probe-probe interactions.

The filtering pipeline is organized into two main checks:

1. `2nd_structure/`  
   Filters probes based on predicted secondary structures.

2. `probe_interaction/`  
   Filters probes based on predicted interactions between probes.

This pipeline can be used either for general probe pool filtration or, when needed, for selecting a final probe set in which probes targeting the same gene carry distinct additional barcodes.

## Recommended workflow for probes with unique additional barcodes

For applications where users want 4 probes per gene and the final selected probes for each gene to be distinguishable by the additional barcode, we recommend the following two-pass design and filtering strategy.

### Overview

```text
Design excess candidate probes (e.g. 8 probes per gene)
        |
        v
Assign additional barcodes to candidate probes
        |
        v
Run probe filtering (2nd structure & probe interaction)
        |
        v
Select the desired number of probes per gene (up to 4 probes)
        |
        v
Re-assign or edit additional barcodes so selected probes are unique within each gene
        |
        v
Run probe filtering again (2nd structure & probe interaction)
        |
        v
Finalize probe pool
```

### Step-by-step procedure

1. Design more candidate probes than the final number needed

- For each gene, first design more candidate probes than will ultimately be used.

- For example, if the final experiment requires 4 probes per gene, we recommend designing more than 4 candidate probes per gene, such as 6-8 or more if possible. This provides backup candidates in case some probes are removed during secondary structure or probe-probe interaction filtering.

2. Assign additional barcodes to the candidate probes

- Assign additional barcodes to candidate probes before the first filtering step.

- The purpose of the additional barcode is to distribute probes targeting the same gene across different color channels in the additional sequencing round. This can help reduce optical crowding when signals are locally clustered or when a gene is highly expressed.

- For a 4-color additional barcode round, there are 4 possible additional barcode identities. If the number of probes per gene is greater than 4, some candidate probes for the same gene will necessarily share the same additional barcode at this stage.

3. Run the probe filtering pipeline

- Run the secondary-structure and probe-probe interaction filters on the candidate probe pool.

- This step removes probes with potentially strong secondary structures or strong interactions with other probes in the pool.

- After filtering, some genes may have fewer remaining candidate probes than originally designed.

4. Select the desired number of probes per gene

- After the first filtering step, select the desired number of probes per gene from the surviving candidates.

- For example, if the target design is 4 probes per gene, select 4 filtered probes for each gene whenever possible.

- At this stage, check whether the selected probes for each gene have unique additional barcodes.

5. Re-assign additional barcodes to make them unique within each gene

- If the selected probes for a gene contain repeated additional barcodes, edit the additional barcode assignments so that the final selected probes have unique additional barcodes whenever possible.

- For example, for a 4-probe-per-gene design using a 4-color additional barcode round, assign the 4 selected probes for each gene to 4 distinct additional barcode identities.

- This step is important because probe filtering may remove some probes unevenly. Even if the original candidate pool was designed with balanced additional barcode usage, the final filtered subset may still contain repeated additional barcodes for the same gene.

6. Run the probe filtering pipeline again

- After changing the additional barcode sequences, run the secondary-structure and probe-probe interaction filtering steps again.

- This second filtering step is necessary because changing the additional barcode changes the full padlock probe sequence. The edited probes may therefore have different secondary-structure or probe-probe interaction properties from the original filtered probes.

7. Final quality-control checks

* Before finalizing the probe pool, we recommend checking the following:

  * Each gene has the intended number of probes whenever possible.
  * For each gene, the selected probes have unique additional barcodes whenever possible.
  * The final edited probe sequences pass secondary-structure filtering.
  * The final edited probe sequences pass probe-probe interaction filtering.
  * The final probe table contains the expected gene names, target sequences, padlock sequences, additional barcodes, and decoding information.
