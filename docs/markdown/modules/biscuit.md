---
title: BISCUIT
displayed_sidebar: multiqcSidebar
description: >
  Maps bisulfite converted DNA sequence reads and determines cytosine methylation states.
---

<!--
~~~~~ DO NOT EDIT ~~~~~
This file is autogenerated from the MultiQC module python docstring.
Do not edit the markdown, it will be overwritten.

File path for the source of this content: multiqc/modules/biscuit/biscuit.py
~~~~~~~~~~~~~~~~~~~~~~~
-->

:::note
Maps bisulfite converted DNA sequence reads and determines cytosine methylation states.

[https://github.com/huishenlab/biscuit](https://github.com/huishenlab/biscuit)
:::

The module parses logs generated by BISCUIT and the quality control script, QC.sh, included with
the BISCUIT software.

**Note**: As of MultiQC v1.9, the module supports only BISCUIT version v0.3.16 and onwards.
If you have BISCUIT data from before this, please use MultiQC v1.8.

#### Insert Size Distribution

The second tab of this plot uses the config option `read_count_multiplier`,
so if millions of reads is not useful for your data you can customise this.

See [Number base (multiplier)](https://docs.seqera.io/multiqc/#number-base-multiplier)
in the documentation.

### File search patterns

```yaml
biscuit/align_isize:
  contents: BISCUITqc Insert Size Table
  fn: "*_isize_table.txt"
  num_lines: 3
biscuit/align_mapq:
  contents: BISCUITqc Mapping Quality Table
  fn: "*_mapq_table.txt"
  num_lines: 3
biscuit/align_strand:
  contents: BISCUITqc Strand Table
  fn: "*_strand_table.txt"
  num_lines: 3
biscuit/base_avg_retention_rate:
  fn: "*_totalBaseConversionRate.txt"
biscuit/covdist_all_base:
  fn: "*_covdist_all_base_table.txt"
biscuit/covdist_all_base_botgc:
  fn: "*_covdist_all_base_botgc_table.txt"
biscuit/covdist_all_base_topgc:
  fn: "*_covdist_all_base_topgc_table.txt"
biscuit/covdist_all_cpg:
  fn: "*_covdist_all_cpg_table.txt"
biscuit/covdist_all_cpg_botgc:
  fn: "*_covdist_all_cpg_botgc_table.txt"
biscuit/covdist_all_cpg_topgc:
  fn: "*_covdist_all_cpg_topgc_table.txt"
biscuit/covdist_q40_base:
  fn: "*_covdist_q40_base_table.txt"
biscuit/covdist_q40_base_botgc:
  fn: "*_covdist_q40_base_botgc_table.txt"
biscuit/covdist_q40_base_topgc:
  fn: "*_covdist_q40_base_topgc_table.txt"
biscuit/covdist_q40_cpg:
  fn: "*_covdist_q40_cpg_table.txt"
biscuit/covdist_q40_cpg_botgc:
  fn: "*_covdist_q40_cpg_botgc_table.txt"
biscuit/covdist_q40_cpg_topgc:
  fn: "*_covdist_q40_cpg_topgc_table.txt"
biscuit/cpg_retention_readpos:
  fn: "*_CpGRetentionByReadPos.txt"
biscuit/cph_retention_readpos:
  fn: "*_CpHRetentionByReadPos.txt"
biscuit/dup_report:
  contents: BISCUITqc Read Duplication Table
  fn: "*_dup_report.txt"
  num_lines: 3
biscuit/qc_cv:
  contents: BISCUITqc Uniformity Table
  fn: "*_cv_table.txt"
  num_lines: 3
biscuit/read_avg_retention_rate:
  fn: "*_totalReadConversionRate.txt"
```
