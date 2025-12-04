<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-ampliseq_logo_dark.png">
    <img alt="nf-core/ampliseq" src="docs/images/nf-core-ampliseq_logo_light.png">
  </picture>
</h1>

ArkeaBio's implementation of [nf-core/ampliseq](https://github.com/nf-core/ampliseq).

# Changes from nf-core
1. [Parameter](arkea-params.json) and [configuration](arkea.config) files included.
2. Utility script for creating a samplesheet
3. Test data
4. Beta diversity distance matrices exported to results

# Development
This version of ampliseq lives permanently on sandbox at `/data/ampliseq/ampliseq-updated`. Updates to this pipeline must be pulled there:

```bash
ssh sandbox
cd /data/ampliseq/ampliseq-updated
git pull
chmod -R +x .
```

A test dataset is included in this repo. You can run the pipeline on the test data with

```bash
nextflow main.nf \
  -profile docker \
  --input arkea-test-data/test_samplesheet.csv \
  --metadata arkea-test-data/test_metadata.tsv \
  -params-file arkea-params.json \
  -c arkea.config
```