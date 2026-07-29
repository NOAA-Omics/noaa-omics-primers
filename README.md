# NOAA 'Omics Metabarcoding Assays

Metabarcoding assays used for DNA metabarcoding by NOAA 'Omics researchers, using terms from the [FAIRe](https://fair-edna.github.io) standard. This repository provides assay information used by the [Ocean DNA Explorer](https://oceandnaexplorer.org).

## Assays

An assay is defined by the unique combination of forward and reverse PCR primers, which have fixed data associated with them (taxonomic target, gene, subfragment, amplicon size, primer sequences, primer names, and primer references).

The file [assays.tsv](https://github.com/NOAA-Omics/noaa-omics-metabarcoding-assays/blob/main/assays.tsv) lists all assays fully documented and approved to work with Ocean DNA Explorer.

Assay field descriptions:

- `assay_name`: A brief, concise identifier for the assay with no spaces or special characters, ensuring machine readability. Must be unique. Suggest including the targeted taxonomic group, gene, subfragment (if applicable), and author or commonly used name of assay. Source: [FAIRe](https://github.com/FAIR-eDNA/FAIRe_checklist).
- `assay_name_alternate`: Other names the assay may be called.
- `targetTaxonomicAssay`: Taxa or species name targeted by the primers. Source: [FAIRe](https://github.com/FAIR-eDNA/FAIRe_checklist).
- `targetTaxonomicScope`: The taxonomic group(s) targeted in the study. This can differ from the targetTaxonomicAssay. For example, the targetTaxonomicAssay may be "Chordata" while the targetTaxonomicScope is "Chondrichthyes" (bony fish, sharks and rays). Source: [FAIRe](https://github.com/FAIR-eDNA/FAIRe_checklist).
- `target`: Target taxa, gene and/or molecule for the protocol. Use terms from [NCBITAXON](https://www.ebi.ac.uk/ols/ontologies/ncbitaxon) or [NCIT](https://www.ebi.ac.uk/ols/ontologies/ncit). Source: [MIOP](https://github.com/BeBOP-OBON/miop).
- `target_gene`: Targeted gene or locus name. Use terms from FAIRe guidelines. Source: [FAIRe](https://github.com/FAIR-eDNA/FAIRe_checklist).
- `target_subfragment`: Name of subfragment of a gene or locus. Used to identify specific regions on marker genes like V6 of 16S rRNA. Use terms from FAIRe guidelines. Source: [FAIRe](https://github.com/FAIR-eDNA/FAIRe_checklist).
- `pcr_primer_[forward|reverse]`: Primer sequence (5'->3'). The primer sequence should NOT contain MIDs and adapter sequences and should be reported in uppercase letters. In primers containing an "N" (any base), the published sequence may contain an "I" (inosine, pairs with any base). Source: [FAIRe](https://github.com/FAIR-eDNA/FAIRe_checklist).
- `pcr_primer_name_[forward|reverse]`: Standardized primer name, using published name if possible but adding the name of the target gene, target subfragment (if applicable), and in some cases author information. Source: [FAIRe](https://github.com/FAIR-eDNA/FAIRe_checklist).
- `pcr_primer_name_published_[forward|reverse]`: Published primer name from original publication.
- `pcr_primer_reference_[forward|reverse]`: Link to original publication (DOI if available) of primer. Source: [FAIRe](https://github.com/FAIR-eDNA/FAIRe_checklist).
- `assay_reference`: Link to original publication and/or protocol (DOI if available) of assay.

## Assay Preps

An assay prep contains additional fields whose values can vary from lab to lab. This metadata should be contained in a protocol available online, with the DOI in the field `nucl_acid_amp`. This DOI should be a Zenodo DOI of a BeBOP-formatted protocol, a protocols.io DOI, or some other stable DOI or URL of the protocol.

The assay prep fields are described here:

- `nucl_acid_amp`: Link to assay protocol. Example: https://doi.org/10.5281/zenodo.15636182
- `ampliconSize`: The length of the amplicon in basepairs *excluding* the primers, adapters, and MIDs. A range can be entered separated by a bar "|" (e.g., "140 | 160"). Units in basepairs.
- `thermocycler`
- `commercial_mm`
- `custom_mm`
- `pcr_cond`
- `barcoding_pcr_appr`
- `pcr2_commercial_mm`
- `pcr2_cond`
- `sequencing_location`
- `platform`
- `instrument`
- `seq_kit`
- `adapter_forward`
- `adapter_reverse`
- `lib_screen`
- `checksum_method`
- `seq_method_additional`
