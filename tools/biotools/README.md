# bio.tools schema

`biotools.xsd` is the unmodified stable biotoolsSchema from
<https://github.com/bio-tools/biotoolsSchema/blob/c31233af4e136f985e83a58f3ca11b02628348f9/stable/biotools.xsd>.
It is maintained by the bio.tools / ELIXIR community and distributed under
CC BY-SA 4.0; the upstream license is included in `LICENSE`.

The local copy makes generation and validation reproducible and usable offline.
When updating it, record the upstream revision here and run
`Rscript tools/test-biotools.R` before regenerating `inst/biotools.xml`.

The initial annotations were checked against EDAM ontology revision
`3a6e2ca4f5eba34d189e737808909f3c96270c84`:
<https://github.com/edamontology/edamontology/blob/3a6e2ca4f5eba34d189e737808909f3c96270c84/EDAM_dev.owl>.
Schema validation checks EDAM URI syntax, not ontology membership, labels or
deprecation; check those against EDAM when changing annotations.
