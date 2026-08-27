# PDC000614 manual input

Do **not** commit the source payloads in this directory.

The from-zero launcher imports the PDC000614 files from this inbox into the project-local data tree. Required source semantics are:

- one TMT18 proteome table matching `*tmt18.tsv`
- one sample-map text file matching `*sample.txt`
- one biospecimen TSV matching `*biospecimen*.tsv`

Optional contextual files:

- `*label.txt`
- `*summary.tsv`

The current analysis uses same-patient, same-plex tumour-normal pairing. Obtain the source files from the official PDC PDC000614 study page and retain the original filenames/checksums where possible.
