# Output structure

The repository does not contain generated outputs. During a local run the workflow writes project-local artifacts including:

```text
output/
├── objects/
├── tables/
└── figures/
logs/
recovery/
└── audit/
```

Part9 consumes the frozen upstream result tables and generates publication-oriented main figures, supplementary figures, supplementary tables, exact figure-source tables, and audit files. Generated outputs are ignored by Git because they can be large and are reproducible from the documented source inputs.

For manuscript submission, supplementary files may be deposited separately with the journal and/or attached to a tagged GitHub/Zenodo release as release assets.
