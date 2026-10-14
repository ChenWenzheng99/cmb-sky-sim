The theory builder verifies primary-source files against `source_manifest.json`.
Retrieve these files from their original providers before rebuilding the theory:

```bash
# From the repository root:
python3 scripts/fetch_assets.py --theory-sources
```

The precomputed spectra in `../outputs/` can be used for remapping without
re-running CAMB. The helper checks exact file sizes and SHA256 values; a provider
change is reported rather than silently accepted as the recorded source.
