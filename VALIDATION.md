# Distribution validation

The GitHub distribution was checked on 14 October 2026 using fresh installations
of the supplied base runtime and the CMB extension runtime. All 57 existing
tests passed: 32 noise/observation, 3 theory, 13 CMB-production and 9 CMB-validation
tests. The numerical sources remain byte-identical to the validated packages.

Five MPI ranks completed a gather operation with Open MPI 5.0.8. The local host
emitted network-interface pairing warnings during this successful check. The
example configurations, command help, shell syntax, theory contracts, downloaded
asset hashes and all seven CMB wheels were also checked. The bootstrap downloader
accepted a valid local asset and rejected invalid hashes, path traversal and an
invalid wheel archive.

These are software and packaging checks. They did not regenerate the production
maps or repeat a scientific ensemble assessment. Existing small synthetic test
fixtures are not downgraded production maps. The archived 50-realization results
are in [validation_50/](validation_50/README.md); older validation records retain
their original scopes. The GitHub distribution carries source and runtime assets,
while large sky maps remain external and PySM templates are downloaded from the
provider. See [the machine-readable report](provenance/github_distribution_validation.json).
