# Implementation provenance

The source snapshots document the numerical noise and extragalactic SED
conventions. They are reference records and are not imported by the release
drivers. The executable code has no dependencies on the original private
utility modules.

The current default G frame for missing-header warehouse inputs follows the
dataset owner declaration of 11 October 2026. Historical provider/header
tracing establishes C/ICRS for upstream mmDL products; that evidence remains
in the manuscript package. The default G assignment is not independently
verified physical-frame metadata. Input FITS headers and explicit source-frame
configuration remain available for correct alignment.

Independent software fixture tests, actual serial/MPI noise controls and
native2048 integration are recorded separately. The source snapshots cannot
recover unrecorded seeds for the inherited archive.
