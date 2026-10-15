# Additional root finding after round-2 review

The original minimap2 command omitted `-c`; default PAF reports approximate
mapping identities, not full base-level alignment identities. Opus verified
coordinate extraction but did not identify this distinction. The original
results remain preserved under `culture_nucleotide_submitted/`. Root is rerunning
with `-c` to obtain aligned identities and CIGAR strings before promoting the
nucleotide fractions or joining to Bakta gene coordinates. The final revision
must use the new output, not silently present approximate identities as exact.
