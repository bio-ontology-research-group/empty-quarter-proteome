#!/bin/bash
# Login-node file handling only (tar, no compression of large sets): archive models under archives/.
# usage: archive.sh monomer|complex   ; small laptop bundle: archive.sh bundle LIST_OF_FILES
set -euo pipefail
D=/mnt/data/u/hohndor/eq-exported-dark-20261111
cd $D; mkdir -p archives
case "$1" in
  monomer)
    nice tar -cf archives/monomer_models_all.tar --exclude='*.png' predictions/monomer predictions/monomer_tail predictions/monomer_tail_shallow predictions/monomer_tail_deep
    nice tar -cf archives/monomer_foldseek_queries_and_hits.tar foldseek/queries foldseek/hits results ;;
  complex)
    nice tar -cf archives/complex_real_models_all.tar --exclude='*.png' predictions/complex_real_shallow predictions/complex_real_deep
    nice tar -cf archives/complex_control_models_all.tar --exclude='*.png' predictions/complex_control_shallow predictions/complex_control_deep ;;
  bundle)
    nice tar -czf archives/laptop_bundle_qualified_monomers_and_above_control_complexes.tar.gz -T "$2" ;;
esac
ls -la archives; sha256sum archives/* > archives/SHA256SUMS
