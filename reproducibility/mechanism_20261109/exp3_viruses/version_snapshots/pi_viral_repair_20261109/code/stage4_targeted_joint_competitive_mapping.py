#!/usr/bin/env python3
"""Compatibility entry point for the repaired unified competitive mapper.

The former hand-picked per-pair mini-reference implementation was scientifically
invalid and has been removed.  Use ``--sample`` (and the same options accepted by
stage4_map_sample_competitive.py); all reviewed pairs are evaluated against the
one validated full background reference.
"""
from stage4_map_sample_competitive import main

if __name__ == "__main__":
    raise SystemExit(main())
