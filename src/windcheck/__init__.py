"""windcheck — cross-winding consistency checking for scroll segmentations.

Layers, bottom-up:

    catalog.py   read the public S3 catalog; enumerate samples/volumes/segments
    fetch.py     download tifxyz surface maps; write data/MANIFEST.json
    tifxyz.py    read a tifxyz surface: (u,v) grid -> (x,y,z) volume coords
    atlas.py     build the winding reference atlas from the w052-w095 segments
    detect.py    the invariant: winding index along a trace must be monotone
                 and continuous; a discontinuity certifies a sheet-switch
    cli.py       command-line entry point

Nothing here depends on the upstream `vesuvius` package. That is deliberate:
the tool audits that pipeline, so it must be able to disagree with it.
"""

__version__ = "0.1.0"

S3_BUCKET = "vesuvius-challenge-open-data"
CATALOG_KEY = "metadata.json"  # gzipped JSON at the bucket root
