#!/bin/sh
# Copy fragit share data into the pixi env so the fragit CLI can find them.
# Required because editable pip installs do not copy data_files from setup.py.
FRAGIT_SHARE="$PIXI_PROJECT_ROOT/../fragit-main/share"
if [ -d "$FRAGIT_SHARE" ]; then
    cp -ru "$FRAGIT_SHARE/." "$CONDA_PREFIX/share/"
fi
