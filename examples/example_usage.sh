#! /bin/bash

python examples/example_usage.py \
    '/nfs/winstor/margrie/SimonWeiler/RawData/Invivo_imaging/3photon_rotation/shared/' \
    '/ceph/margrie/laura/cimaut/' \
    --folder_read_pattern '2*' \
    --file_read_pattern 'rotation_00001.tif' \
    --file_read_pattern '*.bin' \
