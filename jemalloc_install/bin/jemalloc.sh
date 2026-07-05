#!/bin/sh

prefix=/home/okcut/FUSOR-hibff/jemalloc_install
exec_prefix=/home/okcut/FUSOR-hibff/jemalloc_install
libdir=${exec_prefix}/lib

LD_PRELOAD=${libdir}/libjemalloc.so.2
export LD_PRELOAD
exec "$@"
