#!/usr/bin/env python

import tempfile
import os
import glob
import sys
import subprocess
import shutil

if len(sys.argv) < 3:
    print >>sys.stderr, "Usage: tagAndMergeLanes.py output.bam lanefile1.bam [lanefile2.bam ...]"
    sys.exit(1)

workdir = tempfile.mkdtemp()

fifos = []
cmds = []

for f in sys.argv[2:]:

    bn = os.path.basename(f)
    bits = bn.split("_")
    try:
        sampleid = "_".join(bits[:-3])
        laneid = bits[-3]
    except IndexError:
        print >>sys.stderr, "Filename", f, "is not in the expected samplename_L???_R?_001.bam format"
        sys.exit(1)
    fifo = os.path.join(workdir, bn)
    fifos.append(fifo)
    cmds.append(["/usr/bin/java", "-jar", "/apps/modules/pkg/apps/picardtools/1.96/noarch/AddOrReplaceReadGroups.jar", "I=%s" % f, "O=%s" % fifo, "RGID=%s" % laneid, "RGSM=%s" % sampleid, "RGPL=illumina", "RGLB=%s" % sampleid, "RGPU=unit1", "COMPRESSION_LEVEL=0"])
    
for f in fifos:
    os.mkfifo(f)

cmds.append(["/home/csmowton/nonroot/bin/samtools", "merge", "-f", sys.argv[1]] + fifos)

procs = []

for c in cmds:
    print >>sys.stderr, c
    procs.append(subprocess.Popen(c))

rc = 0

for p in procs:
    if p.wait() != 0:
        rc = p.returncode

if rc != 0:
    print >>sys.stderr, "At least one child failed"

shutil.rmtree(workdir)

sys.exit(rc)
