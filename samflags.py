#!/usr/bin/env python

import sys
import errno

flag_code = {
    1: 'M',
    2: '',
    4: 'U',
    8: 'u',
    16: 'R',
    32: 'r',
    64: '1',
    128: '2',
    256: 'n',
    512: 'Q',
    1024: 'D',
    2048: 'N'
}

try:

    for l in sys.stdin:

        bits = l.split("\t")
        flags = int(bits[1])
        bit = 1
        flagstr = ""
        while flags != 0:
            if flags & bit:
                flagstr = flagstr + flag_code[bit]
                flags = flags & ~bit
            bit = bit * 2
        bits[1] = flagstr
        sys.stdout.write("\t".join(bits))
        
except IOError as e:
    if e.errno == errno.EPIPE:
        sys.exit(0)
    else:
        raise e
