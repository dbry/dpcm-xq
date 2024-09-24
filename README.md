## DPCM-XQ

Xtreme Quality DPCM Encoder/Decoder

Copyright (c) 2024 David Bryant.

All Rights Reserved.

Distributed under the [BSD Software License](https://github.com/dbry/dpcm-xq/blob/master/license.txt).

## What is this?

DPCM is a variation on PCM audio where regular 16-bit samples are stored
in just 8 bits. This provides 2:1 compression and has generally very good
quality. The encoded values are essentially the deltas (differences)
between one sample and the next, although there are many variations with
respect to the scaling of the deltas and how synchronization is implemented
(if it is). Other than its use in canned samples for games and other
applications (because of its high-quality and very low decoding complexity),
DPCM was never popular as a general-purpose audio format and there is
essentially no support for it in conventional containers (other than SDX2
in AIFF).

This encoder combines two different techniques to achieve higher quality
than existing DPCM encoders while remaining fully compatible with standard
decoders. The first is dynamic noise shaping, which shifts the quantization
noise up or down in frequency based on the spectrum of the source signal.
This technique is identical to the algorithm used in WavPack's lossy mode.

The other technique is "lookahead" in which the encoder exhaustively
searches ahead to find the optimum coding sequence based on future samples.
This process can reduce the quantization noise by around 3dB at most frequencies
(and 1 dB at high frequencies). Unfortunately, at its maximum settings this can
be very slow, but this should be relatively irrelevant if the encoder is being
used to generate canned samples.

**DPCM-XQ** consists of two standard C files, [dpcm.c](dpcm.c) and [dns.c](dns.c),
and is a stand-alone command-line program.  Because DPCM is not well-supported in
containers, the input and output is always raw, continuous PCM or DPCM with no
headers or frames. Pipes are fully supported for both input and output.

## Variations

There are several forms and variations of DPCM encoding. The two standard ones
handled by **DPCM-XQ** are referred to as **SDX2** (Squareroot-Delta-Exact) and
**CBD2** (Cuberoot-Delta-Exact), both in stereo and mono. Also I have created
another variation I call **L2XP** (Linear-to-Exponential) that provides better
S/N at lower levels in exchange for worse S/N at the highest levels by employing
a delta scaling table that is linear up to the tangent point of an exponential
curve (but that can still be generated in advance with only 32-bit integer math).

## Building

To build the command-line tool (**DPCM-XQ**) on Linux:

> $ gcc -O3 *.c -lm -o dpcm-xq

on Darwin/Mac:

> $ clang -O3 *.c -lm -o dpcm-xq

on MS Visual Studio:

> C:\cl -O3 dpcm.c dns.c

## Help

```
 DPCM-XQ  Xtreme Quality Raw DPCM Encoder / Decoder  Version 0.1
 Copyright (c) 2024 David Bryant. All Rights Reserved.

 Usage:     DPCM [-d|-e] [-options] infile.raw outfile.raw

 Operation: Specify '-' for either file spec for piped operation.
            Mode is 8-bit raw SDX2 stereo unless overridden below.

 Options:  -[0-16]= encode lookahead samples (default = 3, max = 16)
           -c     = override default mode and use CBD2 table instead
           -d     = decode operation (raw 8-bit DPCM to raw 16-bit PCM)
           -e     = encode operation (raw 16-bit PCM to raw 8-bit DPCM)
           -f     = encode flat noise (no noise shaping, aka -s0.0)
           -h     = display this help message
           -l     = override default mode and use L2XP table instead
           -m     = override stereo default and treat files as mono
           -n     = measure and report quantization noise
           -q     = quiet mode (display errors only)
           -s<n>  = enable noise shaping, (-1.0 < n <= 1.0)
           -v     = verbose (display lots of info)
           -y     = overwrite outfile if it exists

 Web:       Visit www.github.com/dbry/dpcm-xq for latest version and info

```

## Caveat

- In some situations, at high lookahead levels, the operation can get very slow
or even seem to be stuck. The default level 3 should always be fine and then the user
can simply try increasing levels until the time becomes untenable. The quantization
noise option (**-n**) can be used to determine if higher levels are providing
improvement (lower numbers are better). Note that the flat noise option (**-f**)
will provide the lowest *measured* noise, but the default dynamic noise shaping
may provide *less audible* noise.
