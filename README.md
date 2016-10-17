# GenBlastG
genBlast is a software developed by Rong She of the Jack Chen lab. The website describing it (including a web service, reference papers and the original source code) can be found here: http://genome.sfu.ca/genblast/

## Changes

this is based on the GenBlastG 1.39 codebase

* changed system calls to check the return status
* changed some command line parameter munging to use snprintf to prevent buffer overruns
* increased the maximum size of the blast command to 1024 bytes

## Reason
This is a simple patch to prevent GenBlastG from randomly segfaulting as part of the WormBase production pipeline.
