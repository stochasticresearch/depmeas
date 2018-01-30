% make the mex file
mex -v COPTIMFLAGS="-O3 -fwrapv -DNDEBUG"  ktau_numer.c