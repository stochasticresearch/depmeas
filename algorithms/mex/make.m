% make the mex file

% stole the optimization flags from:
% https://github.com/tsogkas/matlab-utils/blob/master/make.m

%mex -v CXXOPTIMFLAGS="-O3 -DNDEBUG -fopenmp -march=native -O3 -ffast-math -pthread -pipe -msse2" LDOPTIMFLAGS="-O3 -fopenmp -march=native -O3 -ffast-math -pthread -pipe -msse2"  ktau_numer.c
mex -v COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ktau_numer.c