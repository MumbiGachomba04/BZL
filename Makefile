CC=clang++
CFLAGS=-O3
LIBS=-lm
AUTOVFLAGS=-mllvm -combiner-store-merging=0 -mllvm -vectorizer-use-vp-strided-load-store -mllvm -enable-mem-access-versioning=0 -mllvm -disable-loop-idiom-memcpy -fno-slp-vectorize
OMPFLAGS=-fopenmp
VECANALYSIS=-Rpass-analysis=loop-vectorize -Rpass=loop-vectorize
VFLAGS =-ffast-math -fno-builtin-sin -mepi

codes10 = lanczosomp lanczosserial lanczossimd lanczosonepointo.x lanczosbothomp
codes07 = lanczosomp lanczosserial lanczossimd lanczosvectorised.x lanczosbothomp

all10: $(codes10)
all07: $(codes07)

lanczosserial: lanczosintegrated.cpp
	    $(CC) $(CFLAGS) $^ -o $@ $(LIBS)

lanczosomp: lanczosintegrated.cpp
	    $(CC) $(CFLAGS) $^ -o $@ -D_OMP -D_OPENMPWORKSHARING $(OMPFLAGS) $(LIBS)

lanczossimd: lanczosintegrated.cpp
	    $(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ -D_OMP -D_OPENMPSIMD $(OMPFLAGS) $(LIBS) $(VECANALYSIS)

lanczosvectorised.x: lanczosintegrated.cpp
	    $(CC) $(CFLAGS)  $(VFLAGS) -mcpu=avispado $(AUTOVFLAGS)  $(SDV_TRACE_INCL) $^ -o $@ -D_RAVETRACE $(LIBS) $(VECANALYSIS) $(SDV_TRACE_C_LINK)

lanczosonepointo.x: lanczosintegrated.cpp
	    $(CC) $(CFLAGS)  $(VFLAGS) -mcpu=spacemit-x60 $(AUTOVFLAGS) $(SDV_TRACE_INCL) $^ -o $@ -D_RAVETRACE $(LIBS) $(VECANALYSIS) $(SDV_TRACE_C_LINK)

lanczosbothomp: lanczosintegrated.cpp
	    $(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ -D_OMP -D_OPENMPBOTH $(OMPFLAGS) $(LIBS) $(VECANALYSIS)

clean:
	    rm -f $(codes10) *.o *.x
	    rm -f $(codes07) *.o *.x

	 

