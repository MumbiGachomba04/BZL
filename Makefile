CC=clang++
CFLAGS=-O3
LIBS=-lm
AUTOVFLAGS=-mllvm -combiner-store-merging=0 -mcpu=avispado -mllvm -vectorizer-use-vp-strided-load-store -mllvm -enable-mem-access-versioning=0 -mllvm -disable-loop-idiom-memcpy -fno-slp-vectorize
OMPFLAGS=-fopenmp
VECANALYSIS=-Rpass-analysis=loop-vectorize -Rpass=loop-vectorize
VFLAGS =-ffast-math -fno-builtin-sin -mepi

codes10 = lanczosomp lanczosserial lanczossimd lanczosonepointo.x
codes07 = lanczosomp lanczosserial lanczossimd lanczosvectorised.x 

all10: $(codes10)
all07: $(codes07)

lanczosserial: lanczosintegrated.cpp
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

lanczosomp: lanczosintegrated.cpp
	$(CC) $(CFLAGS) $^ -o $@ -D_OMP -D_OPENMPWORKSHARING $(OMPFLAGS) $(LIBS)

lanczossimd: lanczosintegrated.cpp
	$(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ -D_OMP -D_OPENMPSIMD $(OMPFLAGS) $(LIBS) $(VECANALYSIS)

lanczosvectorised.x: lanczosintegrated.cpp
	$(CC) $(CFLAGS) $(VFLAGS) $(AUTOVFLAGS) $^ -o $@ $(LIBS) $(VECANALYSIS)

lanczosonepointo.x: lanczosintegrated.cpp
	$(CC) $(CFLAGS) $(VFLAGS) $(AUTOVFLAGS) $^ -o $@ $(LIBS) $(VECANALYSIS)

clean:
	rm -f $(codes10) $(codes07) *.o *.x
	 

