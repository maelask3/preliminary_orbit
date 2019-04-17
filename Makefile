CFLAGS+=-std=c99 -pedantic -Wall
LDFLAGS+=-lm
MATLABUTILS-SRC = MatlabUtils.o rpoly.o
TESTSRC = $(MATLABUTILS-SRC) MatlabUtilsTest.o 
TESTBIN = test
ARTIFACTS = *.o $(TESTBIN)

all: test

test: run_test


run_test: build_test
	./test
build_test: $(TESTSRC)
	$(CC) $(CFLAGS) -o $(TESTBIN) $(TESTSRC) $(LDFLAGS)

clean: 
	rm $(ARTIFACTS)
