CFLAGS+=-std=c99 -pedantic -Wall
LDFLAGS+=-lm
TESTSRC = MatlabUtils.o MatlabUtilsTest.o 
TESTBIN = test
ARTIFACTS = *.o $(TESTBIN)

all: test

test: run_test


run_test: build_test
	./test
build_test: MatlabUtils.o MatlabUtilsTest.o
	$(CC) $(CFLAGS) -o $(TESTBIN) $(TESTSRC) $(LDFLAGS)

clean: 
	rm $(ARTIFACTS)
