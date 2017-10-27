CC = gcc
CXX = g++
CFLAGS = -O3 -ffast-math -std=c99 `pkg-config --cflags glib-2.0` -I../taco/include/
CXXFLAGS = -O3 -std=c++11 -DNDEBUG -DTACO `pkg-config --cflags glib-2.0` -DDECIMAL_DIG=17 -I../taco/include/
LDFLAGS = -I/usr/local/include -L/usr/local/lib

LDLIBS = -L../taco/build/lib -ltaco -lgsl -lgslcblas `pkg-config --libs glib-2.0` -lm -ldl

all: reference asx oski spmv
clean:
	rm -rf reference asx oski *.o *.dSYM *.trace *.pyc

reference: run_fill.o test_fill.o reference.o util.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

asx: run_fill.o test_fill.o asx.o util.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

oski: run_fill.o test_fill.o oski.o util.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

spmv: run_spmv.o test_spmv.o util.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)