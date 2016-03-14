CC=clang++
CFLAGS=-Wall `/opt/alicesw/fastjet/3.0.6/inst/bin/fastjet-config --cxxflags --libs --plugins`
SOURCES=TestFastJet.cpp

TestFastJet: TestFastJet.cpp TestFastJet.h
	$(CC) $(CFLAGS) TestFastJet.cpp -o TestFastJet
