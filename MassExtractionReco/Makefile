CC = g++
ROOTFLAGS = `root-config --cflags --glibs`

all: hist_filler envelope mass MassCalibration

hist_filler: src/hist_filler.C include/hist_filler.h
	$(CC) -g -o hist_filler src/hist_filler.C $(ROOTFLAGS)

envelope: src/envelope.cc include/envelope.h
	$(CC) -g -o envelope src/envelope.cc $(ROOTFLAGS)

mass: src/mass.cc include/mass.h
	$(CC) -g -o mass src/mass.cc $(ROOTFLAGS)

MassCalibration: src/MassCalibration.cc
	$(CC) -g -o MassCalibration src/MassCalibration.cc $(ROOTFLAGS)

clean:
	@rm -f hist_filler envelope mass MassCalibration
