CC = g++
ROOTFLAGS = `root-config --cflags --glibs`


INCDIR   = include
SRCDIR   = src

SOURCES  := $(wildcard $(SRCDIR)/*.cc)
INCLUDES := $(wildcard $(INCDIR)/*.h)

MACROS := AcceptancePuritySubjets SubjetFlavor resolution_subjets_flavorJEC resolution_subjets_flavorJEC_updown tt_envelope tt_masspoints tt_pdf tt_sys

all: $(MACROS) # clean first and build always everything


$(MACROS): % : $(SRCDIR)/%.cc $(INCLUDES)
	@echo 'building' $@
	@$(CC) -o $@ $(SRCDIR)/$@.cc $(ROOTFLAGS)


clean:
	@rm -f $(MACROS)
