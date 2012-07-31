# @(#)makefile	2-99
# $Revision:$

.SUFFIXES: .o .C .c .cpp

# FLAGS=  -O2
FLAGS= -g -std=c++0x

# BINDIR = /homes/daniel/bin
# SRCDIR = /homes/daniel/src
SRCDIR = .
BINDIR = .
OBJDIR = obj/

HEADER = utils.h  obs.h obsSeq.h stateTrans.h plainStateTrans.h gammaProb.h explicitDurationTrans.h obsProb.h discreteObsProb.h  gaussianObsProb.h vectorObsProb.h flexibleObsProb.h initStateProb.h hmm.h snp.h read.h

OBJ = utils.o readConfigFile.o obsSeq.o plainStateTrans.o gammaProb.o explicitDurationTrans.o initStateProb.o hmm.o discreteObsProb.o gaussianObsProb.o vectorObsProb.o flexibleObsProb.o snp.o read.o

OBJ3 = program.o

#===============================================================================


full : $(HEADER) $(OBJ) $(OBJ3)
	   g++ $(FLAGS) -o $(BINDIR)/fullprogram $(OBJ) $(OBJ3) -lm
	   mv *.o $(OBJDIR)

program : $(HEADER) $(OBJ) $(OBJ3)
	   g++ $(FLAGS) -o $(BINDIR)/program $(OBJ) $(OBJ3) -lm
	   # mv *.o $(OBJDIR)

hmmPur : $(HEADER) $(OBJ) $(OBJ3)
	   purify g++ $(FLAGS) -o $(BINDIR)/hmmPur $(OBJ) $(OBJ3) -lm


#===============================================================================

#===============================================================================

snp.o: snp.h

read.o: read.h

utils.o: utils.h

obsSeq.o: obs.h obsSeq.h read.h

plainStateTrans.o: stateTrans.h plainStateTrans.h

gammaProb.o: gammaProb.h

explicitDurationTrans.o: stateTrans.h gammaProb.h explicitDurationTrans.h 

discreteObsProb.o: obs.h obsProb.h discreteObsProb.h read.h

gaussianObsProb.o: obs.h obsProb.h gaussianObsProb.h read.h

vectorObsProb.o: obs.h obsProb.h gaussianObsProb.h vectorObsProb.h read.h discreteObsProb.h

flexibleObsProb.o: obs.h obsProb.h gaussianObsProb.h flexibleObsProb.h read.h discreteObsProb.h

initStateProb.o: initStateProb.h

hmm.o: utils.h obs.h obsSeq.h stateTrans.h plainStateTrans.h gammaProb.h explicitDurationTrans.h obsProb.h discreteObsProb.h gaussianObsProb.h initStateProb.h hmm.h vectorObsProb.h

program.o: utils.h obs.h obsSeq.h stateTrans.h plainStateTrans.h gammaProb.h explicitDurationTrans.h obsProb.h discreteObsProb.h gaussianObsProb.h initStateProb.h hmm.h read.h snp.h

#===============================================================================

#===============================================================================

.cpp.o: $(HEADERI)
	g++ -c  $(FLAGS) $<
.C.o: $(HEADERI)
	g++ -c  $(FLAGS) $<

#checkout:
#	co -l $(REV) $(ALL)

#checkin:
#	ci $(REV)  $(ALL)

clean:
	rm -f *.o core tags

tags: $(SRC)
	ctags *.c

lint:	rm -f lint.out
	lint *.c > lint.out

tar.all: $(SRC) $(DOC) $(HEADERS) makefile
	tar cvf tar.all $(SRC) $(DOC) $(HEADERS) makefile

#===============================================================================
#===============================================================================
