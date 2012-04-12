# $* is prefix shared by target and dependent;  $@ is name of target file

CFLAGS = -c -O4 -Isamtools -Iboost -Ibamtools/src
OBJS= main.o  
NAME= CNV

$(NAME): $(OBJS)
	g++ -o $@ $(OBJS)  fastahack/Fasta.cpp fastahack/split.cpp  -lbamtools -Lbamtools/lib -lz -Lzlib 

# Suffix rules: tell how to  take file with first suffix and make it into
#	file with second suffix
	
.cpp.o:
	g++ $(CFLAGS) $*.cpp
	
clean:
	rm *.o
