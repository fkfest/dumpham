# program name
MAIN = ../dumpham
OBJ = readdump

test : $(OBJ)

$(OBJ) :
	@test -e $@ || $(MAIN) $@.dh
	diff -q $@.FCIDUMP $@.ref 
	rm -f $@.FCIDUMP
	[[ ! -e $@.refo ]] || diff -q $@.ORBDUMP $@.refo
	rm -f $@.ORBDUMP
