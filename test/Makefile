# program name
MAIN = ../dumpham
OBJ = readdump

test : $(OBJ)

$(OBJ) :
	@test -e $@ || $(MAIN) $@.dh > $@.log
	@echo "check $@.FCIDUMP" && diff -q $@.FCIDUMP $@.ref 
	@rm -f $@.FCIDUMP
	@[[ ! -e $@.refo ]] || ( echo "check $@.ORBDUMP" && diff -q $@.ORBDUMP $@.refo)
	@rm -f $@.ORBDUMP $@.log
