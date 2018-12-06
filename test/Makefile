# program name
MAIN = ../dumpham
OBJ = readdump symuhfdump rdm

test : $(OBJ)

$(OBJ) :
	@test -e $@ || $(MAIN) $@.dh > $@.log
	@echo "check $@.FCIDUMP" && diff -q $@.FCIDUMP $@.ref 
	@rm -f $@.FCIDUMP
	@[[ ! -e $@.refo ]] || ( echo "check $@.ORBDUMP" && diff -q $@.ORBDUMP $@.refo)
	@rm -f $@.ORBDUMP $@.log
	@[[ ! -e $@.rdm1ref ]] || ( echo "check $@.1rdm" && diff -q $@.1rdm $@.rdm1ref)
	@rm -f $@.1rdm $@.log
