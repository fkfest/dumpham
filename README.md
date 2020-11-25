Dump various model Hamiltonians as FCIDUMP files.

  *  print help:

        dumpham -h

  *  rewrite FCIDUMP removing symmetry and adding all zero elements back:
  
        dumpham -d <FCIDUMP> [<FCIDUMP.NEW>]
        

  *  rewrite FCIDUMP removing symmetry and adding all zero elements back + create a file with orbital coefficients:
  

        dumpham -do <FCIDUMP> [<FCIDUMP.NEW> [<ORBDUMP>]]
        

  *  use dumpham-input file:
  

        dumpham <input.dh>


  *  the input parameters can be set from the command line (will be overwritten by the input file), e.g.,
    
  
        dumpham -i "hubbard,U=2.0" -i "ham,out=hubbard.fcidump" hubbard.dh


