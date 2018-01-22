#include "FCIdump.h"

#ifdef __cplusplus
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#else
#include <stdio.h>
#endif

#ifdef FCIDUMP_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
//int main()
{
  int parallel_size=1, parallel_rank=0;
#ifdef FCIDUMP_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&parallel_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&parallel_rank);
  //  if (parallel_rank > 0) freopen("/dev/null", "w", stdout);
#endif

  int oneelectron[] = {2,1};
  int twoelectron[] = {4,3,2,1};
  std::vector<std::string> files(1,"H2SYM.FCIDUMP"); //files.push_back("uhf.fcidump");
  int ifile;
  if (parallel_rank == 0) {
      std::cout << std::endl<<"Process file "<<files[ifile]<<std::endl;
      FCIdump dump(files[ifile]);
      std::vector<int> NELEC=dump.parameter("NELEC");
      std::cout << "NELEC="<<NELEC[0]<<std::endl;
      std::vector<int> MS2=dump.parameter("MS2");
      std::cout << "MS2="<<MS2[0]<<std::endl;
      std::vector<int> NORB=dump.parameter("NORB");
      std::cout << "NORB="<<NORB[0]<<std::endl;
      std::vector<int> IUHF=dump.parameter("IUHF");
      std::cout << "IUHF="<<IUHF[0]<<std::endl;
      std::vector<int> ORBSYM=dump.parameter("ORBSYM");
      std::cout << "ORBSYM="; for (std::vector<int>::const_iterator s=ORBSYM.begin(); s<ORBSYM.end(); s++) std::cout <<*s<<","; std::cout<<std::endl;
      int i,j,k,l;
      double value;
      int nn = NORB[0];
      int n1el = nn*(nn+1)/2;
      int n2el = n1el*(n1el+1)/2;
      std::cout << "n2el: " << n2el << std::endl;
        
      
      std::vector<double> twoel(n2el,0.0),oneel(n1el,0.0);
      double Ecore = 0.0;
      FCIdump::integralType type;
      dump.rewind();
      
      while ((type = dump.nextIntegral(i,j,k,l,value)) != FCIdump::endOfFile) {
        if ( type == FCIdump::endOfRecord ) {
            // store
            for (i = 1; i <= nn; ++i )
                for (j = 1; j <= i; ++j )
                    for (k = 1; k <= nn; ++k )
                        for (l = 1; l <= k; ++l ) {
                            int ij = (i-1)*i/2 + j;
                            int kl = (k-1)*k/2 + l;
                            if ( kl <= ij ){
                                int ijkl = (ij-1)*ij/2 + kl - 1;
                                std::cout << twoel[ijkl] << ": " << i << " " << j << " " << k << " " << l << std::endl;
                            }
                        }
        }
        if ( k != 0 && l !=0 ) {
            // Two-electron integrals
            int ij = (i-1)*i/2 + j;
            int kl = (k-1)*k/2 + l;
            int ijkl = (ij-1)*ij/2 + kl - 1;
            twoel[ijkl] = value;
        } else if ( i != 0 && j != 0 ) {
            // One-electron integrals
            int ij = (i-1)*i/2 + j - 1;
            oneel[ij] = value;
        } else if ( type == FCIdump::I0 ){
            Ecore = value;
        }
    }
    if (dump.write("new.fcidump",FCIdump::FileFormatted,false))
        std::cout << "will be written to new file"<<std::endl;
    else
        std::cout << "failure to write to new file"<<std::endl;
    for (i = 1; i <= nn; ++i)
        for (j = 1; j <= i; ++j )
            for (k = 1; k <= nn; ++k )
                for (l = 1; l <= k; ++l ) {
                    int ij = (i-1)*i/2 + j;
                    int kl = (k-1)*k/2 + l;
                    if ( kl <= ij ){
                        int ijkl = (ij-1)*ij/2 + kl - 1;
                        dump.writeIntegral(i,j,k,l,twoel[ijkl]);
//                         std::cout<< twoel[ijkl] << ": " << i << " " << j << " " << k << " " << l << std::endl;
                    }
                }
    for (i = 1; i <= nn; ++i)
        for (j = 1; j <= i; ++j ) {
            int ij = (i-1)*i/2 + j - 1;
            dump.writeIntegral(i,j,0,0,oneel[ij]);
        }
    dump.writeIntegral(0,0,0,0,Ecore); 
  }

#ifdef FCIDUMP_MPI
  MPI_Finalize();
#endif
  return 0;
}
