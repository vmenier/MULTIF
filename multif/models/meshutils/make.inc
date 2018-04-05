OPTIM=oui

# Rq:
# Le -statique n'a pas de sens sous mac 

ifeq  ($(ARCHI),ppcg4) 
 
  ifndef BIN
    BIN=$(HOME)/bin/ppcg4
  endif  
   
  ifeq ($(COMP),ibm)
  
    CC=xlc
    F77=xlf
    F90=xlf95
    LINK=$(CC)
    LIBF=/opt/ibmcmp/xlf/8.1/lib/
    LIB= -static $(LIBF)/*.a
  
    ifeq ($(OPTIM),oui)
      ifdef IPA
         CFLAGS=   -O3 -D __ibm__     -qhot -qipa -qinline  
         F90FLAGS= -O3 -qsuffix=f=f90 -qhot -qipa -qinline  
         F77FLAGS= -O3                -qhot -qipa -qinline 
       else
         CFLAGS=   -O3 -D __ibm__     -qhot -qinline  
         F90FLAGS= -O3 -qsuffix=f=f90 -qhot -qinline  
         F77FLAGS= -O3                -qhot -qinline 
       endif        
          	 
    else
      CFLAGS=    -O -D __ibm__     
      F90FLAGS=  -O -qsuffix=f=f90 
      F77FLAGS=  -O  
    endif
    
  else
  
    F77=g95
    F90=g95
    CC=gcc
    LINK     = $(F90)
    LIB=   -lSystemStubs 
    
    ifeq ($(OPTIM),oui)
      CFLAGS=   -O5     
      F90FLAGS= -O5 -cpp   -fsloppy-char -fno-underscoring
      F77FLAGS= -O5  -cpp   -fsloppy-char -fno-underscoring 
    else
      CFLAGS=    -g 
      F90FLAGS=  -cpp  -fsloppy-char -fno-underscoring  -fbounds-check -ftrace=full
      F77FLAGS=  -cpp  -fsloppy-char -fno-underscoring  -fbounds-check -ftrace=full  
    endif  
  
  
  endif 
endif     
   
ifeq  ($(ARCHI), ppcg5) 
 
  ifndef BIN
    BIN=$(HOME)/bin/ppcg5
  endif 
   
  ifeq ($(COMP), ibm)
  
    CC=xlc
    F77=xlf
    F90=xlf95
    LINK=$(CC)
    LIBF=/opt/ibmcmp/xlf/8.1/lib/
    LIB= -static $(LIBF)/*.a
  
    ifeq ($(OPTIM),oui)
      CFLAGS=   -O3 -D __ibm__     -qhot -qipa -qinline -qtune=g5 -qarch=g5
      F90FLAGS= -O3 -qsuffix=f=f90 -qhot -qipa -qinline -qtune=g5 -qarch=g5
      F77FLAGS= -O3                -qhot -qipa -qinline -qtune=g5 -qarch=g5
    else
      CFLAGS=    -O -D __ibm__      -qtune=g5 -qarch=g5
      F90FLAGS=  -O -qsuffix=f=f90  -qtune=g5 -qarch=g5
      F77FLAGS=  -O  
    endif
    
  else
  
    F77=g95
    F90=g95
    CC=gcc-4.2
    LINK     = $(F77)
    LIB=  -lSystemStubs 
    
    ifeq ($(OPTIM),oui)
      CFLAGS=   -fast -m64     
      F90FLAGS= -fast -m64  -cpp  -fsloppy-char -fno-underscoring
      F77FLAGS= -fast -m64  -cpp  -fsloppy-char -fno-underscoring 
    else
      CFLAGS=    -g 
      F90FLAGS=  -cpp  -fsloppy-char -fno-underscoring  -fbounds-check -ftrace=full
      F77FLAGS=  -cpp  -fsloppy-char -fno-underscoring  -fbounds-check -ftrace=full  
    endif  
  
  
  endif 
endif     
   

ifeq ($(ARCHI), macintel)
#  MAC G5 COMPILO GNU g95
  
  ifndef BIN
    BIN=$(HOME)/bin/macintel/
  endif  
  
  
  ifeq ($(COMP), intel)
   
    CC=icc
    F77=ifort
    F90=ifort
    LIBF=/opt/intel/fce/10.1.015/lib
    LINK=$(CC)
    LIB=  $(LIBF)/*.a 
    
    ifeq ($(OPTIM),oui) 
#      CFLAGS= -fast -DF77_NO_UNDER_SCORE -mp1
#      CFLAGS= -O3 -DF77_NO_UNDER_SCORE -mp1 -ipo
			 CFLAGS=   -O3 -m64  -finline-limit=10000 --param inline-unit-growth=1000 --param large-function-growth=1000 -Wall  -Wno-unused-function -Wno-unused-but-set-variable -Wno-inline
      F90FLAGS= -O3 -mp1 -ipo  -nus
      F77FLAGS= -O3 -mp1 -ipo  -nus
    else
      CFLAGS= -g -DF77_NO_UNDER_SCORE
      F90FLAGS= -check all  -nus 
      F77FLAGS= -check all  -nus
    endif  
   
   
  else
   
    F77=gfortran
    F90=gfortran
    CC=gcc
    LINK= $(F90)
    LIB=  -lSystemStubs 
    ifeq ($(OPTIM),oui)
      CFLAGS=   -O3 -m64  -Winline -finline-limit=5000 --param inline-unit-growth=200 --param large-function-growth=200 -Wall
#      CFLAGS=   -fast -m64  -Winline -finline-limit=5000 --param inline-unit-growth=200 --param large-function-growth=200 -Wall
      F90FLAGS= -fast -m64  -cpp   -fsloppy-char -fno-underscoring
      F77FLAGS= -fast -m64  -cpp   -fsloppy-char -fno-underscoring 
    else
      CFLAGS=    -g
      F90FLAGS=  -cpp  -fsloppy-char -fno-underscoring  -fbounds-check -ftrace=full
      F77FLAGS=  -cpp  -fsloppy-char -fno-underscoring  -fbounds-check -ftrace=full  
    endif  
   
    
   endif
endif



ifeq ($(ARCHI), linux32)
#  MAC G5 COMPILO GNU g95

  ifndef BIN
    BIN=$(HOME)/bin/linux32
  endif

  ifeq ($(COMP), intel )
    CC=iccbin
    F77=ifortbin
    F90=ifortbin
    LIBF=/opt/intel/fc/9.0/lib
    LINK=$(CC)
    LIB= -static $(LIBF)/*.a 
    
    ifeq ($(OPTIM),oui) 
      CFLAGS= -O3 -DF77_NO_UNDER_SCORE -mp1 -ipo 
      F90FLAGS= -O3 -mp1 -ipo   -nus
      F77FLAGS= -O3 -mp1 -ipo   -nus
    else
      CFLAGS= -g -DF77_NO_UNDER_SCORE
      F90FLAGS= -check all  -nus 
      F77FLAGS= -check all  -nus
    endif  
  
  else
  
    F77=g95
    F90=g95
    CC=gcc
    LINK= $(F90)
    LIB= -static  
    ifeq ($(OPTIM),oui) 
      CFLAGS= -O5 -DF77_NO_UNDER_SCORE  
      F90FLAGS= -O5 -cpp   -fsloppy-char -fno-underscoring
      F77FLAGS= -O5 -cpp   -fsloppy-char -fno-underscoring
    else
      CFLAGS= -g -DF77_NO_UNDER_SCORE
      F90FLAGS= -g -fbounds-check -ftrace=full  -fsloppy-char -fno-underscoring
      F77FLAGS= -g -fbounds-check -ftrace=full  -fsloppy-char -fno-underscoring
    endif  
  
  endif
endif 
  

ifeq ($(ARCHI), linux64)
#  MAC G5 COMPILO GNU g95

  ifndef BIN
    BIN=$(HOME)/bin/linux64
  endif  


  ifeq ($(COMP), intel)
    CC=icc
    F77=ifort
    F90=ifort
    LIBF=/opt/intel/fce/9.1.036/lib
    LINK=$(CC)
    LIB= -static $(LIBF)/*.a 
    
    ifeq ($(OPTIM),oui) 
      CFLAGS=   -fast -DF77_NO_UNDER_SCORE -mp1 -ipo 
      F90FLAGS= -fast -mp1 -nus
      F77FLAGS= -fast -mp1 -nus
    else
      CFLAGS= -g -DF77_NO_UNDER_SCORE
      F90FLAGS= -check all  -nus 
      F77FLAGS= -check all  -nus
    endif  
  
  else
  
    F77=g95
    F90=g95
    CC=gcc
    LINK= $(F90)
    LIB= -static  
    ifeq ($(OPTIM),oui) 
      CFLAGS=   -O5 -m64 -DF77_NO_UNDER_SCORE  
      LDFLAGS=  -lpthread
      F90FLAGS= -fast -m64 -cpp -fsloppy-char -fno-underscoring
      F77FLAGS= -fast -m64 -cpp -fsloppy-char -fno-underscoring
    else
      CFLAGS= -g -DF77_NO_UNDER_SCORE
      F90FLAGS= -g -fbounds-check -ftrace=full  -fsloppy-char -fno-underscoring
      F77FLAGS= -g -fbounds-check -ftrace=full  -fsloppy-char -fno-underscoring
    endif  
  
  endif
endif 
  
