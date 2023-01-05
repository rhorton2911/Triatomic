      module pathpot

      PARAMETER (NDIM12=4,NDIM22=20000+1,NDIM32=100,LLMAX2=20) 

      real*8:: RFULL(NDIM22),RHALF(2*NDIM22),Y(NDIM12,NDIM22)
      real*8:: EGUESS,CORE(2*NDIM22),RDMASS

      end module pathpot
