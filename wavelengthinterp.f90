!=======================================================================
!                               Begin Program
!                               wavelengthtrim
!=======================================================================
!
      program wavelengthtrim
!
!     written by: Aaron Gilli!h                    Tuesday, Feb 6, 2007
!
!     purpose: This file is adapted from Catherine's wavelength2.f code.
!              Trims a Phoenix .70 intensity field file and interpolates
!              it into smooth wavelength spacing so that it can
!              be run through Catherine's code.  Outputs the 
!              lte_atm00# files.
!
!-----------------------------------------------------------------------
!
      IMPLICIT none
       
      CHARACTER*40 flnm(200),  flnm1
      CHARACTER*40 flnm2
      CHARACTER*1  ten,  one
      
      integer  i,j,k,ii,ihun,iten,ione,icount,l,ntheta,nindex,n
      integer  iarray
      
      parameter (iarray=2000000)
      
      real*8     lamint(iarray)         !interpolated lambda array
      real*8     intensint(iarray,32)    !32 must match ntheta-careful
      real*8     lamold, lamnew
      real*8     intold(40), intnew(40)
      real*8     blambda,tlambda,x,f
      real*8     cosa(40)
      NAMELIST     /lambdas/ blambda,tlambda,x

!-----------------------------------------------------------------------
! Open the file wavelength that contains the boundaries between which
! we are trimming.  nml is namelist of the three values.

      OPEN(UNIT=10,FILE='wavelength.nml')
      READ(10,nml=lambdas)
      CLOSE(10)

!-----------------------------------------------------------------------
! Set up interpolated lambda array at desired wavelength spacing

      nindex = (tlambda - blambda)/x + 1.d0
      do i=1,nindex
         lamint(i) = blambda + (i-1)*x
!         write (*,*) lamint(i)
      enddo

!-----------------------------------------------------------------------
! Read fname3 namelist and assign files to flnm array

      icount = 0
         OPEN(UNIT=12,FILE='tmp/atmfilelist',STATUS='old')
         
         DO i = 1, 300
           READ(12,*,END=107) flnm(i)
           icount = icount+1 
         ENDDO
!-----------------------------------------------------------------------
         
107   continue
      DO ii = 1,icount     !begin loop over 'icount' number of files

         !--------------------------------------------------------------
         !This first part assigns a file name and number lte_atm00# for
         !working with more than one file.
         
         ihun = ii/100
         iten = (ii-ihun*100)/10
         ione = ii-iten*10-ihun*100
         flnm1 = flnm(ii)
         !flnm2 = 'lte_atm'//char(ihun+48)//char(iten+48)//char(ione+48)
         flnm2 = trim(flnm1)//'.i'
         write(*,*) trim(flnm1)//'.i'
         write(*,*) flnm1,flnm2
            
         OPEN(UNIT=10,FILE=flnm1,STATUS='OLD')
         OPEN(UNIT=11,FILE=trim(flnm2))

         !--------------------------------------------------------------
         !Write number of angles for intensity, usually 32
         
         READ(10,*) ntheta
         WRITE(11,101) ntheta    
101      FORMAT(5x,i4)
         
         !--------------------------------------------------------------
         !Read through first block.  The cosine of ntheta(32) angles
         
         DO j = 1, ntheta, 10    !begin read cosines loop
            l = j+9
               
            IF (l.GT.ntheta) l = ntheta
               
            IF (j.LE.ntheta) THEN
               READ(10,*) (cosa(k),k=j,l)
               WRITE(11,102) (cosa(k),k=j,l)
102            FORMAT(1x,10d12.4)
            ENDIF
            
         ENDDO   ! end read cosines loop
         
         !-------------------------------------------------------------
         !Take in first wavelength and intensity set, call it new
            
         READ(10,104) lamnew
         write(*,104) lamnew
          
         DO i = 1,ntheta,10
            j = i+9
            IF (j.GT.ntheta) j = ntheta
            IF (i.LE.ntheta) READ(10,*) (intnew(k),k=i,j)
            !intnew(17:32) = dlog10(intnew(17:32))
!            write(*,105) (intnew(k),k=i,j)
         ENDDO
            
            
         !--------------------------------------------------------------
         !Begin loop through all wavelengths from blambda to tlambda
         
         DO n = 1,nindex    !begin loop over wavelengths
500         continue
!             write (*,*) 'inside the index loop', n     
            !-----------------------------------------------------------
            !Take in next wavelength and intensity set, call it new and
            
            
            IF (lamint(n) .ge. lamnew) THEN
               
               lamold = lamnew
               DO i=1,ntheta
                  intold(i) = intnew(i)
               ENDDO
               
               
               
               READ(10,*) lamnew
               
!               write(*,*) lamold, lamnew
               
               DO i = 1,ntheta,10
                  j = i+9
                  IF (j.GT.ntheta) j = ntheta
                  IF (i.LE.ntheta) READ(10,*) (intnew(k),k=i,j)
                  !intnew(17:32) = dlog10(intnew(17:32))
               ENDDO
            ENDIF
            
            !-----------------------------------------------------------
            !Keep moving through lamnew until we get just past lamint(n)
            IF (lamnew .le. lamint(n)) THEN
               goto 500
            ENDIF
            
            
            !-----------------------------------------------------------
            !If point lamint is between lamold and lamnew: interpolate
            
!            IF (lamold .lt. lamint(n) .and.lamnew .gt. lamint(n)) THEN
               
               f = (lamint(n)-lamold)/(lamnew - lamold)
               DO j=1,ntheta
				  if(j>=17)then
					intensint(n,j)=(intnew(j)-intold(j))*f + intold(j)
					if(intensint(n,j).ne.0.d0) then
					intensint(n,j) = dlog10(intensint(n,j))
					else
					intensint(n,j) = 1d-20
					endif
				  else
					intensint(n,j)=0.d0
				  endif
               ENDDO
               
               write(11,104) lamint(n)
               
               DO i = 1,ntheta,10
                  j = i+9
                  IF (j.GT.ntheta) j = ntheta
                  WRITE(11,105) (intensint(n,k),k=i,j)
               ENDDO
               
104            FORMAT (1d30.16)
105            FORMAT (1x,10d12.4)
               
!            ENDIF

         ENDDO
             
1001     WRITE(*,*) 'End of file reached'
         CLOSE(11)
         CLOSE(10)
         
      ENDDO

      END
