      PROGRAM Main


c  A program to take an outputflux file such as those produced by CLIC
c  and generate apparent magnitudes in a few of the common filter systems.


      IMPLICIT NONE
      
      INTEGER nlambda,i,j,nwave,nincl
      PARAMETER(nlambda = 4549,nincl = 10)
c 8846, 1048992
      DOUBLE PRECISION lambda(nlambda),flux(nlambda,nincl)
      

      OPEN(UNIT=11,FILE='outputflux')

c  I'm not quite sure what the best way to do this is....I don't think I want
c  to read all 300 000 lines of my input file at once (although I suppose that
c  is feasible).  But I think I need to have enough wavelengths to process
c  an entire band at once.  I suppose that's not entirely true, I could 
c  use bins the same way I do in CLIC when I do my doppler shifting.
     
c  I think I will start by reading everything in and then processing each
c  filter set in it's own subroutine.


      DO i = 1,nlambda
         READ(11,*,END=104) lambda(i),(flux(i,j),j=1,nincl)
 101     FORMAT(1d30.16,11d14.6)
      ENDDO
 104  nwave = i
      write(*,*) 'read in ',i,' wavelengths'
      call johns(flux,lambda,nlambda,nwave)

      call strom(flux,lambda,nlambda,nwave)

      call walrav(flux,lambda,nlambda,nwave)

      END


      SUBROUTINE johns(flux,lambda,nlambda,nwave)

      IMPLICIT NONE

      INTEGER nwave,i,j,nlambda,jj,kk,m,k,n
c      PARAMETER(nlambda = 400000)
      DOUBLE PRECISION flux(nlambda,10),lambda(nlambda)
      DOUBLE PRECISION uwave(35),bwave(35),vwave(35),rwave(35)
      DOUBLE PRECISION iwave(35)
      DOUBLE PRECISION uflux(35),bflux(35),vflux(35),rflux(35)
      DOUBLE PRECISION iflux(35)
      DOUBLE PRECISION umag(10),bmag(10),vmag(10),rmag(10),imag(10)
      DOUBLE PRECISION fract,f


      
      OPEN(UNIT=12,FILE='johnson.dat')
      READ(12,*)
      READ(12,*)
      DO i = 1,26
         READ(12,*),uwave(i),uflux(i),bwave(i),bflux(i),vwave(i),
     &        vflux(i),rwave(i),rflux(i),iwave(i),iflux(i)
c  converts input wavelengths from microns to angstroms
         uwave(i) = uwave(i)*1e4
         bwave(i) = bwave(i)*1e4
         vwave(i) = vwave(i)*1e4
         rwave(i) = rwave(i)*1e4
         iwave(i) = iwave(i)*1e4
      ENDDO
c      write(*,*) 'read Johnson responses'
c  Calculate magnitudes:
      j = 1
      jj = 1
      k = 1
      kk = 1
      m = 1
c  This is just a first pass - produces a step function across each
c  section - I will need to do something to smooth the function - possibly
c  linear interpolation.


      DO i = 1,10
         umag(i) = 0
         bmag(i) = 0
         vmag(i) = 0
         rmag(i) = 0
         imag(i) = 0
      ENDDO
c  also, the question arises - these bands are all given at wavelengths 
c  spacings of 10 angstroms - is it really necessary to have such fine
c  spacing in my atmospheric models.
c      OPEN(UNIT=33,FILE='trans_U')
      DO i = 2,nwave
c  uband:
c      write(*,*) j,i
         IF (lambda(i).GE.uwave(j).AND.lambda(i).LT.uwave(j+1)) THEN
            fract = (lambda(i)-uwave(j))/(uwave(j+1)-uwave(j))
            f = fract*uflux(j+1)+(1-fract)*uflux(j)
            DO n = 1,10
            umag(n)=umag(n)+flux(i,n)*f*(lambda(i+1)-lambda(i-1))/2.
            ENDDO
c            WRITE(33,*)lambda(i),flux(i,1)*f,(lambda(i+1)-lambda(i-1))/2.
         ENDIF
         IF (lambda(i).GE.uwave(j+1)) j = j+1
c  b-band     

         IF (lambda(i).GE.bwave(jj).AND.lambda(i).LT.bwave(jj+1)) THEN
            fract = (lambda(i)-bwave(jj))/(bwave(jj+1)-bwave(jj))
            f = fract*bflux(jj+1)+(1-fract)*bflux(jj)
            DO n = 1,10
            bmag(n)=bmag(n)+flux(i,n)*f*(lambda(i+1)-lambda(i-1))/2.
            ENDDO
         ENDIF
         IF (lambda(i).GE.bwave(jj+1)) jj = jj+1
c  v-band         
         IF (lambda(i).GE.vwave(k).AND.lambda(i).LT.vwave(k+1)) THEN
            fract = (lambda(i)-vwave(k))/(vwave(k+1)-vwave(k))
            f = fract*vflux(k+1)+(1-fract)*vflux(k)
            DO n = 1,10
            vmag(n)=vmag(n)+flux(i,n)*f*(lambda(i+1)-lambda(i-1))/2.
            ENDDO
         ENDIF
         IF (lambda(i).GE.vwave(k+1)) k = k+1
c r - band
         IF (lambda(i).GE.rwave(kk).AND.lambda(i).LT.rwave(kk+1)) THEN
            fract = (lambda(i)-rwave(kk))/(rwave(kk+1)-rwave(kk))
            f = fract*rflux(kk+1)+(1-fract)*rflux(kk)
            DO n = 1,10
            rmag(n)=rmag(n)+flux(i,n)*f*(lambda(i+1)-lambda(i-1))/2.
            ENDDO
         ENDIF
         IF (lambda(i).GE.rwave(kk+1)) kk = kk+1
c  i-band           
         IF (lambda(i).GE.iwave(m).AND.lambda(i).LT.iwave(m+1)) THEN
            fract = (lambda(i)-iwave(m))/(iwave(m+1)-iwave(m))
            f = fract*iflux(m+1)+(1-fract)*iflux(m)
            DO n = 1,10
            imag(n)=imag(n)+flux(i,n)*f*(lambda(i+1)-lambda(i-1))/2.
            ENDDO
         ENDIF
         IF (lambda(i).GT.iwave(m+1)) m = m+1
         IF (j.GT.25) j = 25
         IF (jj.GT.25) jj = 25
         IF (k.GT.25) k = 25
         IF (kk.GT.25) kk = 25
         IF (m.GT.25) m = 25
      ENDDO

      OPEN(UNIT=13,FILE='jmagout')
      
      DO i = 1,10
         WRITE(13,102) umag(i),bmag(i),vmag(i),rmag(i),imag(i)
      ENDDO
 102  FORMAT(5d13.5)
      RETURN
      END


      SUBROUTINE strom(flux,lambda,nlambda,nwave)

      IMPLICIT NONE

      INTEGER nwave,i,j,nlambda,jj,kk,m,k,n
c      PARAMETER(nlambda = 400000)
      DOUBLE PRECISION flux(nlambda,10),lambda(nlambda)
      DOUBLE PRECISION uwave(50),bwave(50),vwave(50),ywave(50)
      DOUBLE PRECISION hbwave(50)
      DOUBLE PRECISION uflux(50),bflux(50),vflux(50),yflux(50)
      DOUBLE PRECISION hbflux(50)
      DOUBLE PRECISION umag(10),bmag(10),vmag(10),ymag(10),hbmag(10)
      DOUBLE PRECISION fract,f


      
      OPEN(UNIT=12,FILE='stromgren.dat')
      READ(12,*)
      READ(12,*)
      DO i = 1,50
         READ(12,*),uwave(i),uflux(i),vwave(i),vflux(i),bwave(i),
     &        bflux(i),ywave(i),yflux(i),hbwave(i),hbflux(i)
      ENDDO
c  Calculate magnitudes:
      j = 1
      jj = 1
      k = 1
      kk = 1
      m = 1

      DO i = 1,10
         umag(i) = 0
         bmag(i) = 0
         vmag(i) = 0
         ymag(i) = 0
         hbmag(i) = 0
      ENDDO

      DO i = 1,nwave
c  uband:
c      write(*,*) j,i
         IF (lambda(i).GE.uwave(j).AND.lambda(i).LT.uwave(j+1)) THEN
            fract = (lambda(i)-uwave(j))/(uwave(j+1)-uwave(j))
            f = fract*uflux(j+1)+(1-fract)*uflux(j)
            DO n = 1,10
               umag(n) = umag(n) + flux(i,n)*f
            ENDDO
         ENDIF
         IF (lambda(i).GE.uwave(j+1)) j = j+1
c  b-band     

         IF (lambda(i).GE.bwave(jj).AND.lambda(i).LT.bwave(jj+1)) THEN
            fract = (lambda(i)-bwave(jj))/(bwave(jj+1)-bwave(jj))
            f = fract*bflux(jj+1)+(1-fract)*bflux(jj)
            DO n = 1,10
               bmag(n) = bmag(n) + flux(i,n)*f
            ENDDO
         ENDIF
         IF (lambda(i).GE.bwave(jj+1)) jj = jj+1
c  v-band         
         IF (lambda(i).GE.vwave(k).AND.lambda(i).LT.vwave(k+1)) THEN
            fract = (lambda(i)-vwave(k))/(vwave(k+1)-vwave(k))
            f = fract*vflux(k+1)+(1-fract)*vflux(k)
            DO n = 1,10
               vmag(n) = vmag(n) + flux(i,n)*f
            ENDDO
         ENDIF
         IF (lambda(i).GE.vwave(k+1)) k = k+1
c y - band
         IF (lambda(i).GE.ywave(kk).AND.lambda(i).LT.ywave(kk+1)) THEN
            fract = (lambda(i)-ywave(kk))/(ywave(kk+1)-ywave(kk))
            f = fract*yflux(kk+1)+(1-fract)*yflux(kk)
            DO n = 1,10
               ymag(n) = ymag(n) + flux(i,n)*f
            ENDDO
         ENDIF
         IF (lambda(i).GE.ywave(kk+1)) kk = kk+1
c  i-band           
         IF (lambda(i).GE.hbwave(m).AND.lambda(i).LT.hbwave(m+1)) THEN
            fract = (lambda(i)-hbwave(m))/(hbwave(m+1)-hbwave(m))
            f = fract*hbflux(m+1)+(1-fract)*hbflux(m)
            DO n = 1,10
               hbmag(n) = hbmag(n) + flux(i,n)*f
            ENDDO
         ENDIF
         IF (lambda(i).GT.hbwave(m+1)) m = m+1
         IF (j.GT.49) j = 49
         IF (jj.GT.49) jj = 49
         IF (k.GT.49) k = 49
         IF (kk.GT.49) kk = 49
         IF (m.GT.49) m = 49
      ENDDO

      OPEN(UNIT=13,FILE='smagout')
      
      DO i = 1,10
         WRITE(13,102) umag(i),bmag(i),vmag(i),ymag(i),hbmag(i)
      ENDDO
 102  FORMAT(5d13.5)
      RETURN
      END


      SUBROUTINE walrav(flux,lambda,nlambda,nwave)

      IMPLICIT NONE

      INTEGER nwave,i,j,nlambda,jj,kk,m,k,n
c      PARAMETER(nlambda = 400000)
      DOUBLE PRECISION flux(nlambda,10),lambda(nlambda)
      DOUBLE PRECISION uwave(60),bwave(60),vwave(60),wwave(60)
      DOUBLE PRECISION lwave(60)
      DOUBLE PRECISION uflux(60),bflux(60),vflux(60),wflux(60)
      DOUBLE PRECISION lflux(60)
      DOUBLE PRECISION umag(10),bmag(10),vmag(10),lmag(10),wmag(10)
      DOUBLE PRECISION fract,f


      
      OPEN(UNIT=12,FILE='walraven.dat')
      READ(12,*)
      READ(12,*)
      DO i = 1,59
         READ(12,*),wwave(i),wflux(i),uwave(i),uflux(i),lwave(i),
     &        lflux(i),bwave(i),bflux(i),vwave(i),vflux(i)
      ENDDO
c  Calculate magnitudes:
      j = 1
      jj = 1
      k = 1
      kk = 1
      m = 1

      DO i = 1,10
         umag(i) = 0
         bmag(i) = 0
         vmag(i) = 0
         lmag(i) = 0
         wmag(i) = 0
      ENDDO

      DO i = 1,nwave
c  w-band:
c      write(*,*) j,i
         IF (lambda(i).GE.wwave(j).AND.lambda(i).LT.wwave(j+1)) THEN
            fract = (lambda(i)-wwave(j))/(wwave(j+1)-wwave(j))
            f = fract*wflux(j+1)+(1-fract)*wflux(j)
            DO n = 1,10
               wmag(n) = wmag(n) + flux(i,n)*f
            ENDDO
         ENDIF
         IF (lambda(i).GE.wwave(j+1)) j = j+1
c  u-band     

         IF (lambda(i).GE.uwave(jj).AND.lambda(i).LT.uwave(jj+1)) THEN
            fract = (lambda(i)-uwave(jj))/(uwave(jj+1)-uwave(jj))
            f = fract*uflux(jj+1)+(1-fract)*uflux(jj)
            DO n = 1,10
               umag(n) = umag(n) + flux(i,n)*f
            ENDDO
         ENDIF
         IF (lambda(i).GE.uwave(jj+1)) jj = jj+1
c  l-band         
         IF (lambda(i).GE.lwave(k).AND.lambda(i).LT.lwave(k+1)) THEN
            fract = (lambda(i)-lwave(k))/(lwave(k+1)-lwave(k))
            f = fract*lflux(k+1)+(1-fract)*lflux(k)
            DO n = 1,10
               lmag(n) = lmag(n) + flux(i,n)*f
            ENDDO
         ENDIF
         IF (lambda(i).GE.lwave(k+1)) k = k+1
c b - band
         IF (lambda(i).GE.bwave(kk).AND.lambda(i).LT.bwave(kk+1)) THEN
            fract = (lambda(i)-bwave(kk))/(bwave(kk+1)-bwave(kk))
            f = fract*bflux(kk+1)+(1-fract)*bflux(kk)
            DO n = 1,10
               bmag(n) = bmag(n) + flux(i,n)*f
            ENDDO
         ENDIF
         IF (lambda(i).GE.bwave(kk+1)) kk = kk+1
c  v-band           
         IF (lambda(i).GE.vwave(m).AND.lambda(i).LT.vwave(m+1)) THEN
            fract = (lambda(i)-vwave(m))/(vwave(m+1)-vwave(m))
            f = fract*vflux(m+1)+(1-fract)*vflux(m)
            DO n = 1,10
               vmag(n) = vmag(n) + flux(i,n)*f
            ENDDO
         ENDIF
         IF (lambda(i).GT.vwave(m+1)) m = m+1
         IF (j.GT.58) j = 58
         IF (jj.GT.58) jj = 58
         IF (k.GT.58) k = 58
         IF (kk.GT.58) kk = 58
         IF (m.GT.58) m = 58
      ENDDO

      OPEN(UNIT=13,FILE='wmagout')
      
      DO i = 1,10
         WRITE(13,102) wmag(i),umag(i),lmag(i),bmag(i),vmag(i)
      ENDDO
 102  FORMAT(5d13.5)
      RETURN
      END
