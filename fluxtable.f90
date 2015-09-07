!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!

subroutine fluxtable(ftable,wav,angles,atmfname,nfnames,ntheta,nwav)
implicit none
!Pythonic array indices, from 0 to n-1.
integer (kind=4) :: nwav,ntheta,nfnames
!f2py depend(nzth) nzth
!f2py depend(nzphi) nzphi
integer (kind=4) :: i,j,k,n
real (kind=8), dimension(0:nwav-1), intent(inout) :: wav
!f2py depend(nwav) wav
!real (kind=8), dimension(0:nwav-1,0:1) :: tempflux
! depend(nwav) tempflux
real (kind=8), dimension(0:nwav-1,0:ntheta-1,0:nfnames-1), intent(inout) :: ftable
!f2py depend(nwav,ntheta,nfnames) fluxtable
!real (kind=8), dimension(0:nwav-1,0:nzphi-1,0:nzth-1), intent(inout) :: fout
! depend(nwav,nzphi,nzth) fout
real (kind=8), dimension(0:ntheta-1) :: angles
!f2py depend(ntheta) angles
character (len=40), intent(in) :: atmfname
character (len=40) :: tempfname
!f2py intent(in,out) :: ftable
!f2py intent(in,out) :: wav


open(unit=66,file="tmp/"//trim(atmfname),status="OLD")
do i = 0,nfnames-1
    read(66,*) tempfname
    tempfname = trim(tempfname)//'.i'
    open(unit=7,file=tempfname,status="OLD")
    read(7,*) n
    read(7,*) (angles(j),j=0,ntheta-1)
    do j = 0,nwav-1
        read(7,*) wav(j)
        read(7,*) (ftable(j,k,i),k=0,ntheta-1)
    enddo
enddo


close(66)
close(7)
return
end subroutine  ! fluxtable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!1!!1!!!!!1!!!!!!!!!!!!!!!!!
subroutine interp(tempflux,ftable,wav,angles,ext,dist,nfnames,ntheta,nwav,nzth,nzphi,incl)
implicit none
integer (kind=4) :: i,j,k,l,m,theta1,theta2,n
integer (kind=4) :: nwav,ntheta,nfnames,nzth,nzphi
real (kind=8) :: incl,p1,p2,p3,p4,ainte,dist,f11,f12,f21,f22,f31,f32,f41,f42
real (kind=8) :: res,deltaw,wavn,temp_iout
real (kind=8), dimension(0:nwav-1), intent(in) :: wav
!f2py depend(nwav) wav
real (kind=8), dimension(0:nwav-1,0:ntheta-1,0:nfnames-1), intent(in) :: ftable
!f2py depend(nwav,ntheta,nfnames) ftable
real (kind=8), dimension(0:nwav-1,0:nzphi-1,0:nzth-1) :: iout
!f2py depend(nwav,nzphi,nzth) iout
real (kind=8), dimension(0:nwav-1,0:1), intent(out) :: tempflux
!f2py depend(nwav) tempflux
real (kind=8), dimension(0:ntheta-1) :: angles
!f2py depend(ntheta) angles
real (kind=8), dimension(0:nzphi-1,0:nzth-1) :: cossquiggle
!f2py depend(nzphi,nzth) cossquiggle
real (kind=8), dimension(0:nzphi-1,0:nzth-1) :: doppler_arr
!f2py depend(nzphi,nzth) doppler_arr
integer (kind=4), dimension(0:nzth-1,0:3) :: atmkeys
!f2py depend(nzth) atmkeys
real (kind=8), dimension(0:nzth-1,0:5) :: tandg
!f2py depend(nzth) tandg
character (len=20) :: ext
!f2py intent(in) :: wav,ftable
!f2py intent(out) :: tempflux

!! Read in all the support files
open(unit=747,file="tmp/atmkeys")
open(unit=748,file="tmp/tandg")
open(unit=749,file="tmp/cossquiggle_i"//trim(ext))
open(unit=750,file="tmp/doppler_arr_i"//trim(ext))
do i=0,nzth-1
    read(747,*) (atmkeys(i,j),j=0,3)
    read(748,*) (tandg(i,j),j=0,5)
enddo
do i=0,nzphi-1
    read(749,*) (cossquiggle(i,j),j=0,nzth-1)
enddo
do i=0,nzphi-1
    read(750,*) (doppler_arr(i,j),j=0,nzth-1)
enddo
close(747)
close(748)
close(749)
close(750)

!resolution:
res = wav(1)-wav(0)

iout(:,:,:)=0.

do i=0,nzphi-1
    do j=0,nzth-1
        theta1=0
        theta2=1
        do k=0,ntheta-2
        if(cossquiggle(i,j)>=angles(k)) then
            theta1 = k
            theta2 = k+1
        endif
        enddo
        ainte = (cossquiggle(i,j)-angles(theta1))/(angles(theta2)-angles(theta1))
        if(cossquiggle(i,j)>=0) then
        do k=0,nwav-1
	!apply doppler shift:
            deltaw = doppler_arr(i,j)*wav(k)
            wavn = wav(k)-deltaw
            n = NINT((wavn-wav(0))/res)
	!found n to put the intensity
		if(n>=0 .and. n < nwav) then
	            p1 = (1.-ainte)*ftable(k,theta1,atmkeys(j,0))+ainte*ftable(k,theta2,atmkeys(j,0))
	            p2 = (1.-ainte)*ftable(k,theta1,atmkeys(j,1))+ainte*ftable(k,theta2,atmkeys(j,1))
	            p3 = (1.-ainte)*ftable(k,theta1,atmkeys(j,2))+ainte*ftable(k,theta2,atmkeys(j,2))
	            p4 = (1.-ainte)*ftable(k,theta1,atmkeys(j,3))+ainte*ftable(k,theta2,atmkeys(j,3))
	            temp_iout = tandg(j,5)*(tandg(j,4)*p1+(1.-tandg(j,4))*p2)+ &
	                        (1-tandg(j,5))*(tandg(j,4)*p3+(1.-tandg(j,4))*p4)
	            temp_iout = 10**(temp_iout)
	            iout(n,i,j) = iout(n,i,j)+temp_iout
		endif
        enddo
        endif
    enddo
enddo

call sumup(iout,tempflux,wav,dist,ext,ntheta,nwav,nzth,nzphi,incl)

return
end subroutine ! interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!1!!1!!!!!1!!!!!!!!!!!!!!!!!
subroutine interp_nodopp(tempflux,ftable,wav,angles,ext,dist,nfnames,ntheta,nwav,nzth,nzphi,incl)
implicit none
integer (kind=4) :: i,j,k,l,m,theta1,theta2
integer (kind=4) :: nwav,ntheta,nfnames,nzth,nzphi
real (kind=8) :: incl,p1,p2,p3,p4,ainte,dist,f11,f12,f21,f22,f31,f32,f41,f42
real (kind=8), dimension(0:nwav-1), intent(in) :: wav
!f2py depend(nwav) wav
real (kind=8), dimension(0:nwav-1,0:ntheta-1,0:nfnames-1), intent(in) :: ftable
!f2py depend(nwav,ntheta,nfnames) ftable
real (kind=8), dimension(0:nwav-1,0:nzphi-1,0:nzth-1) :: iout
!f2py depend(nwav,nzphi,nzth) iout
real (kind=8), dimension(0:nwav-1,0:1), intent(out) :: tempflux
!f2py depend(nwav) tempflux
real (kind=8), dimension(0:ntheta-1) :: angles
!f2py depend(ntheta) angles
real (kind=8), dimension(0:nzphi-1,0:nzth-1) :: cossquiggle
!f2py depend(nzphi,nzth) cossquiggle
integer (kind=4), dimension(0:nzth-1,0:3) :: atmkeys
!f2py depend(nzth) atmkeys
real (kind=8), dimension(0:nzth-1,0:5) :: tandg
!f2py depend(nzth) tandg
character (len=20) :: ext
!f2py intent(in) :: wav,ftable
!f2py intent(out) :: tempflux

!! Read in all the support files
open(unit=747,file="tmp/atmkeys")
open(unit=748,file="tmp/tandg")
open(unit=749,file="tmp/cossquiggle_i"//trim(ext))
do i=0,nzth-1
    read(747,*) (atmkeys(i,j),j=0,3)
    read(748,*) (tandg(i,j),j=0,5)
enddo
do i=0,nzphi-1
    read(749,*) (cossquiggle(i,j),j=0,nzth-1)
enddo
close(747)
close(748)
close(749)

do i=0,nzphi-1
    do j=0,nzth-1
        theta1=0
        theta2=1
        do k=0,ntheta-2
        if(cossquiggle(i,j)>=angles(k)) then
            theta1 = k
            theta2 = k+1
        endif
        enddo
        ainte = (cossquiggle(i,j)-angles(theta1))/(angles(theta2)-angles(theta1))
        if(cossquiggle(i,j)>=0) then
        do k=0,nwav-1
            p1 = (1.-ainte)*ftable(k,theta1,atmkeys(j,0))+ainte*ftable(k,theta2,atmkeys(j,0))
            p2 = (1.-ainte)*ftable(k,theta1,atmkeys(j,1))+ainte*ftable(k,theta2,atmkeys(j,1))
            p3 = (1.-ainte)*ftable(k,theta1,atmkeys(j,2))+ainte*ftable(k,theta2,atmkeys(j,2))
            p4 = (1.-ainte)*ftable(k,theta1,atmkeys(j,3))+ainte*ftable(k,theta2,atmkeys(j,3))
            iout(k,i,j)=tandg(j,5)*(tandg(j,4)*p1+(1.-tandg(j,4))*p2)+ &
                        (1-tandg(j,5))*(tandg(j,4)*p3+(1.-tandg(j,4))*p4)
			iout(k,i,j) = 10**(iout(k,i,j))
        enddo
        else
            iout(:,i,j)=0.
        endif
    enddo
enddo

call sumup(iout,tempflux,wav,dist,ext,ntheta,nwav,nzth,nzphi,incl)

return
end subroutine ! interp_nodopp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sumup(iout,tempflux,wav,dist,ext,ntheta,nwav,nzth,nzphi,incl)
implicit none
integer (kind=4) :: i,j,k,l,m
integer (kind=4) :: nwav,ntheta,nfnames,nzth,nzphi
real (kind=8) :: incl,d,dist,temp
real (kind=8), dimension(0:nwav-1), intent(in) :: wav
!f2py depend(nwav) wav
real (kind=8), dimension(0:nwav-1,0:1), intent(inout) :: tempflux
!f2py depend(nwav) tempflux
real (kind=8), dimension(0:nwav-1,0:nzphi-1,0:nzth-1), intent(in) :: iout
!f2py depend(nwav,nzphi,nzth) iout
real (kind=8), dimension(0:nzphi-1,0:nzth-1) :: cossquiggle,darea
!f2py depend(nzphi,nzth) cossquiggle, darea
character (len=20) :: ext
!f2py intent(in,out) tempflux
!f2py intent(in) wav,fout

d = dist

!! Read in all the support files
open(unit=749,file="tmp/cossquiggle_i"//trim(ext))
open(unit=750,file="tmp/darea")
do i=0,nzphi-1
    read(749,*) (cossquiggle(i,j),j=0,nzth-1)
    read(750,*) (darea(i,j),j=0,nzth-1)
enddo
close(749)
close(750)

do i=0,nwav-1
    temp = sum(iout(i,:,:)*cossquiggle(:,:)*darea(:,:),MASK=iout(i,:,:).GT.0.0)/((d**2)*1.0e8)
    tempflux(i,1) = temp
    tempflux(i,0) = wav(i)
enddo

return
end subroutine ! sumup
