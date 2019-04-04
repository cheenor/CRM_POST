      PARAMETER (im=202, km=52, inn=480, itt=2880, ii=121, ifi=6)
c      dimension qc(im,km),qr(im,km),qa(im,km),qb(im,km),rho(km)
      dimension q1e(im,km),q2e(im,km),fx(im,km)
      dimension q1d(im,km),q2d(im,km)
      dimension q1s(im,km),q2s(im,km),rlw(im,km),rsw(im,km)
      dimension e1flux(im,km),e2flux(im,km),eflux(im,km)
c      dimension u(im,km),w(im,km),t(im,km),q(im,km),rh(im,km)
      dimension smq1e(itt,km),smq2e(itt,km),smfx(itt,km),smq1d(itt,km)
      dimension smq2d(itt,km),smq1s(itt,km),smq2s(itt,km)
      dimension sme1flux(itt,km),sme2flux(itt,km),smeflux(itt,km)
      dimension smrlw(itt,km),smrsw(itt,km),tmrlw(ii,km),tmrsw(ii,km)
      dimension tmq1e(ii,km),tmq2e(ii,km),tmfx(ii,km),tmq1d(ii,km)
      dimension tmq2d(ii,km),tmq1s(ii,km),tmq2s(ii,km)    ! ,tmq(ii,km)
      dimension tme1flux(ii,km),tme2flux(ii,km),tmeflux(ii,km) ! ,tmftc(ii,km)
      character*100 chenm, path ! case name ,must be recorded 
      character casenm*20,runfold*20
      casenm='casename'
      runfold='runfold'
      if(casenm(1:3)=='MLY')then
        path='/home/jhchen/jhchen/ERA_Interim/'//casenm(1:4)//
     +   '/'//trim(runfold)
      else
        path='/home/jhchen/jhchen/ERA_Interim/'//casenm(1:3)//
     +   '/'//trim(runfold)    
      endif
!      path='/home/jhchen/jhchen/ERA_Interim/ETP/'//trim(runfold)
      chenm=trim(path)//'/eddydiffrad_'//trim(casenm)    !'whereisthe_data_tgn2d1_(1-6)'      
       open(81,file=trim(chenm)//'_1'
     * ,form='UNFORMATTED',status='OLD',convert='big_endian') 
      open(82,file=trim(chenm)//'_2'
     *,form='unformatted',status='old',convert='big_endian') 
      open(83,file=trim(chenm)//'_3'
     *,form='unformatted',status='old',convert='big_endian') 
      open(84,file=trim(chenm)//'_4'
     *,form='unformatted',status='old',convert='big_endian') 
      open(85,file=trim(chenm)//'_5'
     *,form='unformatted',status='old',convert='big_endian') 
      open(86,file=trim(chenm)//'_6'
     *,form='unformatted',status='old',convert='big_endian') 
      open(40,file=trim(chenm)//'_All'
     * ,form='unformatted')
      open(50,file=trim(chenm)//'_All.txt')
c      open(51,file=trim(chenm)//'_rho.txt')
c      open(52,file=trim(chenm)//'_Raw_q12edfx.txt')
c      open(53,file=trim(chenm)//'_Raw_q12seflux.txt')
c      open(54,file=trim(chenm)//'_Raw_rlwrsw.txt')
      it=0
      do 999 if=1,ifi
      IH=80+if
      do 100 in=1,inn
      it=it+1
      read(IH) q1e,q2e,fx 
      read(IH) q1d,q2d,q1s,q2s
      read(IH) rlw,rsw
      read(IH) e1flux,e2flux,eflux
C
CC  the following loop output the data of every time stemp and grid. 
c      do i=2,im-1
c      write(52,99)q1e(i,:),q2e(i,:),fx(i,:)
c     +            ,q1d(i,:),q2d(i,:)
c      write(53,99)q1s(i,:),q2s(i,:)
c     +            ,e1flux(i,:),e2flux(i,:)
c     +            ,eflux(i,:)
c      write(54,99)rlw(i,:),rsw(i,:) 
c      end do
C
      do 40 k=1,km
c      write(51,*)k,rho(k)
      smq1e(it,k)=0.
      smq2e(it,k)=0.
      smfx(it,k)=0.
      smq1d(it,k)=0.
      smq2d(it,k)=0.
      smq1s(it,k)=0.
      smq2s(it,k)=0.
      sme1flux(it,k)=0.
      sme2flux(it,k)=0.
      smeflux(it,k)=0.
      smrsw(it,k)=0.
      smrlw(it,k)=0.
      do 30 i=2,im-1
      smq1e(it,k)=smq1e(it,k)+q1e(i,k)/float(im-2) 
      smq2e(it,k)=smq2e(it,k)+q2e(i,k)/float(im-2) 
      smfx(it,k)=smfx(it,k)+fx(i,k)/float(im-2) 
      smq1d(it,k)=smq1d(it,k)+q1d(i,k)/float(im-2) 
      smq2d(it,k)=smq2d(it,k)+q2d(i,k)/float(im-2) 
      smq1s(it,k)=smq1s(it,k)+q1s(i,k)/float(im-2) 
      smq2s(it,k)=smq2s(it,k)+q2s(i,k)/float(im-2) 
      sme1flux(it,k)=sme1flux(it,k)+e1flux(i,k)/float(im-2) 
      sme2flux(it,k)=sme2flux(it,k)+e2flux(i,k)/float(im-2) 
      smeflux(it,k)=smeflux(it,k)+eflux(i,k)/float(im-2)
      smrlw(it,k)=smrlw(it,k)+rlw(i,k)/float(im-2)
      smrsw(it,k)=smrsw(it,k)+rsw(i,k)/float(im-2)
  30  continue
  40  continue

 100  continue   
 999  continue   

      do 200 k=1,km
      do 50 mt=1,ii
      n1=(mt-1)*24-11
      n2=(mt-1)*24+12
      if(mt.eq.1) then
      n1=1
      n2=12
      endif
      if(mt.eq.ii) then
      n1=itt-11
      n2=itt
      endif
      tmq1e(mt,k)=0.
      tmq2e(mt,k)=0.
      tmfx(mt,k)=0.
      tmq1d(mt,k)=0.
      tmq2d(mt,k)=0.
      tmq1s(mt,k)=0.
      tmq2s(mt,k)=0.
      tme1flux(mt,k)=0.
      tme2flux(mt,k)=0.
      tmeflux(mt,k)=0.
      tmrlw(mt,k)=0.
      tmrsw(mt,k)=0.
      do 11 n=n1,n2 
      tmq1e(mt,k)=tmq1e(mt,k)+smq1e(n,k)/float(n2-n1+1)
      tmq2e(mt,k)=tmq2e(mt,k)+smq2e(n,k)/float(n2-n1+1)
      tmq1d(mt,k)=tmq1d(mt,k)+smq1d(n,k)/float(n2-n1+1)
      tmq2d(mt,k)=tmq2d(mt,k)+smq2d(n,k)/float(n2-n1+1)
      tmq1s(mt,k)=tmq1s(mt,k)+smq1s(n,k)/float(n2-n1+1)
      tmq2s(mt,k)=tmq2s(mt,k)+smq2s(n,k)/float(n2-n1+1)
      tmfx(mt,k)=tmfx(mt,k)+smfx(n,k)/float(n2-n1+1)
      tme1flux(mt,k)=tme1flux(mt,k)+sme1flux(n,k)/float(n2-n1+1)
      tme2flux(mt,k)=tme2flux(mt,k)+sme2flux(n,k)/float(n2-n1+1)
      tmrlw(mt,k)=tmrlw(mt,k)+smrlw(n,k)/float(n2-n1+1)
      tmrsw(mt,k)=tmrsw(mt,k)+smrsw(n,k)/float(n2-n1+1)
  11  continue
  50  continue
 200  continue
      write(40) tmq1e,tmq2e,tmq1d,tmq2d
      write(40) tmq1s,tmq2s,tmrlw,tmrsw
      write(40) tmfx
      write(40) tme1flux,tme2flux,tmeflux
      print*,'tmq1e(1,1)= , tmq2e(1,1)= ',tmq1e(1,1),tmq2e(1,1)
      print*,'tmq1e(80,1)= , tmq2e(80,1)= ',tmq1e(80,1),tmq2e(80,1)
c
      do mt=1,ii
      write(50,99) tmq1e(mt,:),tmq2e(mt,:),tmq1d(mt,:),tmq2d(mt,:)
     + ,tmq1s(mt,:),
     +tmq2s(mt,:),tmrlw(mt,:),tmrsw(mt,:), tmfx(mt,:),
     +tme1flux(mt,:),tme2flux(mt,:),tmeflux(mt,:) 
      enddo
99    format(8e12.4)
c  units:              
c     horizontal velocity (tmu  m/s)
c     vertical velocity (tmw  m/s)
c     temperature (tmt  K)
c     water vapor mixing ratio (tmq  g/kg)
c     relative humidity (tmrh  %)
c     cloud water mixing ratio (tmqc  g/kg)
c     rain water mixing ratio  (tmqr  g/kg)
c     type A ice water mixing ratio (tmqa  g/kg)
c     type B ice water mixing ratio (tmqb  g/kg)
      stop
      end
