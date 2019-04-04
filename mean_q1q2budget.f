      PARAMETER (im=202, km=52, inn=480,itt=480*73,ii=365*4,ifi=73)
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
      dimension con(im,km),eva(im,km),dep(im,km)     !!readin 
      dimension sub(im,km),fus(im,km)                    !!readin 
      dimension smcon(itt,km),smeva(itt,km),smdep(itt,km)     !
      dimension smsub(itt,km),smfus(itt,km)                    
      dimension tmcon(ii,km),tmeva(ii,km),tmdep(ii,km)     !! output
      dimension tmsub(ii,km),tmfus(ii,km)                    !! output
      character*100 chenm, path ! case name ,must be recorded 
      character casenm*20,fold*20,strid*3
      casenm='ETP2D0'
      fold='run2'
      if(casenm(1:3)=='MLY')then
        path='/public/home/chenjh/Models/CRM/ERA/Year/'//casenm(1:4)//
     +   '/'//trim(fold)
      else
        path='/public/home/chenjh/Models/CRM/ERA/Year/'//casenm(1:3)//
     +   '/'//trim(fold)
      endif
!     path='/home/jhchen/jhchen/ERA_Interim/ETP/'//trim(fold)!
       chenm=trim(path)//'/eddydiffradcon_'//trim(casenm)    !'whereisthe_data_tgn2d1_(1-6)'      
!       open(81,file=trim(chenm)//'_1'
!     * ,form='UNFORMATTED',status='OLD',convert='big_endian') 
!      open(82,file=trim(chenm)//'_2'
!     *,form='unformatted',status='old',convert='big_endian') 
!      open(83,file=trim(chenm)//'_3'
!     *,form='unformatted',status='old',convert='big_endian') 
!      open(84,file=trim(chenm)//'_4'
!     *,form='unformatted',status='old',convert='big_endian') 
!      open(85,file=trim(chenm)//'_5'
!     *,form='unformatted',status='old',convert='big_endian') 
!      open(86,file=trim(chenm)//'_6'
!     *,form='unformatted',status='old',convert='big_endian') 
      open(40,file=trim(chenm)//'q1q2budget_Allcomps'
     * ,form='unformatted')
      open(888,file=trim(chenm)//'_ALLcomps.txt')
c      open(51,file=trim(chenm)//'_rho.txt')
c      open(52,file=trim(chenm)//'_Raw_q12edfx.txt')
c      open(53,file=trim(chenm)//'_Raw_q12seflux.txt')
c      open(54,file=trim(chenm)//'_Raw_rlwrsw.txt')
      sec=10.           !!!Time Step
      it=0
      do 999 if=1,ifi
      IH=80+if
      write(strid,'(I3.3)')if
      chenm=trim(path)//'/eddydiffradcon_'//trim(casenm)//'_'//strid
       open(IH,file=trim(chenm)
     * ,form='UNFORMATTED',status='OLD',convert='big_endian')
      do 100 in=1,inn
      it=it+1
      read(IH) q1e,q2e,fx 
      read(IH) q1d,q2d,q1s,q2s
      read(IH) rlw,rsw
      read(IH) e1flux,e2flux,eflux
      read(IH) con
      read(IH) eva
      read(IH) dep
      read(IH) sub
      read(IH) fus

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
      smcon(it,k)=0.
      smeva(it,k)=0.
      smdep(it,k)=0.
      smsub(it,k)=0.
      smfus(it,k)=0.
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
      smcon(it,k)=smcon(it,k)+con(i,k)*3600.*24./float(im-2)/sec
      smeva(it,k)=smeva(it,k)+eva(i,k)*3600.*24./float(im-2)/sec
      smdep(it,k)=smdep(it,k)+dep(i,k)*3600.*24./float(im-2)/sec
      smsub(it,k)=smsub(it,k)+sub(i,k)*3600.*24./float(im-2)/sec
      smfus(it,k)=smfus(it,k)+fus(i,k)*3600.*24./float(im-2)/sec
  30  continue
  40  continue

 100  continue 
      close(IH)  
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
      tmcon(mt,k)=0.
      tmeva(mt,k)=0.
      tmdep(mt,k)=0.
      tmsub(mt,k)=0.
      tmfus(mt,k)=0.      
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
      tmcon(mt,k)=tmcon(mt,k)+smcon(n,k)/float(n2-n1+1)
      tmeva(mt,k)=tmeva(mt,k)+smeva(n,k)/float(n2-n1+1)
      tmdep(mt,k)=tmdep(mt,k)+smdep(n,k)/float(n2-n1+1)
      tmsub(mt,k)=tmsub(mt,k)+smsub(n,k)/float(n2-n1+1)
      tmfus(mt,k)=tmfus(mt,k)+smfus(n,k)/float(n2-n1+1)
  11  continue
  50  continue
 200  continue
      write(40) tmq1e,tmq2e,tmq1d,tmq2d
      write(40) tmq1s,tmq2s,tmrlw,tmrsw
      write(40) tmfx
      write(40) tme1flux,tme2flux,tmeflux
      write(40) tmcon,tmeva,tmdep
      write(40) tmsub,tmfus
      print*,'tmq1e(1,1)= , tmq2e(1,1)= ',tmq1e(1,1),tmq2e(1,1)
      print*,'tmq1e(80,1)= , tmq2e(80,1)= ',tmq1e(80,1),tmq2e(80,1)
c
      do mt=1,ii
       write(888,99) tmq1e(mt,:),tmq2e(mt,:)
     +     ,  tmq1d(mt,:),tmq2d(mt,:)
     + , tmq1s(mt,:), tmq2s(mt,:),tmrlw(mt,:),tmrsw(mt,:)
     + , tmfx(mt,:), tme1flux(mt,:),tme2flux(mt,:),tmeflux(mt,:)
     + , tmcon(mt,:),tmeva(mt,:),tmdep(mt,:),tmsub(mt,:),tmfus(mt,:)
      enddo
99    format(8e12.4)
c  units:              
c     all the variables are in K/day
      stop
      end
