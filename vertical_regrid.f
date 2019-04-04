      PARAMETER (in=365*4, kn=33)
      PARAMETER (im=365*4, km=52, NZM=km-1)
      PARAMETER (it=365*4)
      PARAMETER (nx=im,nz=km-1,nxg=nx,nzg=33)
      dimension XI(km)
      dimension uo(in,kn),to(in,kn),qo(in,kn)
      dimension umi(im,km),tmi(im,km),emi(im,km),qmi(im,km)
      dimension q1e(im,km),q2e(im,km),dmi(im,km)
      dimension q1d(im,km),q2d(im,km)
      dimension q1b(im,km),q2b(im,km)
      dimension rlw(im,km),rsw(im,km)
      dimension aq1m(nx,nz),aq2m(nx,nz)
      dimension tm(nxg,nzg),qm(nxg,nzg)
      dimension bto(nzg),bqo(nzg),btm(nzg),bqm(nzg)
      dimension btqm(nzg)
      dimension x(nx),z(nz),zm(nz),xg(nxg),zg(nzg)
      dimension work1(nx,nz),work2(nx,nz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension eddydiffrad(12,it,km)
!      dimension tmq1e(it,km),tmq2e(it,km),tmq1d(it,km),tmq2d(it,km)
!      dimension tmq1s(it,km),tmq2s(it,km),tmrlw(it,km),tmrsw(it,km),tmfx(it,km)
!      dimension tme1flux(it,km),tme2flux(it,km),tmeflux(it,km)      
      dimension qabcr(13,it,km),tmpwork(nx,nz),dat43(6,it,km)
!      tmt(mt,:),tmq(mt,:),tmu(mt,:),tmw(mt,:),tmqc(mt,:),
!     +tmqr(mt,:),tmqa(mt,:),tmqb(mt,:), tmrh(mt,:),
!     +tmfqv(mt,:),tmftc(mt,:),tmrat(mt,:)
      dimension bugetqc(5,it,km)
!     tmc(it,:),tme(it,:),tmd(it,:),tms(it,:),tmf(it,:)
      dimension obsq12(2,it,km),obsftq(2,it,km)
      dimension simq12(2,it,km)
!     q1    q2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      character path*100,fold*20,casenm*30,obsq12nm*20
      
      casenm='ETP2D0'
      obsq12nm='ETP_20100101_365d'
!!!!!!!!!!! file 1
      fold='run2'
      if(casenm(1:3)=='MLY')then
       path='/public/home/chenjh/Models/CRM/ERA/Year/'//casenm(1:4)//'/'
      else
       path='/public/home/chenjh/Models/CRM/ERA/Year/'//casenm(1:3)//'/'
      endif      
!      open(30,file=trim(path)//trim(fold)//'/eddydiffradcon_'//
!     *       trim(casenm)//'_ALLcomps.txt')
      open(30,file=trim(path)//trim(fold)//'/eddydiffradcon_'
     +  //trim(casenm)//'_ALLcomps.txt') !eddydiffradcon_MLYR2D0q1q2budget_Allcomps
      open(31,file=trim(path)//trim(fold)//'/regrid/'//trim(casenm)//
     * '_regrid_eddydiffrad_All.txt' )
      open(40,file=trim(path)//trim(fold)//'/'//trim(casenm)//
     * '_All.txt' )
      open(41,file=trim(path)//trim(fold)//'/regrid/'//trim(casenm)//
     * '_regrid_qabcr_All.txt' )
!      open(50,file=trim(path)//trim(fold)//'/postdata/'//trim(casenm)//
!     * '_micro_202_6hour.txt' )
      open(51,file=trim(path)//trim(fold)//'/regrid/'//trim(casenm)//
     * '_micro_regrid_6hour.txt' )
      open(60,file=trim(path)//trim(fold)//'/'//trim(obsq12nm)//
     +     '_ERA.99')
      open(61,file=trim(path)//trim(fold)//'/regrid/'//trim(obsq12nm)
     +       //'_regrid_ERA.99')
      open(62,file=trim(path)//trim(fold)//'/regrid/'//trim(casenm)
     +       //'_q1q2_regrid.txt')
      open(70,file=trim(path)//trim(fold)//'/regrid/'//
     +            '/gridded_height_33levels.txt')
      open(80,file=trim(path)//trim(fold)//'/'//trim(obsq12nm)//
     +   '_LSFORCING_ERA.37')
      open(81,file=trim(path)//trim(fold)//'/regrid/'//trim(obsq12nm)
     +       //'_lsforcing_regrid_ERA.37')
c 
      open(82,file=trim(path)//trim(fold)//'/'//trim(obsq12nm)//
     +   '_ERA.43')
      open(83,file=trim(path)//trim(fold)//'/regrid/'//trim(obsq12nm)
     +       //'_regrid_ERA.43')
c     open(31,file='/hiba/xiaoqing/gate/tqpcedsf_cs13.dat'
c    *,form='unformatted',status='old')
c     open(32,file='/hiba/xiaoqing/gate/eddydiffrad_cs13.dat'
c    *,form='unformatted',status='old')
c      open(31,file='/hiba/xiaoqing/gate_dat/tqpcesf_450_os1.dat'
c     *,form='unformatted',status='old')
c      open(32,file='/hiba/xiaoqing/gate_dat/eddydiffrad_450_os1.dat'
c     *,form='unformatted',status='old')
c     open(42,file=outfile,status='new')
c
c      read(30) to,qo,uo    
c      read(31)     
c      read(31)  
cc     read(31) tmi,emi,dmi,qmi,umi 
c      read(31) tmi,emi,qmi,umi 
c      read(32) q1e,q2e 
c      read(32) q1d,q2d 
c      read(32) rlw,rsw 
c      read(32) q1b,q2b 
c      read(32) 
c     write(42,23) ((tmi(i,k),i=1,57),k=1,62)
c
cc  setup with 100m increasing to 500 m at 6 km, 500 m to 26 km
cc  1 domain
      rat=15.
      XI(2)=0.
      DO 162 K=2,NZM
      RATZ=RAT
      DEL=100.
      nzm1=nz
      k1=k
      XI(K+1)=XI(K)+((RATZ-1.)/FLOAT(NZM1-2)*FLOAT(K1-2)+1.)*DEL
 162  CONTINUE
      do 153 k=1,nz
 153  z(k)=XI(K+1)/1000.
      zm(1)=0.
      do 154 k=2,nz
 154  zm(k)=0.5*(z(k-1)+z(k))
      do 150 i=1,nx
 150  x(i)=float(i-1)
      do 252 k=1,nzg
 252  zg(k)=float(k-1)/2.
      do 250 i=1,nxg
 250  xg(i)=float(i-1)
      do  k=1,nzg
      write(70,79)zg(k)
      enddo
C      print*,zm
C      print*,zg
C      print*,x,xg
79    format(1X,F12.4)
c      do 10 i=1,nx
c      rlw(i,1)=0.
c      rsw(i,1)=0.
c      do 10 k=1,nz
c     aq1m(i,k)=tmi(i,k)+emi(i,k)+dmi(i,k)+qmi(i,k)+umi(i,k)
c      aq1m(i,k)=tmi(i,k)+emi(i,k)+qmi(i,k)+umi(i,k)
c     *         +rlw(i,k)+rsw(i,k)+q1e(i,k)+q1d(i,k)
c     aq2m(i,k)=tmi(i,k)+emi(i,k)+dmi(i,k)+qmi(i,k)*2.5e10/2.834e10
c      aq2m(i,k)=tmi(i,k)+emi(i,k)+qmi(i,k)*2.5e10/2.834e10
c     *         +q2e(i,k)+q2d(i,k)
c      aq1m(i,k)=tmi(i,k)+emi(i,k)+dmi(i,k)+qmi(i,k)+umi(i,k)
c      aq2m(i,k)=tmi(i,k)+emi(i,k)+dmi(i,k)+qmi(i,k)*2.5e10/2.834e10
c     aq1m(i,k)=rlw(i,k)
c     aq2m(i,k)=rsw(i,k)
c     aq1m(i,k)=q1e(i,k)
c     aq2m(i,k)=q2e(i,k)
c     aq1m(i,k)=q1d(i,k)
c     aq2m(i,k)=q2d(i,k)
c  10  continue
c     
cc interpolate into regular grid:
      do mt=1,it
       read(30,99)eddydiffrad(1,mt,:),eddydiffrad(2,mt,:),
     +                              eddydiffrad(3,mt,:),
     +      eddydiffrad(4,mt,:),eddydiffrad(5,mt,:),eddydiffrad(6,mt,:),
     +      eddydiffrad(7,mt,:),eddydiffrad(8,mt,:),eddydiffrad(9,mt,:),
     +    eddydiffrad(10,mt,:),eddydiffrad(11,mt,:),eddydiffrad(12,mt,:)
     +    , bugetqc(1,mt,:),bugetqc(2,mt,:),bugetqc(3,mt,:),    !!!!bugetqc
     +      bugetqc(4,mt,:),bugetqc(5,mt,:)
       read(40,99)qabcr(1,mt,:),qabcr(2,mt,:),qabcr(3,mt,:),
     +      qabcr(4,mt,:),qabcr(5,mt,:),qabcr(6,mt,:),
     +      qabcr(7,mt,:),qabcr(8,mt,:),qabcr(9,mt,:),
     +      qabcr(10,mt,:),qabcr(11,mt,:),qabcr(12,mt,:)
     +     , qabcr(13,mt,:)
!      read(50,99)bugetqc(1,mt,:),bugetqc(2,mt,:),bugetqc(3,mt,:),    !!!!bugetqc
!     +      bugetqc(4,mt,:),bugetqc(5,mt,:)
      enddo
      do i=1,nx
         simq12(1,i,1)=0.   !rlw  long wave
         simq12(2,i,1)=0.   !rsw  short wave
         do k=2,nz
           simq12(1,i,k)=bugetqc(1,i,k)+bugetqc(2,i,k)+bugetqc(3,i,k)
     +       + bugetqc(4,i,K)+bugetqc(5,i,k)
     +       +eddydiffrad(8,i,k)+eddydiffrad(7,i,k)+eddydiffrad(1,i,k) ! qr(long+short) qe
     +       +eddydiffrad(3,i,k)+eddydiffrad(5,i,k)              !q1d q1s
           simq12(2,i,k)=bugetqc(1,i,K)+bugetqc(2,i,k)
     +         + bugetqc(3,i,K)+bugetqc(4,i,k)*2.5e10/2.834e10
     +         +eddydiffrad(2,i,k)+eddydiffrad(4,i,k)          !q2e q2d
     +         +eddydiffrad(6,i,k)                             ! q2s
         enddo
      enddo
!############################
      do iv=1,12
!      tmpwork(:,:)=eddydiffrad(iv,:,:)
       do mt=1,it
          tmpwork(mt,1)=0. 
          do k=2,nz
           tmpwork(mt,k)=eddydiffrad(iv,mt,k)
          enddo
        enddo
!      tmpwork(:,1)=0.0   !!!!! the bottom level is zero
       call gridint(tmpwork,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
      do  k=1,nzg
        do  i=1,nxg
          eddydiffrad(iv,i,k)=tmpwork(i,k)
c     to(i,k)=to(i,k)-uo(i,k)
        enddo
       enddo
      enddo  
      do mt=1,it
C       write(31,99)eddydiffrad(1,mt,:),eddydiffrad(2,mt,:),
C     +                              eddydiffrad(3,mt,:),
c     +      eddydiffrad(4,mt,:),eddydiffrad(5,mt,:),eddydiffrad(6,mt,:),
C     +      eddydiffrad(7,mt,:),eddydiffrad(8,mt,:),eddydiffrad(9,mt,:),
C     +    eddydiffrad(10,mt,:),eddydiffrad(11,mt,:),eddydiffrad(12,mt,:)
      write(31,99)(eddydiffrad(1,mt,k),k=1,nzg),
     +                 (eddydiffrad(2,mt,k),k=1,nzg)
     +    ,(eddydiffrad(3,mt,k),k=1,nzg),(eddydiffrad(4,mt,k),k=1,nzg)
     +    ,(eddydiffrad(5,mt,k),k=1,nzg),(eddydiffrad(6,mt,k),k=1,nzg)
     +    ,(eddydiffrad(7,mt,k),k=1,nzg),(eddydiffrad(8,mt,k),k=1,nzg)
     +    ,(eddydiffrad(9,mt,k),k=1,nzg),(eddydiffrad(10,mt,k),k=1,nzg)
     +   ,(eddydiffrad(11,mt,k),k=1,nzg),(eddydiffrad(12,mt,k),k=1,nzg)
!      read(40,99)qabcr(1,mt,:),qabcr(2,mt,:),qabcr(3,mt,:),
!     +      qabcr(4,mt,:),qabcr(5,mt,:),qabcr(6,mt,:),
!     +      qabcr(7,mt,:),qabcr(8,mt,:),qabcr(9,mt,:),
!     +      qabcr(10,mt,:),qabcr(11,mt,:),qabcr(12,mt,:)
!     +     , qabcr(13,mt,:)
       enddo
!!!!!!!!!!!!!!!!!!!!!!!! file 40
      do iv=1,13
       do mt=1,it
!          tmpwork(mt,1)=0
          do k=1,nz
            tmpwork(mt,k)=qabcr(iv,mt,k)
          enddo
        enddo
!       tmpwork(:,:)=qabcr(iv,:,:)
!       tmpwork(:,1)=0.0   !!!!! the bottom level is zero
       call gridint(tmpwork,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
       do  k=1,nzg
         do  i=1,nxg
          qabcr(iv,i,k)=tmpwork(i,k)
c     to(i,k)=to(i,k)-uo(i,k)
        enddo
        enddo
       enddo 
      do mt=1,it
      write(41,99)(qabcr(1,mt,k),k=1,nzg),(qabcr(2,mt,k),k=1,nzg)
     +              ,(qabcr(3,mt,k),k=1,nzg),(qabcr(4,mt,k),k=1,nzg)
     +              ,(qabcr(5,mt,k),k=1,nzg),(qabcr(6,mt,k),k=1,nzg)
     +              ,(qabcr(7,mt,k),k=1,nzg),(qabcr(8,mt,k),k=1,nzg)
     +              ,(qabcr(9,mt,k),k=1,nzg),(qabcr(10,mt,k),k=1,nzg)
     +              ,(qabcr(11,mt,k),k=1,nzg),(qabcr(12,mt,k),k=1,nzg)
     +              ,(qabcr(13,mt,k),k=1,nzg)
!      read(50,99)bugetqc(1,mt,:),bugetqc(2,mt,:),bugetqc(3,mt,:),    !!!!bugetqc
!     +      bugetqc(4,mt,:),bugetqc(5,mt,:)
      enddo
!  !!!!!!!!!!!!!!!!!!!!!!!! file 50
!      do i=1,nx      
!         eddydiffrad(7,i,1)=0.   !rlw  long wave
!         eddydiffrad(8,i,1)=0.   !rsw  short wave
!         do k=1,nz
!           simq12(1,i,k)=bugetqc(1,i,k)+bugetqc(2,mt,k)+bugetqc(3,mt,k)
!     +       + bugetqc(4,mt,K)+bugetqc(5,mt,k)
!     +       +eddydiffrad(8,i,k)+eddydiffrad(7,i,k)+eddydiffrad(1,i,k) ! qr(long+short) qe
!     +       +eddydiffrad(3,i,k)+eddydiffrad(5,i,k)              !q1d q1s
!           simq12(2,i,k)=bugetqc(1,mt,K)+bugetqc(2,mt,k)
!     +         + bugetqc(3,mt,K)+bugetqc(4,mt,k)*2.5e10/2.834e10
!     +         +eddydiffrad(2,i,k)+eddydiffrad(4,i,k)          !q2e q2d
!     +         +eddydiffrad(6,i,k)                             ! q2s
!         enddo
!      enddo
      do iv=1,2
!      tmpwork(:,:)=simq12(iv,:,:)
       do mt=1,it
          do k=1,nz
            tmpwork(mt,k)=simq12(iv,mt,k)
          enddo
        enddo
!      tmpwork(:,1)=0.0   !!!!! the bottom level is zero
       call gridint(tmpwork,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
        do  k=1,nzg
          do  i=1,nxg
            simq12(iv,i,k)=tmpwork(i,k)
c     to(i,k)=to(i,k)-uo(i,k)
           enddo
         enddo
       enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      do iv=1,5
!      tmpwork(:,:)=bugetqc(iv,:,:)
       do mt=1,it
          tmpwork(mt,1)=0.
          do k=2,nz
            tmpwork(mt,k)=bugetqc(iv,mt,k)
          enddo
        enddo
!      tmpwork(:,1)=0.0   !!!!! the bottom level is zero
      call gridint(tmpwork,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
        do  k=1,nzg
          do  i=1,nxg
            bugetqc(iv,i,k)=tmpwork(i,k)
c     to(i,k)=to(i,k)-uo(i,k)
           enddo
         enddo
      enddo 
      do mt=1,it
       write(51,99)(bugetqc(1,mt,k),k=1,nzg),(bugetqc(2,mt,k),k=1,nzg)
     +              ,(bugetqc(3,mt,k),k=1,nzg),(bugetqc(4,mt,k),k=1,nzg)
     +              ,(bugetqc(5,mt,k),k=1,nzg)
        read(60,99)timeid,obsq12(1,mt,:),obsq12(2,mt,:)
        read(82,98)timeid,dat43(1,mt,:),dat43(2,mt,:),dat43(3,mt,:)
     +            ,dat43(4,mt,:),dat43(5,mt,:),dat43(6,mt,:)
        write(62,99)(simq12(1,mt,k),k=1,nzg),(simq12(2,mt,k),k=1,nzg)
      enddo
      do iv=1,2
!      tmpwork(:,:)=obsq12(iv,:,:)
       do mt=1,it
          do k=1,nz
            tmpwork(mt,k)=obsq12(iv,mt,k)
          enddo
        enddo
!      tmpwork(:,1)=0.0   !!!!! the bottom level is zero
      call gridint(tmpwork,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
         do  k=1,nzg
           do  i=1,nxg
             obsq12(iv,i,k)=tmpwork(i,k)
c     to(i,k)=to(i,k)-uo(i,k)
            enddo
           enddo
      enddo
!
      do iv=1,6
!      tmpwork(:,:)=obsq12(iv,:,:)
       do mt=1,it
          do k=1,nz
            tmpwork(mt,k)=dat43(iv,mt,k)
          enddo
        enddo
!      tmpwork(:,1)=0.0   !!!!! the bottom level is zero
      call gridint(tmpwork,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
         do  k=1,nzg
           do  i=1,nxg
              dat43(iv,i,k)=tmpwork(i,k)
c     to(i,k)=to(i,k)-uo(i,k)
            enddo
           enddo
      enddo
      do mt=1,it
         write(83,99)(dat43(1,mt,k),k=1,nzg),(dat43(2,mt,k),k=1,nzg),
     +               (dat43(3,mt,k),k=1,nzg),(dat43(4,mt,k),k=1,nzg),
     +               (dat43(5,mt,k),k=1,nzg),(dat43(6,mt,k),k=1,nzg)
         write(61,99)(obsq12(1,mt,k),k=1,nzg),(obsq12(2,mt,k),k=1,nzg)
         read(80,99)timeid,obsftq(1,mt,:),obsftq(2,mt,:)
      enddo
      do iv=1,2
!      tmpwork(:,:)=obsftq(iv,:,:)
       do mt=1,it
          do k=1,nz
            tmpwork(mt,k)=obsftq(iv,mt,k)
          enddo
        enddo
!      tmpwork(:,1)=0.0   !!!!! the bottom level is zero
       call gridint(tmpwork,nx,nz,x,zm,nxg,nzg,xg,zg,work1,work2)
         do  k=1,nzg
           do  i=1,nxg
             obsftq(iv,i,k)=tmpwork(i,k)
c     to(i,k)=to(i,k)-uo(i,k)
            enddo
           enddo
      enddo
      do mt=1,it
         write(81,99)(obsftq(1,mt,k),k=1,nzg),(obsftq(2,mt,k),k=1,nzg)
      enddo
C      (ain,nx,nz,xx,zz,nx1,nz1,xx1,zz1,x0,z0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc this subroutine performs interpolation of the data given in the
cc input array ain(nx,nz) on the grid defined by xx(nx) and zz(nz)
cc into the grid given by xx1(nx1) and zz1(nz1). NIETHER OF THE
cc GRIDS HAS TO BE REGULAR. Data is returned in the ain(nx1,nz1)
cc part of the input array.
cc 
cc    levels in the input array are given in zz(nz), 
cc    levels in the output array are given in zz1(nz1)
cc      x-coordinate in the input array are in xx(nx)
cc      x-coordinate in the output array are in xx1(nx1)
cc        x0(nx,nz) and z0(nx,nz) are working arrays
cc
cc NOTE that nx1 (or nz1) must be smaller than nx (or nz) and xx1 (zz1)
cc  must be a subdomain of xx (zz)
  99  format(8e12.4)
  98  format(8e13.5)    
      STOP
      END
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine gridint
     1 (ain,nx,nz,xx,zz,nx1,nz1,xx1,zz1,x0,z0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc this subroutine performs interpolation of the data given in the
cc input array ain(nx,nz) on the grid defined by xx(nx) and zz(nz)
cc into the grid given by xx1(nx1) and zz1(nz1). NIETHER OF THE
cc GRIDS HAS TO BE REGULAR. Data is returned in the ain(nx1,nz1)
cc part of the input array.
cc 
cc    levels in the input array are given in zz(nz), 
cc    levels in the output array are given in zz1(nz1)
cc      x-coordinate in the input array are in xx(nx)
cc      x-coordinate in the output array are in xx1(nx1)
cc        x0(nx,nz) and z0(nx,nz) are working arrays
cc
cc NOTE that nx1 (or nz1) must be smaller than nx (or nz) and xx1 (zz1)
cc  must be a subdomain of xx (zz)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension ain(nx,nz),zz(nz),xx(nx)
      dimension zz1(nz1),xx1(nx1) 
      dimension x0(nx,nz),z0(nx,nz)

c      check consistency of the input data
      ier=0
cc nx1,nz1 not larger than nx,nz
      if(nx1.gt.nx) ier=1
      if(nz1.gt.nz) ier=2
cc limits (zz1(nz1).le.zz(nz) ?) 
      if(zz1(1).lt.zz(1)) ier=3
      if(zz1(nz1).gt.zz(nz)) ier=4
cc limits (xx1(nx1).le.xx(nx) ?) 
      if(xx1(1).lt.xx(1)) ier=5
      if(xx1(nx1).gt.xx(nx)) ier=6
      if(ier.ne.0) then
      print 999,ier
 999  format(2x,' ** problems with input data. will stop.'/
     1 ' ier = ',i3,'. see code why stoped.')
      stop
      endif
cc
      nxz=nx*nz
      do 99 i=1,nxz
      z0(i,1)=1.
  99  x0(i,1)=1.
cc  map vertical grid positions:
      do 1 k1=1,nz1
      zzh=zz1(k1)
      do 2 k=1,nz
      kk=k
      if(zz(k).ge.zzh) go to 6
  2   continue
  6   kkm=max0(1,kk-1)
      z0(1,k1)=float(kkm)+(zzh-zz(kkm))/(zz(kk)-zz(kkm)+1.e-6)
  1   continue
      do 3 i1=2,nx1
      do 3 k1=1,nz1
  3   z0(i1,k1)=z0(1,k1)
c
cc  map horizontal grid positions:
      do 11 i1=1,nx1
      xxh=xx1(i1)
      do 12 i=1,nx
      ii=i
      if(xx(i).ge.xxh) go to 16
 12   continue
 16   iim=max0(1,ii-1)
      x0(i1,1)=float(iim)+(xxh-xx(iim))/(xx(ii)-xx(iim)+1.e-6)
 11   continue
      do 13 i1=1,nx1
      do 13 k1=2,nz1
 13   x0(i1,k1)=x0(i1,1)
cc
cc  call Piotr's interpolation routine
c      print*,x0(:,1)
c      print*,z0(1,:)
      call inter2(ain,x0,z0,nx,nz)
cc
      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INTER2(XF,XD1,XD2,NX,NZ)
C IOR=ORDER OF ACCURACY/2; ONLY EVEN ORDER TRMBACK SCHEMES ARE CONSIDERED
      PARAMETER(IOR=2)
      PARAMETER(LINER=0)
      DIMENSION XF(*),XD1(*),XD2(*)
CC  N1 - HORIZONTAL INDEX, N2 - VERTICAL INDEX
      PARAMETER (N1=365*4, N2=51,NN=N1*N2)  !!!!!!! N1 is im in, N2 is km-1
      DIMENSION Z(NN,-IOR:IOR)
      DATA  EP/ 1.E-10/
      PARAMETER(NONOS=1)
      REAL      MX,MN
      PARAMETER(IBC=0)
      COMMON // IG0(NN),JG0(NN),X(-IOR+1:N1+IOR,-IOR+1:N2+IOR)
C  next is for shavano: 
C      DONOR(Y1,Y2,A)=CVMGM(Y2,Y1,A)*A
C  next is for workstation:
      DONOR(Y1,Y2,A)=AMAX1(0.,A)*Y1 + AMIN1(0.,A)*Y2
      TR2(Y1,Y2,A)=A*.5*(Y1+Y2)-A**2*.5*(Y2-Y1)
      TR4(YM1,Y0,YP1,YP2,A)=A/12.*(7.*(YP1+Y0)-(YP2+YM1))
     1 -A**2/24.*(15.*(YP1-Y0)-(YP2-YM1))-A**3/12.*((YP1+Y0)
     2 -(YP2+YM1))+A**4/24.*(3.*(YP1-Y0)-(YP2-YM1))
      TR6(YM2,YM1,Y0,YP1,YP2,YP3,A)=-A/60.*(-YM2+8.*YM1-37.*Y0
     1                                     -37.*YP1+8.*YP2-YP3)
     2-A**2/360.*(-2.*YM2+25.*YM1-245.*Y0+245.*YP1-25.*YP2+2.*YP3)
     3-A**3/48.*(YM2-7.*YM1+6.*Y0+6.*YP1-7.*YP2+YP3)
     4-A**4/144.*(YM2-11.*YM1+28.*Y0-28.*YP1+11.*YP2-YP3)
     5-A**5/240.*(-YM2+3.*YM1-2.*Y0-2.*YP1+3.*YP2-YP3)
     6-A**6/720.*(-YM2+5.*YM1-10.*Y0+10.*YP1-5.*YP2+YP3)
      PP(XI)=AMAX1(0.,XI)
      PN(XI)=AMIN1(0.,XI)
C
CC CHECK COSISTENCY OF THE DATA:
      IF(NX.NE.N1.OR.NZ.NE.N2) THEN
      PRINT 777
 777  FORMAT(2X,'!!! CALLS TO INTER2 WITH NON-MATCHING DIMENSIONS.'
     1  ,' STOP.')
      STOP
      ENDIF
CC
      DO 1 K=1,NN
      IG0(K)=NINT(XD1(K))
    1 JG0(K)=NINT(XD2(K))

C  GRID EXTENSION FOR BC REMOVAL 
      DO 508 I=1,N1      
      DO 509 J=1,N2
      II=(J-1)*N1+I
  509 X(I,J)=XF(II) 
      DO 5091 IS=1-IOR,0
C     II=(1-1)*N1+I
 5091 X(I,IS)=XF(I)
      DO 5092 IS=1,IOR
      II=(N2-1)*N1+I
 5092 X(I,N2+IS)=XF(II)
  508 CONTINUE
      DO 507 J=-IOR+1,N2+IOR
      DO 5071 IS=-IOR+1,1
 5071 X(IS,J)=X(1,J)*(1-IBC)+IBC*X(N1+IS-1,J)
      DO 5072 IS=0,IOR
 5072 X(N1+IS,J)=X(N1,J)*(1-IBC)+IBC*X(1+IS,J)
  507 CONTINUE
C  END OF GRID EXTENSION
C                     
C
C  HERE STARTS REZIDUAL ADVECTION
C
                     DO 50 J=-IOR,IOR
C
      IF(LINER.EQ.1) THEN
      DO 211 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
  211 Z(II,J)=Y0-(FL1-FL0) 
      GO TO 50
      ENDIF
C
      IF(IOR.EQ.1) THEN
        IF(NONOS.EQ.1) THEN
      DO 311 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      F0=TR2(YM1, Y0,U)
      F1=TR2(Y0 ,YP1,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  311 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 321 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      F0=TR2(YM1, Y0,U)
      F1=TR2(Y0 ,YP1,U)
  321 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
      IF(IOR.EQ.2) THEN
        IF(NONOS.EQ.1) THEN
      DO 312 II=1,NN
      U=IG0(II)-XD1(II)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      F0=TR4(YM2,YM1,Y0 ,YP1,U)
      F1=TR4(YM1,Y0 ,YP1,YP2,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  312 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 322 II=1,NN
      U=IG0(II)-XD1(II)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      F0=TR4(YM2,YM1,Y0 ,YP1,U)
      F1=TR4(YM1,Y0 ,YP1,YP2,U)
  322 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
      IF(IOR.EQ.3) THEN
        IF(NONOS.EQ.1) THEN
      DO 313 II=1,NN
      U=IG0(II)-XD1(II)
      YM3=X(IG0(II)-3,JG0(II)+J)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      YP3=X(IG0(II)+2,JG0(II)+J)
      F0=TR6(YM3,YM2,YM1,Y0 ,YP1,YP2,U)
      F1=TR6(YM2,YM1,Y0 ,YP1,YP2,YP3,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  313 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 323 II=1,NN
      U=IG0(II)-XD1(II)
      YM3=X(IG0(II)-3,JG0(II)+J)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      YP3=X(IG0(II)+2,JG0(II)+J)
      F0=TR6(YM3,YM2,YM1,Y0 ,YP1,YP2,U)
      F1=TR6(YM2,YM1,Y0 ,YP1,YP2,YP3,U)
  323 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
C
   50 CONTINUE
C  
      IF(LINER.EQ.1) THEN
      DO 212 II=1,NN
      U=JG0(II)-XD2(II)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
  212 XF(II)=Z(II,0)-(FL1-FL0) 
      RETURN
      ENDIF
C
      IF(IOR.EQ.1) THEN
        IF(NONOS.EQ.1) THEN
      DO 411 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR2(Z(II,-1),Z(II,0),U)
      F1=TR2(Z(II, 0),Z(II,1),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  411 CONTINUE
        ELSE
      DO 421 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR2(Z(II,-1),Z(II,0),U)
      F1=TR2(Z(II, 0),Z(II,1),U)
  421 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF

      IF(IOR.EQ.2) THEN
        IF(NONOS.EQ.1) THEN
      DO 412 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR4(Z(II,-2),Z(II,-1),Z(II,0),Z(II,1),U)
      F1=TR4(Z(II,-1),Z(II, 0),Z(II,1),Z(II,2),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  412 CONTINUE
        ELSE
      DO 422 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR4(Z(II,-2),Z(II,-1),Z(II,0),Z(II,1),U)
      F1=TR4(Z(II,-1),Z(II, 0),Z(II,1),Z(II,2),U)
  422 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF

      IF(IOR.EQ.3) THEN
        IF(NONOS.EQ.1) THEN
      DO 413 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR6(Z(II,-3),Z(II,-2),Z(II,-1),Z(II,0),
     1                     Z(II, 1),Z(II, 2),U)
      F1=TR6(Z(II,-2),Z(II,-1),Z(II, 0),Z(II,1),
     1                     Z(II, 2),Z(II, 3),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  413 CONTINUE
        ELSE
      DO 423 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR6(Z(II,-3),Z(II,-2),Z(II,-1),Z(II,0),
     1                     Z(II, 1),Z(II, 2),U)
      F1=TR6(Z(II,-2),Z(II,-1),Z(II, 0),Z(II,1),
     1                     Z(II, 2),Z(II, 3),U)
  423 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF
      RETURN
      END   
