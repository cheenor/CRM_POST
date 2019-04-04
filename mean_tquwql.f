      PARAMETER (im=202, km=52, inn=480, itt=35040,ii=4*365,ifi=73)
      dimension qc(im,km),qr(im,km),qa(im,km),qb(im,km),rho(km)
      dimension te(km),ps(km),potf(km),pote(km)
      dimension u(im,km),w(im,km),t(im,km),q(im,km),rh(im,km)
      dimension smqc(itt,km),smqr(itt,km),smqa(itt,km),smqb(itt,km)
      dimension smtc(itt,km),tmtc(ii,km)
      dimension smu(itt,km),smw(itt,km),smt(itt,km),smq(itt,km)
      dimension smrh(itt,km),smrat(itt,km),smfqv(itt,km),smftc(itt,km)
      dimension tmqc(ii,km),tmqr(ii,km),tmqa(ii,km),tmqb(ii,km)
      dimension tmu(ii,km),tmw(ii,km),tmt(ii,km),tmq(ii,km)
      dimension tmrh(ii,km),tmrat(ii,km),tmfqv(ii,km),tmftc(ii,km)
      character*100 chenm, path ! case name ,must be recorded 
      character casenm*20,fold*20,idstr*3
      casenm='ETP2D0'
      fold='run2'
      if(casenm(1:3)=='MLY')then
        path='/public/home/chenjh/Models/CRM/ERA/Year/'//casenm(1:4)//
     +   '/'//trim(fold)
      else
        path='/public/home/chenjh/Models/CRM/ERA/Year/'//casenm(1:3)//
     +   '/'//trim(fold)
      endif
!      path='/home/jhchen/jhchen/ERA_Interim/ETP/'//trim(fold)
      chenm=trim(path)//'/'//trim(casenm)    !'whereisthe_data_tgn2d1_(1-6)'      
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
!      open(40,file=trim(chenm)//'_All'
!     * ,form='unformatted')
      open(50,file=trim(chenm)//'_All.txt')
      open(51,file=trim(chenm)//'_rho.txt')
      open(52,file=trim(chenm)//'_Raw_qcqaqbqr.txt')
      open(53,file=trim(chenm)//'_Raw_tqrhuvw.txt')
      it=0
      do 999 if=1,ifi
      IH=80+if
      write(idstr,'(I3.3)')if
      open(IH,file=trim(chenm)//'_'//idstr
     * ,form='UNFORMATTED',status='OLD',convert='big_endian')
      do 100 in=1,inn
      it=it+1
      read(IH) qc,qr 
      read(IH) qa,qb 
      read(IH) t,q
      read(IH) rh,te,ps
      read(IH) u,w
      read(IH) 
      read(IH) rho,potf,pote
      read(IH) 
C
      do i=2,im-1
      write(52,99)qc(i,:)*1000,qa(i,:)*1000
     +            ,qb(i,:)*1000,qr(i,:)*1000
      write(53,99)t(i,:),q(i,:)*1000,rh(i,:)*100
     +            ,u(i,:)*10,w(i,:)*10 
      end do
C
      do 40 k=1,km
      write(51,*)k,rho(k)
      smqc(it,k)=0.
      smqr(it,k)=0.
      smqa(it,k)=0.
      smqb(it,k)=0.
      smrh(it,k)=0.
      smt(it,k)=0.
      smq(it,k)=0.
      smu(it,k)=0.
      smw(it,k)=0.
      smtc(it,k)=0.
      do 30 i=2,im-1
      smqc(it,k)=smqc(it,k)+qc(i,k)*1000./float(im-2) 
      smqr(it,k)=smqr(it,k)+qr(i,k)*1000./float(im-2) 
      smqa(it,k)=smqa(it,k)+qa(i,k)*1000./float(im-2) 
      smqb(it,k)=smqb(it,k)+qb(i,k)*1000./float(im-2) 
      smrh(it,k)=smrh(it,k)+rh(i,k)*100./float(im-2)
      tmpetc=t(i,k)*te(k)/potf(k)/(1.+pote(k))-273.16  !!!!  C degree
      smtc(it,k)=smtc(it,k)+tmpetc/float(im-2) 
      smt(it,k)=smt(it,k)+t(i,k)/float(im-2) 
      smq(it,k)=smq(it,k)+q(i,k)*1000./float(im-2) 
      smu(it,k)=smu(it,k)+0.5*(u(i-1,k)+u(i,k))*10./float(im-2) 
      smw(it,k)=smw(it,k)+w(i,k)*10./float(im-2) 
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
      tmqc(mt,k)=0.
      tmqr(mt,k)=0.
      tmqa(mt,k)=0.
      tmqb(mt,k)=0.
      tmrh(mt,k)=0.
      tmtc(mt,k)=0.
      tmt(mt,k)=0.
      tmq(mt,k)=0.
      tmu(mt,k)=0.
      tmw(mt,k)=0.
      do 11 n=n1,n2 
      tmqc(mt,k)=tmqc(mt,k)+smqc(n,k)/float(n2-n1+1)
      tmqr(mt,k)=tmqr(mt,k)+smqr(n,k)/float(n2-n1+1)
      tmqa(mt,k)=tmqa(mt,k)+smqa(n,k)/float(n2-n1+1)
      tmqb(mt,k)=tmqb(mt,k)+smqb(n,k)/float(n2-n1+1)
      tmrh(mt,k)=tmrh(mt,k)+smrh(n,k)/float(n2-n1+1)
      tmt(mt,k)=tmt(mt,k)+smt(n,k)/float(n2-n1+1)
      tmq(mt,k)=tmq(mt,k)+smq(n,k)/float(n2-n1+1)
      tmu(mt,k)=tmu(mt,k)+smu(n,k)/float(n2-n1+1)
      tmw(mt,k)=tmw(mt,k)+smw(n,k)/float(n2-n1+1)
      tmtc(mt,k)=tmtc(mt,k)+smtc(n,k)/float(n2-n1+1)
  11  continue
  50  continue
 200  continue
      write(40) tmt,tmq,tmu,tmw
      write(40) tmqc,tmqr,tmqa,tmqb
      write(40) tmrh
      write(40) tmfqv,tmftc,tmrat
      print*,'tmt(1,1)= , tmq(1,1)= ',tmt(1,1),tmq(1,1)
      print*,'tmt(80,1)= , tmq(80,1)= ',tmt(80,1),tmq(80,1)
c
      do mt=1,ii
      write(50,99) tmt(mt,:),tmq(mt,:),tmu(mt,:),tmw(mt,:),tmqc(mt,:),
     +tmqr(mt,:),tmqa(mt,:),tmqb(mt,:), tmrh(mt,:),
     +tmfqv(mt,:),tmftc(mt,:),tmrat(mt,:),tmtc(mt,:)
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
