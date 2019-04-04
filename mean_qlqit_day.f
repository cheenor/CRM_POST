      PARAMETER (im=202, km=52, inn=480, itt=480*73,ii=20*73,ifi=73) ! every file stores 5 days' data
      PARAMETER (ndays=ifi*5)
      parameter (imm=200,kmm=34)
      dimension qc(im,km),qr(im,km),qa(im,km),qb(im,km),rho(km)
      dimension u(im,km),w(im,km),th(im,km),q(im,km),rh(im,km)
!      dimension te(im,km)
      dimension te(km),ps(km),potf(km),pote(km)
c
      dimension smql(ndays,kmm),smqi(ndays,kmm),smtc(ndays,kmm)
      dimension den(ndays,kmm),smqc(ndays,kmm)
      dimension smqa(ndays,kmm),smqb(ndays,kmm),smqr(ndays,kmm)
      !dimension tmql(ii,kmm),tmqi(ii,kmm),tmtc(ii,kmm)
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
c      open(40,file=trim(chenm)//'_All'
c     * ,form='unformatted')
!       open(45,file=trim(chenm)//'_qlqit_41-73',form='unformatted')
       open(50,file=trim(chenm)//'_qlqit_day.txt')
c      open(51,file=trim(chenm)//'_rho.txt')
cc      open(52,file=trim(chenm)//'_Raw_qcqaqbqr.txt')
c      open(53,file=trim(chenm)//'_Raw_tqrhuvw.txt')
      it=0
      iday=1
      smql=0.0
      smqi=0.0
      smtc=0.0
      den=0.0
      smqc=0.0
      smqa=0.0
      smqb=0.0
      smqr=0.0      
      do 999 if=1,ifi
      IH=80+if
      write(idstr,'(I3.3)')if
      chenm=trim(path)//'/'//trim(casenm)//'_'//idstr    !'whereisthe_data_tgn2d1_(1-6)'
      open(IH,file=trim(chenm),form='unformatted',
     +  status='old',convert='big_endian')

      do 100 in=1,inn
      it=it+1
      read(IH) qc,qr 
      read(IH) qa,qb 
      read(IH) th,q
      read(IH) rh,te,ps
      read(IH) u,w
      read(IH) 
      read(IH) rho,potf,pote
      read(IH) 
C
C      
      do 40 k=1,kmm
      den(iday,k)=den(iday,k)+rho(k)/96.0
      do 30 i=2,im-1
      ij=ij+1
      tmql=qc(i,k)*1000.+ qr(i,k)*1000.
      tmqi=qa(i,k)*1000.+qb(i,k)*1000.
      tmptc=th(i,k)*te(k)/potf(k)/(1.+pote(k))-273.16 ! convert theta to temp in C degree
      smtc(iday,k)=smtc(iday,k)+tmptc/(96.0*imm)
      smql(iday,k)=smql(iday,k)+tmpql/(96.0*imm)
      smqi(iday,k)=smqi(iday,k)+tmpqi/(96.0*imm)
      !!!
      smqc(iday,k)=smqc(iday,k)+qc(i,k)*1000./(96.0*imm)
      smqr(iday,k)=smqr(iday,k)+qr(i,k)*1000./(96.0*imm)
      smqa(iday,k)=smqa(iday,k)+qa(i,k)*1000./(96.0*imm)
      smqb(iday,k)=smqb(iday,k)+qb(i,k)*1000./(96.0*imm)
  30  continue
  40  continue
      if (it==96)then
        iday=iday+1
        it=0
      endif
 100  continue   
 999  continue
      print*, iday
      do iday=1,ndays:
        write(50,'I3')iday
        do k=1,kmm
          write(50,99)den(iday,k),smqc(iday,k),smqr(iday,k),smqa(iday,k),smqb(iday,k)
    +          smql(iday,k),smqi(iday,k),smtc(iday,k) 
        enddo
      enddo  
 
!      print*,'tmql(1,1,1)= , tmqi(1,1,1)= ',tmt(1,1,1),tmq(1,1,1)
!      print*,'tmtc(80,5,3)=  ',tmtc(80,5,3)
c
!      do mt=1,ii
!      write(50,99) tmql(mt,:),tmqi(mt,:),tmtc(mt,:),tmden(mt,:) !,tmqc(mt,:),
!     +tmqr(mt,:),tmqa(mt,:),tmqb(mt,:), tmrh(mt,:),
!     +tmfqv(mt,:),tmftc(mt,:),tmrat(mt,:) 
!      enddo
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
