      PARAMETER (nday=30, im=201, itt=35040, ii=8*365, ifi=73)
      dimension upsws(im),dnsws(im),uplws(im),dnlws(im)
      dimension upswt(im),dnswt(im),uplwt(im),dnlwt(im)
       dimension upswtc(im),dnswtc(im),uplwtc(im),dnlwtc(im)
       dimension upswsc(im),dnswsc(im),uplwsc(im),dnlwsc(im)
      dimension upss(itt),dnss(itt),upls(itt),dnls(itt)
      dimension upst(itt),dnst(itt),uplt(itt),dnlt(itt)
      dimension aupss(ii),adnss(ii),aupls(ii),adnls(ii)
      dimension aupst(ii),adnst(ii),auplt(ii),adnlt(ii)
      real rsw(itt),sw(itt)
      real arsw(ii),asw(ii)
      real brsw(nday)
      character casenm*20,path*100,fold*20,idstr*3
      casenm='ETP2D0'
      fold='run2'
      if(casenm(1:3)=='MLY')then
        path='/public/home/chenjh/Models/CRM/ERA/Year/'//casenm(1:4)//
     +   '/'//trim(fold)//'/postdata/'
      else
        path='/public/home/chenjh/Models/CRM/ERA/Year/'//casenm(1:3)//
     +   '/'//trim(fold)//'/postdata/'
      endif
!      path ='/home/jhchen/jhchen/ERA_Interim/ETP/'//trim(fold)//
!     +'/postdata/'
      open(81,file=
     *trim(path)//'radflux_'//trim(casenm)
     * ,form='UNFORMATTED',status='OLD') 
c     open(60,file=
c    1'/mnt/raid51/wuxq/newarm/2darm_data/arm2d1_rfl_r30'
c    1,  form='unformatted',status='old')
c     open(65,file=
c    1'/mnt/raid51/wuxq/newarm/2darm_data/arm2d1_rfl_r30_cl'
c    1,  form='unformatted',status='old')
      open(20,file=trim(path)//'rad_3hr_'//trim(casenm),
     1 form='formatted')
      open(40,file=
     1 trim(path)//'racsw97_'//trim(casenm)//'_3hr.dat',
     1 form='formatted')
      open(45,file=
     1 trim(path)//'racsw97_'//trim(casenm)//'_day.dat',
     1 form='formatted')
      it=0
      ars=0.
      ars2=0.
      tt=0.
      do 999 jf=1,ifi
      IH=81
      inn=480
c     if(jf.eq.ifi) inn=576
      do 100 in=1,inn
      it=it+1
      read(81) (upsws(j),j=1,im),
     .         (dnsws(j),j=1,im),
     .         (uplws(j),j=1,im),
     .         (dnlws(j),j=1,im)
      read(81) (upswt(j),j=1,im),
     .         (dnswt(j),j=1,im),
     .         (uplwt(j),j=1,im),
     .         (dnlwt(j),j=1,im)
c      read(60)  (upsws(j),j=1,im),
c    .           (dnsws(j),j=1,im),
c    .           (uplws(j),j=1,im),
c    .           (dnlws(j),j=1,im)
c      read(60)  (upswt(j),j=1,im),
c    .           (dnswt(j),j=1,im),
c    .           (uplwt(j),j=1,im),
c    .           (dnlwt(j),j=1,im)
c      read(60) 
c      read(60) 
c      read(65)  (upswsc(j),j=1,im),
c    .           (dnswsc(j),j=1,im),
c    .           (uplwsc(j),j=1,im),
c    .           (dnlwsc(j),j=1,im)
c      read(65)  (upswtc(j),j=1,im),
c    .           (dnswtc(j),j=1,im),
c    .           (uplwtc(j),j=1,im),
c    .           (dnlwtc(j),j=1,im)
c      read(65) 
c      read(65)
c       print*,dnsws 
c       print*,dnsws
      upss(it)=0.
      dnss(it)=0.
      upls(it)=0.
      dnls(it)=0.
      upst(it)=0.
      dnst(it)=0.
      uplt(it)=0.
      dnlt(it)=0.
      rsw(it)=0.
      sw(it)=0.
      do 30 i=2,im-1
      upss(it)=upss(it)+upsws(i)*1.e-3/float(im-2)
      dnss(it)=dnss(it)+dnsws(i)*1.e-3/float(im-2)
      upls(it)=upls(it)+uplws(i)*1.e-3/float(im-2)
      dnls(it)=dnls(it)+dnlws(i)*1.e-3/float(im-2)
      upst(it)=upst(it)+upswt(i)*1.e-3/float(im-2)
      dnst(it)=dnst(it)+dnswt(i)*1.e-3/float(im-2)
      uplt(it)=uplt(it)+uplwt(i)*1.e-3/float(im-2)
      dnlt(it)=dnlt(it)+dnlwt(i)*1.e-3/float(im-2)
c       print*,dnswsc(i)
      if(dnswsc(i).ne.0.) then
      rsw(it)=rsw(it)+dnsws(i)/dnswsc(i)
      sw(it)=sw(it)+dnsws(i)*1.e-3
      print*,sw(i),rsw(i)
      endif
  30  continue
      rsw(it)=rsw(it)/float(im-2)
      sw(it)=sw(it)/float(im-2)
c     print*, rsw(it),sw(it)

      if(rsw(it).gt.0. .and. rsw(it).le.1.) then
      ars=ars+rsw(it)
      ars2=ars2+rsw(it)*rsw(it)
      tt=tt+1.
      endif
 100  continue   
 999  continue  
      print*, tt,ars,ars2 
      ars=ars/tt
      ars2=ars2/tt
      std=sqrt(abs(ars2-ars**2))
      print*,'30-day mean and std dnsws from CRM', ars,std

      do 50 mt=1,ii
      n1=(mt-1)*12-7
      n2=(mt-1)*12+4
      if(mt.eq.1) then
      n1=1
      n2=4
      endif
      if(mt.eq.ii) then
      n1=itt-7
      n2=itt
      endif
      aupss(mt)=0.
      adnss(mt)=0.
      aupls(mt)=0.
      adnss(mt)=0.
      aupst(mt)=0.
      adnst(mt)=0.
      auplt(mt)=0.
      adnst(mt)=0.
      arsw(mt)=0.
      asw(mt)=0.
      tt=0.
      do 11 n=n1,n2 
      aupss(mt)=aupss(mt)+upss(n)/float(n2-n1+1)
      adnss(mt)=adnss(mt)+dnss(n)/float(n2-n1+1)
      aupls(mt)=aupls(mt)+upls(n)/float(n2-n1+1)
      adnls(mt)=adnls(mt)+dnls(n)/float(n2-n1+1)
      aupst(mt)=aupst(mt)+upst(n)/float(n2-n1+1)
      adnst(mt)=adnst(mt)+dnst(n)/float(n2-n1+1)
      auplt(mt)=auplt(mt)+uplt(n)/float(n2-n1+1)
      adnlt(mt)=adnlt(mt)+dnlt(n)/float(n2-n1+1)
c     arsw(mt)=arsw(mt)+rsw(n)/float(n2-n1+1)
c     asw(mt)=asw(mt)+sw(n)/float(n2-n1+1)
      if(rsw(n).gt.0. .and. rsw(n).le.1.) then
      arsw(mt)=arsw(mt)+rsw(n)
      asw(mt)=asw(mt)+sw(n)
      tt=tt+1.
      endif
  11  continue
      if(tt.ne.0.) then
      arsw(mt)=arsw(mt)/tt
      asw(mt)=asw(mt)/tt
      endif
      write(20,776) aupst(mt),adnst(mt),auplt(mt),adnlt(mt) 
      write(20,776) aupss(mt),adnss(mt),aupls(mt),adnls(mt) 
 776  format(1x,4e18.8)
      write(40,876) arsw(mt),asw(mt),tt
 876  format(1x,3e18.8)
  50  continue
      do it=1,nday      ! to save 26-day daily toa fluxes, nday should be 26
       brsw(it)=0.
       bt=0.
       do j=1,96
       if(rsw((it-1)*96+j).gt.0. .and. rsw((it-1)*96+j).le.1.) then
       brsw(it)=brsw(it)+rsw((it-1)*96+j)
       bt=bt+1.
       endif
       enddo
       brsw(it)=brsw(it)/bt
       write(45,676) brsw(it)
 676  format(1x,e18.8)
       enddo
c  units:                
c     rad flux (W/m**2)
      stop
      end
