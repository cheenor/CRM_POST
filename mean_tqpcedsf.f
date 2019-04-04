      PARAMETER (im=202, km=52, inn=72, itt=2880, ii=121) !!!! 4*30+1
      dimension t(im),q(im),p(im)
      dimension c(im,km),e(im,km),d(im,km),s(im,km),f(im,km)
      dimension smt(itt),smq(itt),smp(itt),smd(itt,km)
      dimension smc(itt,km),sme(itt,km),sms(itt,km),smf(itt,km)
      dimension tmt(ii),tmq(ii),tmp(ii),tmd(ii,km)
      dimension tmc(ii,km),tme(ii,km),tms(ii,km),tmf(ii,km)
      character casenm*20, runfold*20
      casenm='casename'
      runfold='runfold'
      open(21,file=
     *'../'//trim(runfold)//'/postdata/micro_202_'//trim(casenm)//'.txt'
     *)
      open(30,file=
     *'../'//trim(runfold)//'/postdata/'//trim(casenm)//
     *'_micro_202_6hour.txt')
      it=0
      do 999 ia=1,1 !7
c     if(ia.le.5) then
      sec=15.
c     else
c     sec=12.
c     endif
      IH=21
      do 100 in=1,itt   ! inn
      it=it+1
      do i=1,im
!      READ(IH,99,end=998) ktype,(t(i),i=1,im)
!      READ(IH,99,end=998) ktype,(q(i),i=1,im)
!      READ(IH,99,end=998) (p(i),i=1,im)
      READ(IH,99,end=998) (c(i,k),k=1,km)
      READ(IH,99,end=998) (e(i,k),k=1,km)
      READ(IH,99,end=998) (d(i,k),k=1,km)
      READ(IH,99,end=998) (s(i,k),k=1,km)
      READ(IH,99,end=998) (f(i,k),k=1,km)
      enddo
      smt(it)=0.
      smq(it)=0.
      smp(it)=0.
      do 31 i=2,im-1
      smt(it)=smt(it)+t(i)/float(im-2) 
      smq(it)=smq(it)+q(i)/float(im-2) 
      smp(it)=smp(it)+p(i)/float(im-2) 
  31  continue
      do 40 k=1,km
      smc(it,k)=0.
      sme(it,k)=0.
      smd(it,k)=0.
      sms(it,k)=0.
      smf(it,k)=0.
      do 30 i=2,im-1
      smc(it,k)=smc(it,k)+c(i,k)/sec/float(im-2) 
      sme(it,k)=sme(it,k)+e(i,k)/sec/float(im-2) 
      smd(it,k)=smd(it,k)+d(i,k)/sec/float(im-2) 
      sms(it,k)=sms(it,k)+s(i,k)/sec/float(im-2) 
      smf(it,k)=smf(it,k)+f(i,k)/sec/float(im-2) 
  30  continue
  40  continue
 100  continue   
 999  continue
!      do 51 mt=2,ii
!      n1=1+(mt-2)*24 !!! 15 mins every raw timestep
!      n2=(mt-1)*24   !!!!
!      tmt(mt)=0.
!      tmq(mt)=0.
!      tmp(mt)=0.
!      do 12 n=n1,n2 
!      tmt(mt)=tmt(mt)+smt(n)/float(n2-n1+1)
!      tmq(mt)=tmq(mt)+smq(n)/float(n2-n1+1)
!      tmp(mt)=tmp(mt)+smp(n)*1000.*3600./float(n2-n1+1)
!  12  continue
!  51  continue
      tmt(1)=0.
      tmq(1)=0.
      tmp(1)=0.
      do 200 k=1,km
      do 50 mt=1,ii    !!! mean 
      n1=(mt-1)*24-11             !1+(mt-2)*9
      n2=(mt-1)*24+12             !(mt-1)*9
      if(mt.eq.1)then
        n1=1
        n2=12
      end if
      if(mt.eq.ii)then
        n1=itt-11
        n2=itt
      endif
      tmc(mt,k)=0.
      tme(mt,k)=0.
      tmd(mt,k)=0.
      tms(mt,k)=0.
      tmf(mt,k)=0.
      do 11 n=n1,n2 
      tmc(mt,k)=tmc(mt,k)+smc(n,k)*3600.*24./float(n2-n1+1)
      tme(mt,k)=tme(mt,k)+sme(n,k)*3600.*24./float(n2-n1+1)
      tmd(mt,k)=tmd(mt,k)+smd(n,k)*3600.*24./float(n2-n1+1)
      tms(mt,k)=tms(mt,k)+sms(n,k)*3600.*24./float(n2-n1+1)
      tmf(mt,k)=tmf(mt,k)+smf(n,k)*3600.*24./float(n2-n1+1)
  11  continue
  50  continue
 200  continue
      do 400 k=1,km
      tmc(1,k)=0.
      tme(1,k)=0.
      tmd(1,k)=0.
      tms(1,k)=0.
      tmf(1,k)=0.
  400 continue
!      write(30,99) tmt,tmq
!      write(30,99) tmp
!      write(30,99) tmc,tme,tmd,tms,tmf
!      write(30,99) tmt,tmq
!      write(30,99) tmp
       do it=1,ii
       write(30,99) tmc(it,:),tme(it,:),tmd(it,:),tms(it,:),tmf(it,:)
       end do
99    format(8e12.4)
c  units:                
c     sensible and latent heat flux (tmt and tmq  W/m**2)
c     precipitation (tmp  mm/h)
c     con, eva, dep, sub and fus (tmc, tme, tmd, tms and tmf  K/day)
  998 continue
      stop
      end
