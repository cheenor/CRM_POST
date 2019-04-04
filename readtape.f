      program read_tape 
      PARAMETER (MY=202,MZ=52,mym=my-1, im=202,km=52,nt=35040)
      CHARACTER*40 IFNAME,WU_IFNAME,WU2_IFNAME
      CHARACTER*6 JTAPE
      DATA JTAPE/'ETP2D0'/
      COMMON/SLHF/ FSH(MY),FLH(MY),ust(my),vst(my)
      COMMON/PRCRM/PRECI(MY),pbl(my)
      COMMON/Q1Q2/con(MY,MZ),eva(MY,MZ),dep(MY,MZ)
     1 ,sub(MY,MZ),fus(MY,MZ)
      COMMON/WRKTQ/WRK_th(MY,MZ),WRK_qv(MY,MZ)
      real smh,smt,smq,smp
      real smc(km),sme(km),sms(km),smf(km),smd(km)
cc radiative fluxes: they are at positions of vert. velocity,
cc surface fluxes are at k=1
c     common/radflux/ upsw(mym,mz),dnsw(mym,mz),
c    *                uplw(mym,mz),dnlw(mym,mz)
c     dimension upswt(my),dnswt(my),upsws(my),dnsws(my),
c    *          uplwt(my),uplws(my),dnlws(my)
      common/radflux/upsw(mym,mz),dnsw(mym,mz)
     *              ,uplw(mym,mz),dnlw(mym,mz)
c    *              ,upswt(mym),dnswt(mym),upsws(mym),dnsws(mym)
c    *              ,uplwt(mym),dnlwt(mym),uplws(mym),dnlws(mym)
      character casenm*6,path*100,fold*20
      casenm='ETP2D0'
      fold='run2'
      path='/public/home/chenjh/Models/CRM/ERA/Year/ETP/'//
     +trim(fold)//'/postdata/'
      open(60,file=trim(path)//
     *'preci_'//trim(casenm)
     *,FORM='UNFORMATTED')
      open(61,file=trim(path)//
     *'radflux_'//trim(casenm)
     *,FORM='UNFORMATTED')
      open(70,file=trim(path)//
     *'micro_'//trim(casenm),FORM='UNFORMATTED')
      open(71,file=trim(path)//
     *'surface_'//trim(casenm),FORM='UNFORMATTED')
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!      open(63,file=trim(path)//
!     *'preci_'//trim(casenm)//'.txt')
C     *,FORM='UNFORMATTED')
!      open(64,file=trim(path)//
!     *'radflux_'//trim(casenm)//'.txt')
C     *,FORM='UNFORMATTED')
!      open(73,file=trim(path)//
!     *'micro_'//trim(casenm)//'.DAT',FORM='UNFORMATTED')
!      open(74,file=trim(path)//
!     *'surface_'//trim(casenm)//'.DAT',FORM='UNFORMATTED')

      NVLT=1
      DO 8888 NIOSTR=1,nt
c     DO 8888 NIOSTR=1,2592
      it=NIOSTR
      write(IFNAME,27)JTAPE,NVLT,NIOSTR
      IU=10
      JU=47
      KU=49
 27   FORMAT(A6,'m',I2.2,'nio',I5.5)
      write(WU_IFNAME,28)IFNAME
 28   FORMAT('WUTAPE_',A18)
      write(WU2_IFNAME,29)IFNAME
 29   FORMAT('RADTAPE_',A18)
      write(6,*)'IFNAME  ',IFNAME
      OPEN(UNIT=IU,FILE=IFNAME,STATUS='OLD',FORM='UNFORMATTED'
     +,convert='big_endian')
      OPEN(JU,FILE=WU_IFNAME,STATUS='OLD',FORM='UNFORMATTED'
     +,convert='big_endian')
      OPEN(KU,FILE=WU2_IFNAME,STATUS='OLD',FORM='UNFORMATTED'
     +,convert='big_endian')
C
      READ(JU,IOSTAT=IOVAL) (pbl(J),J=1,MY)
C          print*,pbl
      IF (IOVAL) 500,807,510
  807 CONTINUE
      READ(JU,IOSTAT=IOVAL) (FSH(J),J=1,MY)
C       print*,'ioval',ioval
      IF (IOVAL) 500,808,510
  808 CONTINUE
      READ(JU,IOSTAT=IOVAL) (FLH(J),J=1,MY)
      IF (IOVAL) 500,809,510
  809 CONTINUE
      READ(JU,IOSTAT=IOVAL) (PRECI(J),J=1,MY)
      IF (IOVAL) 500,810,510
  810 CONTINUE
      READ(JU,IOSTAT=IOVAL) ((con(J,K),J=1,MY),K=1,MZ)
      IF (IOVAL) 500,815,510
  815 CONTINUE
      READ(JU,IOSTAT=IOVAL) ((eva(J,K),J=1,MY),K=1,MZ)
      IF (IOVAL) 500,820,510
  820 CONTINUE
      READ(JU,IOSTAT=IOVAL) ((dep(J,K),J=1,MY),K=1,MZ)
      IF (IOVAL) 500,830,510
  830 CONTINUE
      READ(JU,IOSTAT=IOVAL) ((sub(J,K),J=1,MY),K=1,MZ)
      IF (IOVAL) 500,840,510
  840 CONTINUE
      READ(JU,IOSTAT=IOVAL) ((fus(J,K),J=1,MY),K=1,MZ)
      IF (IOVAL) 500,850,510
  850 CONTINUE
c     READ(JU,IOSTAT=IOVAL) ((WRK_th(J,K),J=1,MY),K=1,MZ)
c     IF (IOVAL) 500,860,510
c 860 CONTINUE
c     READ(JU,IOSTAT=IOVAL) ((WRK_qv(J,K),J=1,MY),K=1,MZ)
c     IF (IOVAL) 500,870,510
c 870 CONTINUE

      READ(KU,IOSTAT=IOVAL) ((upsw(J,K),J=1,MYM),K=1,MZ)
      IF (IOVAL) 500,915,510
  915 CONTINUE
      READ(KU,IOSTAT=IOVAL) ((dnsw(J,K),J=1,MYM),K=1,MZ)
      IF (IOVAL) 500,920,510
  920 CONTINUE
      READ(KU,IOSTAT=IOVAL) ((uplw(J,K),J=1,MYM),K=1,MZ)
      IF (IOVAL) 500,930,510
  930 CONTINUE
      READ(KU,IOSTAT=IOVAL) ((dnlw(J,K),J=1,MYM),K=1,MZ)
      IF (IOVAL) 500,940,510
  940 CONTINUE

      smh=0.
      smt=0.
      smq=0.
      smp=0.
      do 31 i=2,im-1
      smh=smh+pbl(i)/float(im-2) 
      smt=smt+fsh(i)/float(im-2) 
      smq=smq+flh(i)/float(im-2) 
      smp=smp+preci(i)/float(im-2)*1000.*3600. 
  31  continue
      do 40 k=1,km
      smc(k)=0.
      sme(k)=0.
      smd(k)=0.
      sms(k)=0.
      smf(k)=0.
      do 30 i=2,im-1
      smc(k)=smc(k)+con(i,k)/float(im-2) 
      sme(k)=sme(k)+eva(i,k)/float(im-2) 
      smd(k)=smd(k)+dep(i,k)/float(im-2) 
      sms(k)=sms(k)+sub(i,k)/float(im-2) 
      smf(k)=smf(k)+fus(i,k)/float(im-2) 
  30  continue
  40  continue
 
      write(70) smc
      write(70) sme
      write(70) smd
      write(70) sms
      write(70) smf
      write(71) smh,smt,smq,smp
!      write(73,194) smc, sme, smd, sms, smf
!      write(74,*) smh,smt,smq,smp
  194 format(8e12.4)

c     sensible and latent heat flux (smt and smq  W/m**2)
c     precipitation (smp  mm/h)
c     pbl height (smh  m)
c     con, eva, sub and fus (smc, sme, smd, sms and smf  K/day)

      write(60) preci,pbl 
      write(60) fsh,flh
      write(61) (upsw(j,1),j=1,mym),
     .          (dnsw(j,1),j=1,mym),
     .          (uplw(j,1),j=1,mym),
     .          (dnlw(j,1),j=1,mym)
      write(61) (upsw(j,mz),j=1,mym),
     .          (dnsw(j,mz),j=1,mym),
     .          (uplw(j,mz),j=1,mym),
     .          (dnlw(j,mz),j=1,mym)
!      write(63,195) preci,pbl,fsh,flh 
!      write(64,195) (upsw(j,1),j=1,mym),
!     .          (dnsw(j,1),j=1,mym),
!     .          (uplw(j,1),j=1,mym),
!     .          (dnlw(j,1),j=1,mym),
c      above is surface, below is top        
!     .           (upsw(j,mz),j=1,mym),
!     .          (dnsw(j,mz),j=1,mym),
!     .          (uplw(j,mz),j=1,mym),
!     .          (dnlw(j,mz),j=1,mym)       
195   format(7e12.5)
       axx=-1.e10
       axn= 1.e10
       bxx=-1.e10
       bxn= 1.e10
      pre=0.
      blh=0.
      do j=2,MY-1
      pre=pre+preci(j)/float(my-2)
      blh=blh+pbl(j)/float(my-2)
      axx=amax1(axx,preci(j))
      axn=amin1(axn,preci(j))
      bxx=amax1(bxx,pbl(j))
      bxn=amin1(bxn,pbl(j))
      enddo

      print*,'FSH',FSH(50)
      print*,'FLH',FLH(50)
      print*,'PRECI',pre,axx,axn
c     print*,'PBL',blh,bxx,bxn
      print*,'upsw',upsw(50,42),upsw(50,1)
      print*,'dnsw',dnsw(50,42),dnsw(50,1)
      print*,'uplw',uplw(50,42),uplw(50,1)
      print*,'dnlw',dnlw(50,42),dnlw(50,1)
 8888 CONTINUE

c     RETURN
      STOP 
  500 CONTINUE
      write(6, 501)
  501 FORMAT('ERROR ABORT.EOF OR EOD ENCOUNTERED ON READ')
  510 CONTINUE
      write(6, 511)
  511 FORMAT('ERROR ABORT OR HARDWARE ERROR ENCOUNTERED ON READ')
c     RETURN
      STOP 
      END
