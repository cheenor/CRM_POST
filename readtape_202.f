      program read_tape 
      PARAMETER (MY=202,MZ=52,mym=my-1, im=202,km=52,nt=2880)
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
      character casenm*10
      casenm='ETP2D1'
      open(70,file=
     *'postdata/micro_202_'//trim(casenm) 
     *,FORM='UNFORMATTED')
      open(71,file=
     *'postdata/micro_202_'//trim(casenm)//'.txt'
     *)
      NVLT=1
      DO 8888 NIOSTR=1,2880
c     DO 8888 NIOSTR=1,2592
      it=NIOSTR
      write(IFNAME,27)JTAPE,NVLT,NIOSTR
      IU=10
      JU=47
      KU=49
 27   FORMAT(A6,'m',I2.2,'nio',I4.4)
      write(WU_IFNAME,28)IFNAME
 28   FORMAT('WUTAPE_',A16)
      write(WU2_IFNAME,29)IFNAME
 29   FORMAT('RADTAPE_',A16)
      write(6,*)'IFNAME  ',IFNAME
      OPEN(UNIT=IU,FILE=IFNAME,STATUS='OLD',FORM='UNFORMATTED'
     *, convert='big_endian' )
      OPEN(JU,FILE=WU_IFNAME,STATUS='OLD',FORM='UNFORMATTED'
     *, convert='big_endian' )
      OPEN(KU,FILE=WU2_IFNAME,STATUS='OLD',FORM='UNFORMATTED'
     *, convert='big_endian' )
C
      READ(JU,IOSTAT=IOVAL) (pbl(J),J=1,MY)
      IF (IOVAL) 500,807,510
  807 CONTINUE
      READ(JU,IOSTAT=IOVAL) (FSH(J),J=1,MY)
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

      write(70) con
      write(70) eva
      write(70) dep
      write(70) sub
      write(70) fus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do ix=1,im 
      write(71,99) (con(ix,iz),iz=1,km)
      write(71,99) (eva(ix,iz),iz=1,km)
      write(71,99) (dep(ix,iz),iz=1,km)
      write(71,99) (sub(ix,iz),iz=1,km)
      write(71,99) (fus(ix,iz),iz=1,km)
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     sensible and latent heat flux (smt and smq  W/m**2)
c     precipitation (smp  mm/h)
c     pbl height (smh  m)
c     con, eva, sub and fus (smc, sme, smd, sms and smf  K/day)

 8888 CONTINUE
99    format(8e12.4)
c     RETURN
      stop 
  500 CONTINUE
      write(6, 501)
  501 FORMAT('ERROR ABORT.EOF OR EOD ENCOUNTERED ON READ')
  510 CONTINUE
      write(6, 511)
  511 FORMAT('ERROR ABORT OR HARDWARE ERROR ENCOUNTERED ON READ')
c     RETURN
      stop 
      END
