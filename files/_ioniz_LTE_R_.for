      PROGRAM CALCULHZ
      USE module_TxtRead

      COMMON /PART/ GAMD,DZETAD,TKPD
C      COMMON /BBB/ BB,TB,PB16,DLCM
C      DIMENSION UPLD(NH),PLELD(NH),TEMPD(NH),UCKD(NH),ECED(NH)
C      DIMENSION UPLD(NH),PLELD(NH),FOTIND(NH),FOTPED(NH)
C      DIMENSION UHMD(NH),BIND(NH),APEKD(NH),DEDXD(NH),RED(NH),ALD(NH)
C      DIMENSION UHMD(NH),BIND(NH),APEKD(NH),ALD(NH)
C      COMMON /MAC/ UPLD,PLELD,TEMPD,UHMD,BIND,APEKD
      COMMON /MAC/ UPLD,PLELD,UHMD,BIND,APEKD
      COMMON /MACF/ FOTIND,FOTPED
C      COMMON /MACT/ UCKD,DEDXD,ECED
      COMMON /VIC/ ALD
      COMMON /ABCYM/ EX,ACYM,BCYM,E1,VF
      COMMON /DHD/ DUHMD                                            ! <!> KV
      COMMON /NUMITER/ NITER

      EXTERNAL FEX,JEX,FEXT,FEXH,JEXT
      DOUBLE PRECISION ATOL,RWORK,RTOL,TD,TOUT,Y,TKPD,DUHMD
      DOUBLE PRECISION SRE1D,SRE2D,HXD,DTD,BBD,GAMD,DZETAD
      DOUBLE PRECISION PLELD,UHMD,UPLD,BIND,APEKD,ALD
      DIMENSION Y(NZ),ATOL(NZ),RWORK(1335),IWORK(121),DUHMD(NZ)
      DIMENSION UPLD(NZ),PLELD(NZ),UHMD(NZ),ALD(NZ)
      DIMENSION BIND(NZ),APEKD(NZ),FOTIND(NZ),FOTPED(NZ)
      DIMENSION UPL(NZ),PLEL(NZ),TEMP(NZ),UHM(NZ),DUHM(NZ)
      ! Veniamin Konovalov addon
      DOUBLE PRECISION :: FOTIND,FOTPED                             ! <!> KV
      CHARACTER(128)   :: CaseDir,TDPDir,sbuf                       ! <!> KV
      CHARACTER(128)   :: InputDataFile, InputGeomFile              ! <!> KV
      INTEGER          :: NT,NSTEP
      INTEGER          :: RunRadEvery,RadCounter
      INTEGER          :: omp_get_max_threads


      ! get "running case" directory name, where input data will be taken and computed data will be saved
      ! ... 1. get from program arguments
      !CALL GETARG(1,CaseDir)
      ! ... 2. get from text file 'run_case.txt'
      open(1,file='run_case.txt')
            read(1,*) CaseDir
            write(*,*) trim(CaseDir)
      close(1)


      ! ... initialize radiation-kinetics model: read NT, NSTEP; create mesh, arrays for radiation...
      call init_RKModel(CaseDir,TDPDir,NT,NSTEP,RunRadEvery)
      CU2D = 0.e0
      WR2D = 0.e0
      WZ2D = 0.e0
      !call close_mod_RKModel()

      ! ... read input data files paths: calul.dat, mpk_elec.dat
      open(1,file=trim(CaseDir)//'config.txt')
            call txt_read_value(1,'InputDataFile',InputDataFile,ios)
            call txt_read_value(1,'InputGeomFile',InputGeomFile,ios)
      close(1)

      PI=3.14159265

      GAM=5./3.
      GAM1=GAM-1.
      FHZ0=0.1
      FHZ0=0.

C      FHZ0=0.1

C      DZETA=0.

      ALKPIZ=0.001
      ALKPIZ=0.01

C       TKPIT=50.
C      ALKPIT=1.E-10

       TKPIT=70.
      ALKPIT=1.E-8

      TOKJFI=0.

C      PRINT 8171
C 8171 FORMAT(1X,'control point')
C      PAUSE



*      OPEN(4,FILE=trim(CaseDir)//'mpk_elec.dat')
      OPEN(4,FILE=trim(CaseDir)//trim(InputGeomFile))
*      READ(4,*) R1N,R2N,DLCM,R0,PB16,EXOL      ! <!> commented KV
      READ(4,*) R1,R2,DLCM,R0,PB16,EXOL
      CLOSE(4)                                  ! <!> KV

C****************
      R0=0.33333
C      PB16=0.081
C      PB16=0.1
      DLCM=60.
      DLCM=10.

      DZETA=2.
      DZETA=1.

C      TB=11600.
C      TB=TB*2.
C      TB=TB*5.
C      TB=4000.
      TB=6000.
      TB=2000.
C      TB=3000.
      TB=750.



C      TB=11600.
C      TB=TB*2.

      TKP=2.16/1.38/TB*10**5

C*******************************
C      PTOP=50.
      PTOP=10.
C      PB16=1.33/1.38*PTOP*1000./TB
C*******************************

C      PB16=25.
      PB16=2.5
      PB16=0.1
      PB16=1.
C      PB16=35.
      PB16=40.


      CAX=1./PB16*(9.1*1.38*TB/2./PI/1.05**2)**1.5/10.

      TOKBKA=300.
      TOKBKA=40.
      TOKBKA=50.

      PM=1.
C      PM=131.

C      DLCM=12.
C*****************

      RXAP=R0
      DL=DLCM*10.
      B=2.*PI*1.38*PB16*TB*RXAP**2*DL**2/10**6/TOKBKA**2

      SIT=1./PB16*(9.1*1.38*TB/2./PI/1.05**2)**1.5*
     *EXP(-2.16*10**5/1.38/TB)/10.
      AION=-SIT/2.E0+SQRT((SIT/2.E0)**2+SIT)
C      AION=1.E-10

      TKPTEM=TKP
      IF(TKPTEM.GT.TKPIT) TKPTEM=TKPIT
C      IF(TKPTEM.GT.TKPIT) GO TO 5757

      TAS=CAX*EXP(-TKPTEM)
      TASP=TAS
      QTAS=SQRT(TASP)
      AION2=-TASP/2.E0+QTAS*SQRT(TASP/4.E0+1.E0)

      AION3=SQRT(TAS/(1.E0+TAS))

      SEF=120.
      SKL=7.
      SK7=1.E-7
C      SRE1=4./9.*4.8**2*SQRT(TB/1.67/9.1/PM/B)*DL*PI/SEF
C      SRE2=SK7/3.*(1.38/4.8)**2*SQRT(PI/1.67/9.1/PM/B)*
C     *DL*TB**2/SKL

C      SRE0=AION/(1.-AION)*4./9.*4.8**2*SQRT(TB/1.67/9.1/PM/B)*
C     *DL*PI/SEF
C      SRE1=SRE0/AION*(1.-AION)

      SRE1=4./9.*4.8**2*SQRT(TB/1.67/9.1/PM/B)*DL*PI/SEF

      SRE2=SK7/3.*(1.38/4.8)**2*SQRT(PI/1.67/9.1/PM/B)*
     *DL*TB**2/SKL
C      EXOL=3./4.8/DLCM*SQRT(1.67/4./PI/PB16)
C**************************
      EXOL=0.
C**************************
      WRITE(*,312) TKP,PB16,CAX,AION,B,SRE1,SRE2,EXOL
  312 FORMAT(1X,'TKP=',F6.2,2X,'PB16=',F6.2,2X,'CAX=',F10.2,
     *2X,'AION=',E10.4,2X,'B=',F6.2,2X,'SRE1=',F10.2,2X,'SRE2=',F6.2,
     *2X,'EXOL=',E10.4)
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)
C      RES=1./SRE0

      WRITE(*,7312) AION,AION2,AION3
 7312 FORMAT(1X,'AION=',E10.4,2X,'AION2=',E10.4,2X,'AION3=',E10.4)
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)



      B100=B*100.
      TB104=TB/10**4
      FQB2=(B100/1.38/TB104)**3
      QB2=SQRT(1.67*PM*FQB2/2.)*DLCM/2./PB16

      FWA2=(B/1.38/TB)**3
      WA2=SQRT(1.67*PM*FWA2/2.)/2./PB16/10**4


      HXAP=200.*TOKBKA/RXAP/DLCM

      VALF=HXAP*1.E+4/SQRT(4.*PI*PM*1.67*PB16)
      VALF1=1.E+6*SQRT(2.*1.38*TB/PM/1.67/B/10**4)
      VALF0=VALF/1.E+6

      WRITE(*,7392) VALF,VALF1,HXAP
 7392 FORMAT(1X,'VALF(cm/c)=',E10.4,2X,'VALF1(cm/c)=',E10.4,
     *2X,'HXAP(ersted)=',F10.2)

      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      VF=TOKBKA/RXAP/DL/SQRT(PI*1.67*PB16)

C      HXAP=200.*TOKBKA/RXAP/DLCM

C      VALF=HXAP*1.E+4/SQRT(4.*PI*PM*1.67*PB16)

      E1=13.6
      E2=E1/4.
      EX=1.6*10**4/1.38/TB
      ACYM=0.10632*DLCM*PB16**2/VF/100.
      BCYM=0.10632*6.0275*DLCM*PB16*10**3/VF/2.


      TKP=2.16/1.38/TB*10**5
      ALF=1.
      BB=B
C****************
      BBK=1.
C      BBK=0.
C****************
      BB2=BB/2.

C      E1=EXOL
      EXOL=0.
C      E1=0.
C*******************
       DZ=1./(NZ-1.)
C*********************************!!!!!!!!!!!!!!!!!
C       DZ=1./(81-1.)

       DY=1./(NR-1.)
      DZ2=DZ*2.
      DY2=DY*2.

      NZM1=NZ-1



C      LKP2=77

C      LKP2=309
C      LKP2=303
C      LKP2=153
      LKP2=77

      LKP2=102

      LKP12=LKP2

      LKP2M1=LKP2-1

      LKP2P1=LKP2+1

C      PRINT 707,LKP2
C  707 FORMAT(1X,'LKP2=',I4)
C      PAUSE

C      DO 708 L=1,80                   !LKP2M1
C      R2(4*L-3)=R2N(L)
C      R2(4*L-2)=0.75*R2N(L)+0.25*R2N(L+1)
C      R2(4*L-1)=0.5*R2N(L)+0.5*R2N(L+1)
C      R2(4*L)=0.25*R2N(L)+0.75*R2N(L+1)

C  708 CONTINUE

      LST=40
C      LST=35
      LST=1

       LSTGEO=34
       LSTGEO=30

      DO 708 L=1,71
C      DO 708 L=1,67
C      R2(2*L-1)=R2N(L)
C      R2(L+LSTGEO)=R2N(L)      ! <!> commented KV
C      R2(L)=R2N(L)

  708 CONTINUE

      DO 799 L=1,LSTGEO
C      R2(L)=R2N(1)-0.0002
  799 CONTINUE

C      DO 7999 L=LKP2,NZ
C      R2(L)=0.0
C**************!!!!!!!!!!!!!!!!!!!!!!!!!!????????????????
 7999 CONTINUE

C      APR=2.
C      ZKR=0.8
C      BPR=ZKR
C      CPR=R0

      DO 711 L=1,NZ
C      R1(L)=R0                  ! <!> commented KV
C      ZZ=(L-1.)*DZ
C      IF(ZZ.GE.ZKR) R1(L)=APR*(ZZ-BPR)**2+CPR

  711 CONTINUE

      OPEN(3,FILE=trim(CaseDir)//'par_mpk.dat')
      WRITE(3,*) R1,R2,DLCM,R0,PB16,EXOL
      CLOSE(3)                             ! <!> KV

C      WRITE(3,*) R1,R2

      PRINT 101
  101 FORMAT(7H R1(L)=,15X,6HR2(L)=,15X,5HZ(L)=)
C      DO 103 L=1,NZ,10
      DO 103 L=1,NZ
       ZZ=(L-1.)*DZ
      PRINT 102,R1(L),R2(L),ZZ,L
  102 FORMAT(F15.6,5X,F15.6,5X,F15.6,5X,I3)
  103 CONTINUE
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)


      DO 1 L=2,NZ
      DR1(L)=(R1(L)-R1(L-1))/DZ
    1 DR2(L)=(R2(L)-R2(L-1))/DZ
      DR1(1)=DR1(2)
      DR2(1)=DR2(2)
      NZM1=NZ-1
C      DO 1 L=2,NZM1
C      DR1(L)=(R1(L+1)-R1(L-1))/DZ2
C    1 DR2(L)=(R2(L+1)-R2(L-1))/DZ2
c      DR1(1)=(R1(2)-R1(1))/DZ
c      DR2(1)=(R2(2)-R2(1))/DZ
C      DR1(NZ)=(R1(NZ)-R1(NZM1))/DZ
C      DR2(NZ)=(R2(NZ)-R2(NZM1))/DZ

      DOS=DY/2.

      DO 2 L=1,NZ
      RDY(L)=R1(L)-R2(L)
      DO 2 M=1,NR
C      YYY=(M-1.)*DY+DOS
      YY=(M-1.)*DY-DY+DOS
      RAD(L,M)=(1.-YY)*R2(L)+YY*R1(L)
C      RADS(L,M)=RAD(L,M)
    2      RDZ(L,M)=(1.-YY)*DR2(L)+YY*DR1(L)
      NZM1=NZ-1
      NRM1=NR-1
      NZM2=NZ-2

      RXAP=R1(1)


C      DO 1799 L=LKP2,NZ
C      RAD(L,1)=RAD(L,2)
C      RAD(L,1)=RAD(L,2)
C      RDZ(L,1)=RDZ(L,2)
C 1799 CONTINUE

      PRINT 8101
 8101 FORMAT(7X,'L=',15X,' RAD(L,1)=')
      DO 8103 L=1,NZ,5
C      DO 8103 L=1,322
       ZZ=(L-1.)*DZ
      PRINT 8102,L,RAD(L,1)
 8102 FORMAT(7X,I4,15X,F15.6)
 8103 CONTINUE
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      DO 260 M=1,NR
        YY=(M-1.)*DY
        YY=(M-1.)*DY
        RAD1M=(1.-YY)*R2(1)+YY*R1(1)
        PBX(M)=RXAP**2/RAD1M**2
        PBX(M)=1.
  260 CONTINUE

      PRINT 891
  891 FORMAT(' YY',15X,'PBX(M)')
      DO 193 M=1,NR,4
       YY=(M-1.)*DY
      PRINT 192,YY,PBX(M)
  192 FORMAT(F15.6,5X,F15.6)
  193 CONTINUE
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      LKP2P1=LKP2+1

      LSTBEG=40
C      LSTBEG=35

      LNEW=LKP2-1
      LNEW=NZ
      DO 3447 M=1,NR
      DO 3447 L=LSTBEG,NZ
      PCI(L,M)=FHZ0*RAD(L,M)**2/2.
      VIF(L,M)=0.
C      FHZ(L,M)=FHZ0
      FHZ(L,M)=0.
      FHR(L,M)=0.
      Z=(L-1.)*DZ
C      Z=Z*2.
C      HRFI(L,M)=-RXAP*((1./10.-1.)*Z+1.)
      HRFI(L,M)=-RXAP*((1./30.-1.)*(Z-DZ*LSTBEG)+1.)
C       VZBX=0.3
       VZBX=0.03
      VIZ(L,M)=VZBX*(60.*(Z-DZ*LSTBEG)+1.)
      PLT(L,M)=PBX(M)*((1./10.-1.)*(Z-DZ*LSTBEG)+1.)
C      PLT(L,M)=PBX(M)
      TEM(L,M)=PLT(L,M)**(GAM-1.)
      TEM(L,M)=1.+Z-DZ*LSTBEG
      VIR(L,M)=RDZ(L,M)*VIZ(L,M)

      TAS=CAX*TEM(L,M)**1.5*EXP(-TKP/TEM(L,M))
      TASP=TAS/PLT(L,M)
      QTAS=SQRT(TASP)
      AL(L,M)=-TASP/2.+QTAS*SQRT(TASP/4.+1.)
C**********
C      AL(L,M)=1.
C***********
      DAB(L,M)=BB/2.*TEM(L,M)*PLT(L,M)*(1.+AL(L,M))
      SS(L,M)=ALOG(DAB(L,M)/(PLT(L,M)**GAM))

 3447 CONTINUE

C      LNEW=90
C      LNEW1=LNEW+1
      DO 710 M=1,NR
      DO 710 L=1,LSTBEG
      PCI(L,M)=FHZ0*RAD(L,M)**2/2.
      VIF(L,M)=0.
C      FHZ(L,M)=FHZ0
      FHZ(L,M)=0.
      FHR(L,M)=0.

C      DO 710 L=LNEW1,NZ
C      PCI(L,M)=PCI(LNEW,M)
C      VIF(L,M)=VIF(LNEW,M)
C      FHZ(L,M)=FHZ(LNEW,M)
C      FHR(L,M)=FHR(LNEW,M)
      HRFI(L,M)=-RXAP
      VIZ(L,M)=VZBX
      PLT(L,M)=PBX(M)
      TEM(L,M)=PLT(L,M)**(GAM-1.)
      VIR(L,M)=RDZ(L,M)*VIZ(L,M)

      TAS=CAX*TEM(L,M)**1.5*EXP(-TKP/TEM(L,M))
      TASP=TAS/PLT(L,M)
      QTAS=SQRT(TASP)
      AL(L,M)=-TASP/2.E0+QTAS*SQRT(TASP/4.E0+1.E0)
C**********
C      AL(L,M)=1.
C***********
      DAB(L,M)=BB/2.*TEM(L,M)*PLT(L,M)*(1.E0+AL(L,M))
      SS(L,M)=ALOG(DAB(L,M)/(PLT(L,M)**GAM))

C      PLT(L,M)=PLT(LNEW,M)
C      TEM(L,M)=TEM(LNEW,M)
C      VIR(L,M)=VIR(LNEW,M)
C      VIR(L,M)=0.
C      HRFI(L,M)=0.
  710 CONTINUE



C        LNEW=61
C      LNEW1=LNEW+1
C      DO 5710 M=1,NR
C      DO 5710 L=LNEW1,NZ
C      PCI(L,M)=PCI(LNEW,M)
C      VIF(L,M)=VIF(LNEW,M)
C      FHZ(L,M)=FHZ(LNEW,M)
C      FHR(L,M)=FHR(LNEW,M)
C      HRFI(L,M)=HRFI(LNEW,M)
C      VIZ(L,M)=VIZ(LNEW,M)
C      PLT(L,M)=PLT(LNEW,M)
C      TEM(L,M)=TEM(LNEW,M)
C      VIR(L,M)=VIR(LNEW,M)
C      VIR(L,M)=0.
C      HRFI(L,M)=0.
C 5710 CONTINUE



      OPEN(1,FILE=trim(CaseDir)//InputDataFile)
CC      READ(1,*) PCIN,HRFIN,FHZN,FHRN,PLTN,TEMN,VIFN,VIRN,VIZN
      READ(1,*) PCI,HRFI,FHZ,FHR,PLT,TEM,VIF,VIR,VIZ,SS,AL
      CLOSE(1)

      OPEN(1,FILE=trim(CaseDir)//'calul.dat')
      REWIND 1



      WRITE(*,104) TB,PB16,DLCM,R0,B,SRE1,SRE2,EXOL,FHZ0,TOKBKA
  104 FORMAT(1X,'TB=',F8.1,2X,'PB16=',F7.3,2X,'DL(CM)=',F6.1,
     *2X,'RXAP=',F7.3,2X,'B=',F7.3,2X,'SRE1=',F10.3,2X,'SRE2=',F10.3,
     *2X,'EXOL=',E10.4,2X,'HZ0=',F6.2,2X,'J(KA)=',F6.1)
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)


 8888 CONTINUE


      PRINT 108
  108 FORMAT('          MACCIB PLT')
      DO 105 L=1,NZ
      DO 105 M=1,NR
      FOHZ(L,M)=FHZ(L,M)
      FOHR(L,M)=FHR(L,M)
      FOPCI(L,M)=PCI(L,M)
      FOTEM(L,M)=TEM(L,M)
      FOVZ(L,M)=VIZ(L,M)
      FOVR(L,M)=VIR(L,M)
      FOVF(L,M)=VIF(L,M)
      FOPL(L,M)=PLT(L,M)
      FORH(L,M)=HRFI(L,M)
      FF(L,M)=PLT(L,M)
C      SS(L,M)=ALOG(TEM(L,M)/(PLT(L,M)**GAM1))
C      SS(L,M)=0.
C      TAS=CAX*TEM(L,M)**1.5*EXP(-TKP/TEM(L,M))
C      TASP=TAS/PLT(L,M)
C      QTAS=SQRT(TASP)
C      AL(L,M)=-TASP/2.+QTAS*SQRT(TASP/4.+1.)
CC**********
CC      AL(L,M)=1.
CC***********
C      DAB(L,M)=BB/2.*TEM(L,M)*PLT(L,M)*(1.+AL(L,M))
C      SS(L,M)=ALOG(DAB(L,M)/(PLT(L,M)**GAM))

  105      CONTINUE


      NZK=NZ/10
C      NZK=NZ/5
      DO 106 M=1,NR,4
      MJ=NR-M+1
      PRINT 107,(FF(L,MJ),L=1,NZ,NZK)
  106 CONTINUE
  107 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 328
  328 FORMAT('          MACCIB TEM')
      DO 325 L=1,NZ
      DO 325 M=1,NR
      FF(L,M)=TEM(L,M)
  325      CONTINUE
      NZK=NZ/10
      DO 326 M=1,NR,4
      MJ=NR-M+1
      PRINT 327,(FF(L,MJ),L=1,NZ,NZK)
  326 CONTINUE
  327 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 7328
 7328 FORMAT('          MACCIB AL')
      DO 7325 L=1,NZ
      DO 7325 M=1,NR
      FF(L,M)=AL(L,M)
 7325      CONTINUE
      NZK=NZ/10
      DO 7326 M=1,NR,4
      MJ=NR-M+1
      PRINT 7327,(FF(L,MJ),L=1,NZ,NZK)
 7326 CONTINUE
 7327 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 1328
 1328 FORMAT('          MACCIB SS')
      DO 1325 L=1,NZ
      DO 1325 M=1,NR
      FF(L,M)=SS(L,M)
 1325      CONTINUE
      NZK=NZ/10
      DO 1326 M=1,NR,4
      MJ=NR-M+1
      PRINT 1327,(FF(L,MJ),L=1,NZ,NZK)
 1326 CONTINUE
 1327 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 728
  728 FORMAT('          MACCIB VIR')
      DO 725 L=1,NZ
      DO 725 M=1,NR
      FF(L,M)=VIR(L,M)
  725      CONTINUE
      NZK=NZ/10
      DO 726 M=1,NR,4
      MJ=NR-M+1
      PRINT 727,(FF(L,MJ),L=1,NZ,NZK)
  726 CONTINUE
  727 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 628
  628 FORMAT('          MACCIB VIZ')
      DO 625 L=1,NZ
      DO 625 M=1,NR
      FF(L,M)=VIZ(L,M)
  625      CONTINUE
      NZK=NZ/10
      DO 626 M=1,NR,4
      MJ=NR-M+1
      PRINT 627,(FF(L,MJ),L=1,NZ,NZK)
  626 CONTINUE
  627 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 4218
 4218 FORMAT('          MACCIB HRFI')
      DO 4215 L=1,NZ
      DO 4215 M=1,NR
      FF(L,M)=HRFI(L,M)
 4215      CONTINUE
      NZK=NZ/10
      DO 4216 M=1,NR,4
      MJ=NR-M+1
      PRINT 4217,(FF(L,MJ),L=1,NZ,NZK)
 4216 CONTINUE
 4217 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      REMT1=1./SRE2
      XOTHP1=EXOL/REMT1
        WRITE(*,110) XOTHP1,REMT1
  110 FORMAT(1X,'XOT=',F10.4,12X,'REM=',F10.4)
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      FKYP=0.15

      HX=AMIN1(DZ,DY)
      ALF=1.

      NRM=NR-8

C      PRINT 111
C  111 FORMAT(1X,'HR(1)=',15X,'HR(NR)=',15X,'PLT(NR)=',15X,'L=')
C      DO 114 L=1,NZ,16

C      PRINT 112,HRFI(L,1),HRFI(L,NR),PLT(L,NR),L
C  112 FORMAT(3(5X,F15.6),5X,I4)
C  114 CONTINUE
C      PAUSE



! =============================================================================
! =============================================================================

      ! OpenMP regime will be on the screen:
      sbuf = "disabled *** "
!$    sbuf = "enabled *** "
      N = 1
!$    N = omp_get_max_threads()
      write(*,"(a,a)")
      write(*,"(a,a)",advance='no')  " *** OpenMP is ",trim(sbuf)
      write(*,"(a,i3)") '    MaxNumThreads = ',N

      TIME=0.

      NDT=NT/NSTEP

      ALLOCATE(TEMPNN(NDT),ALNN(NDT),BPEM(NDT))
      ALLOCATE(TEMPN2(NDT),ALN2(NDT),BPEM2(NDT))

      RadCounter = 0
      DO 3 NN=1,NDT
      DO 333 ND=1,NSTEP
      write(*,"('.')",advance="no")


      SUP=0.
      DO 5 L=1,NZ
      DO 5 M=1,NR

      TEM(1,M)=1.
      PLT(1,M)=1.

      TKPTEM=TKP/TEM(L,M)
      IF(TKPTEM.GT.TKPIT) TKPTEM=TKPIT
C      IF(TKPTEM.GT.TKPIT) GO TO 5757

      TAS=CAX*TEM(L,M)**1.5*EXP(-TKPTEM)
      TASP=TAS/PLT(L,M)
      QTAS=SQRT(TASP)
      AL(L,M)=-TASP/2.E0+QTAS*SQRT(TASP/4.E0+1.E0)
C      AL(L,M)=1.
C 5757 CONTINUE

      DAB(L,M)=BB/2.E0*TEM(L,M)*PLT(L,M)*(1.E0+AL(L,M))
      SS(L,M)=ALOG(DAB(L,M)/(PLT(L,M)**GAM))

      CKOP=SQRT(VIF(L,M)**2+VIR(L,M)**2+VIZ(L,M)**2)
      FDAH=(GAM*DAB(L,M)+HFI(L,M)**2)/PLT(L,M)
      SU1=CKOP+SQRT(FDAH)
    5 IF(SUP.LT.SU1) SUP=SU1
      DT=FKYP*HX/SUP
      DT2=DT/2.

      TIME=TIME+DT



      ADVR=0.
      FPC=0.0
      DO 43 L=1,NZ
      DO 43 M=1,NR
      HFI(L,M)=HRFI(L,M)/RAD(L,M)
   43 CONTINUE


      DO 813 L=LKP2,NZ
      HFI(L,1)=0.
c      HFI(L,2)=0.
      HRFI(L,1)=0.
C      HRFI(L,2)=0.

  813 CONTINUE

      ZION=1.

      DO 8181 L=2,NZM1
      DO 8181 M=2,NRM1
C      FGB=2.*(GAM-1.)/BB*RDY(L)
C      FGB=(GAM-1.)/BB*RDY(L)/TEM(L,M)
C      FGB=2.*(GAM-1.)/BB*RDY(L)/TEM(L,M)/(1.+AL(L,M))

      DRHDY=(HRFI(L,M+1)-HRFI(L,M-1))/DY2
      DRHDZ=(HRFI(L+1,M)-HRFI(L-1,M))/DZ2

C      DRHDY=(HRFI(L-1,M+1)+HRFI(L,M+1)+HRFI(L+1,M+1)
C      *-HRFI(L-1,M-1)-HRFI(L,M-1)-HRFI(L+1,M-1))/DY2/3.

C      DRHDZ=(HRFI(L+1,M-1)+HRFI(L+1,M)+HRFI(L+1,M+1)
C      *-HRFI(L-1,M-1)-HRFI(L-1,M)-HRFI(L-1,M+1))/DZ2/3.

      RDZRDY=RDZ(L,M)/RDY(L)
      PR1=(DRHDY/RDY(L))**2
      PR2=(DRHDZ-RDZRDY*DRHDY)**2
      TOKF(L,M)=(PR1+PR2)/RAD(L,M)**2
C      QCB0=FGB*RE(L,M)*FJ2/RAD(L,M)
 8181      CONTINUE

      DO 3181 M=2,NRM1
      TOKF(1,M)=0.
      TOKF(NZ,M)=TOKF(NZM1,M)
 3181      CONTINUE



C      DO 6105 L=1,NZ
C      DO 6105 M=1,NR
C      TAS=CAX*TEM(L,M)**1.5*EXP(-TKP/TEM(L,M))
C      TASP=TAS/PLT(L,M)
C      QTAS=SQRT(TASP)
C      AL(L,M)=-TASP/2.+QTAS*SQRT(TASP/4.+1.)
CC      AL(L,M)=1.

C      DAB(L,M)=BB/2.*TEM(L,M)*PLT(L,M)*(1.+AL(L,M))
C      SS(L,M)=ALOG(DAB(L,M)/(PLT(L,M)**GAM))

 6105      CONTINUE



      DO 18 L=1,NZ
      DO 18 M=1,NR
C      RE1=(1.-ALF)/ALF/SRE1
C      RELM=RE1+1./SRE2/TEM(L,M)**1.5
       ALRE=AL(L,M)
       IF(ALRE.LT.ALKPIT) ALRE=ALKPIT
      RE1=(1.-ALRE)/ALRE/SRE1
      RE(L,M)=RE1+1./SRE2/TEM(L,M)**1.5

C      RELM=ZION/SRE2/TEM(L,M)**1.5
      RELM=ZION*RE(L,M)
C***********************
      RE(L,M)=RELM
C      RE(L,M)=0.
C***********************
C      X=EXOL*H/RELM/ALF/PLT(L,M)
      CNUR(L,M)=1.

      CNUZ(L,M)=1.

      CNURZ(L,M)=0.

      CNUA(L,M)=1.
   18      CONTINUE



C      PRINT 5918
C 5918 FORMAT('          ANOD')

      M=NR

      DO 130 L=2,NZM1
      DZREL=DR1(L)
      DZREL2=1.+DZREL**2

      DHRDZ=(HRFI(L+1,M)-HRFI(L-1,M))/DZ2
      DHR2DZ=(HRFI(L+1,M)**2-HRFI(L-1,M)**2)/DZ2
      DHRDY=(HRFI(L,M)-HRFI(L,M-1))/DY

      RADLM=RAD(L,M)
      RELM=RE(L,M)

      VN2=0.
      PAR1=RDY(L)/(1.+DZREL**2)
      BEL1=DZREL*RELM*DHRDZ

      REDHRY=PAR1*BEL1
      WUP=-DY*REDHRY
      BAN1(L)=WUP

      BEL1=DZREL*DHRDZ

      REDHRY=PAR1*BEL1
      BAN(L)=DY*REDHRY*RDY(L)

      VNAN(L)=VN2

  130 CONTINUE


      BAN1(1)=0.

      BAN1(NZ)=BAN1(NZM1)

      BAN(1)=0.

      BAN(NZ)=BAN(NZM1)

      VNAN(NZ)=VNAN(NZM1)



C      PRINT 6918
C 6918 FORMAT('         KATOD')

      M=1

C      DO 7131 L=LKP2,NZM1
C      VNKA(L)=0.
C 7131 CONTINUE

      DO 131 L=2,NZM1
C      DO 131 L=2,LKP2M1
      DZREL=DR2(L)
      DZREL2=1.+DZREL**2

      DHRDZ=(HRFI(L+1,M)-HRFI(L-1,M))/DZ2
      DHR2DZ=(HRFI(L+1,M)**2-HRFI(L-1,M)**2)/DZ2

      DHRDY=(HRFI(L,M+1)-HRFI(L,M))/DY

      RADLM=RAD(L,M)
      RELM=RE(L,M)


      VN2=0.




      PAR1=RDY(L)/(1.+DZREL**2)
      BEL1=DZREL*RELM*DHRDZ
      REDHRY=PAR1*BEL1
      WUP=-DY*REDHRY
      BKA1(L)=WUP

      BEL1=DZREL*DHRDZ
      REDHRY=PAR1*BEL1
      BKA(L)=-DY*REDHRY*RDY(L)

      VNKA(L)=VN2


  131 CONTINUE


      BKA1(1)=0.

      BKA1(NZ)=BKA1(NZM1)

      BKA(1)=0.

      BKA(NZ)=BKA(NZM1)

      VNKA(NZ)=VNKA(NZM1)



      IF(BBK.EQ.0.) GOTO 214

      ! ... radiation transport:
      RadCounter = RadCounter + 1
      if(RadCounter.eq.RunRadEvery) then
            call rad_transport(TIME,NZ,NR,PLT,TEM,CU2D,WR2D,WZ2D)
            RadCounter = 0
      endif


      DO 181 L=2,NZM1
      DO 181 M=2,NRM1
C      FGB=2.*(GAM-1.)/BB*RDY(L)
      FGB=(GAM-1.)/BB*RDY(L)/TEM(L,M)
      FGB=2.*(GAM-1.)/BB*RDY(L)/TEM(L,M)/(1.+AL(L,M))

      DRHDY=(HRFI(L,M+1)-HRFI(L,M-1))/DY2
      DRHDZ=(HRFI(L+1,M)-HRFI(L-1,M))/DZ2

C      DRHDY=(HRFI(L-1,M+1)+HRFI(L,M+1)+HRFI(L+1,M+1)
C      *-HRFI(L-1,M-1)-HRFI(L,M-1)-HRFI(L+1,M-1))/DY2/3.

C      DRHDZ=(HRFI(L+1,M-1)+HRFI(L+1,M)+HRFI(L+1,M+1)
C      *-HRFI(L-1,M-1)-HRFI(L-1,M)-HRFI(L-1,M+1))/DZ2/3.

      RDZRDY=RDZ(L,M)/RDY(L)
      PR1=(DRHDY/RDY(L))**2
      PR2=(DRHDZ-RDZRDY*DRHDY)**2
      FJ2=PR1+PR2
      QCB0=FGB*RE(L,M)*FJ2/RAD(L,M)

      TEB=TEM(L,M)*TB/1.16/10**4
      FQ=QB2*PB16**2*PLT(L,M)**2/10**5
C      uzing the encyclopedia for T < 4 eV
C      QLIN=6.25*FQ*ZION**4
C      uzing the paper by Kogan
C      QLIN=80.*FQ*ZION**7/SQRT(TEB**3)

      QLIN=80.*FQ*ZION**7/SQRT(TEB**3)*(1.-AL(L,M))

      QPEK=4.4*FQ*ZION**5/SQRT(TEB)
      QTOP=0.154*FQ*ZION**3*SQRT(TEB)

      QFOT=QLIN+QPEK+QTOP
C      QFOT=QPEK+QTOP

C      QCB1=QFOT*FGB*RAD(L,M)
C      QCB1=QFOT*FGB*RAD(L,M)*AL(L,M)

      QCB1=QFOT*FGB*RAD(L,M)*AL(L,M)*AL(L,M)

      DWDY=(WR2D(L,M+1)-WR2D(L,M-1))/DY2
      DWDZ=(WZ2D(L+1,M)-WZ2D(L-1,M))/DZ2
      PW1=DWDY/RDY(L)
      PW2=DWDZ-RDZRDY*DWDY
      DIVW=PW1+PW2
      QCB1=DIVW*WA2*FGB*RAD(L,M)

C      QCB1=0.
C      IF(AL(L,M).LE.ALKPIZ) QCB1=0.

C      NHM1=NH-1
C      DO 768 N=1,NH
      TEMPN=TEM(L,M)
      TEB=TEMPN/EX
      PLOTN=PLT(L,M)

      BPR0=DLCM/VF*PB16*5.6/SQRT(10.)/10**8
      EA1=2.334733
      EA2=0.250621
      EB1=3.330657
      EB2=1.681534
      APR0=0.332462*BPR0/10**5

      XE=E1/TEB
      IF(XE.GT.TKPIT) XE=TKPIT

      E12X=(XE**2+EA1*XE+EA2)/
     *(XE**2+EB1*XE+EB2)/XE
      E1X=EXP(-XE)*E12X
      BPRF=SQRT(TEB)*E1X
      BPR=BPRF*BPR0
      APR=APR0*E12X/TEB


      KCYM=20
      CYM=0.
      DO 1517 K=1,KCYM

      FLK=0.3*ALOG(1.017+0.462*TEB/E1*K**2*(K+1.)**2/(2.*K+1.))

      X=TEB*K**2*(K+1)**2/E1/(2*K+1)
      CYM=CYM+2.*(2*K+1)*EXP(-E1/TEB/(K+1)**2)/K**3/(K+1)**4/FLK
C      /S(X)
 1517 CONTINUE

      APEKN=ACYM/CYM/TEB**2+APR

      BINN=BCYM*EXP(-XE)/SQRT(TEB)/CYM+BPR
C        BIND(N)=DBLE(BINN)
C        APEKD(N)=DBLE(APEKN)
      FTINN=0.312*DLCM/VF*E1**2*EXP(-XE)*TEB
      FTPEN=6.457*DLCM*PB16/VF*E1**2*SQRT(TEB)/10**7
C      FOTIND(N)=DBLE(FTINN)
C        FOTPED(N)=DBLE(FTPEN)

      PLELN=AL(L,M)*PLT(L,M)

C  768 CONTINUE


C      DO 1 N=2,NHM1
C      T3=Y(N)
C      RE1=(UPLD(N)-PLELD(N))/PLELD(N)/SRE1D
C      RE2=1.D0/SRE2D/DSQRT(T3*T3*T3)
C      RELMD=ZIOND*(RE1+RE2)
      FK1=(PLOTN-PLELN)*PLELN*BINN
      FK2=PLELN*PLELN*PLELN*APEKN
      FK3=(PLOTN-PLELN)*FTINN
      FK4=PLELN*PLELN*FTPEN



C      PART2=0.D0
C      PART2=2.D0/BBD/UPLD(N)/(1.D0+ALD(N))*(GAMD-1.D0)*QFOTD(N)
C      PART3=2.D0/BBD/UPLD(N)/(1.D0+ALD(N))*(GAMD-1.D0)*RELMD*DUHMD(N)
C      PART2=0.D0
C      PART3=0.D0

C      YDOT(N)=PART3-(FK1-FK2+FK3-FK4)*(GAMD-1.D0)*DZETAD*TKPD/UPLD(N)/
C      *(1.D0+ALD(N))-PART2

      QCB2=FGB*RAD(L,M)*(FK1-FK2+FK3-FK4)*DZETA*BB/2.*TKP

C      QCB2=0.

      EQ(L,M)=QCB0-QCB1-QCB2
C      EQ(L,M)=0.


C      IF(RE(1,1).EQ.0.) GOTO 181
C      DHDR12=(DHDRLM(L,M+1)+DHDRLM(L,M))/2.
C      DHDZ12=(DHDZLM(L+1,M)+DHDZLM(L,M))/2.
C      PR1=(DHDRLM(L,M+1)/RDY(L))**2
C      PR2=(DHDZLM(L+1,M)-RDZRDY*DHDRLM(L,M+1))**2
C      PR1=(DHDR12/RDY(L))**2
C      PR2=(DHDZ12-RDZRDY*DHDR12)**2

C      FJ2=PR1+PR2
C      QCB0=FGB*RE(L,M)*FJ2/RAD(L,M)

C      EQ(L,M)=QCB0

  181 CONTINUE

      M=2
C      DO 9181 L=LKP2,NZM1
C      FGB=2.*(GAM-1.)/BB*RDY(L)
C      FGB=(GAM-1.)/BB*RDY(L)/TEM(L,M)
C      DRHDY=(HRFI(L,M+1)-HRFI(L,M))/DY
C      DRHDZ=(HRFI(L+1,M)-HRFI(L-1,M))/DZ2
C      RDZRDY=RDZ(L,M)/RDY(L)
C      PR1=(DRHDY/RDY(L))**2
C      PR2=(DRHDZ-RDZRDY*DRHDY)**2
C      FJ2=PR1
C      QCB0=FGB*RE(L,M)*FJ2/RAD(L,M)



C      EQ(L,M)=QCB0
C      EQ(L,M)=0.

C      IF(RE(1,1).EQ.0.) GOTO 9181
C      DHDR12=(DHDRLM(L,M+1)+DHDRLM(L,M))/2.

C      PR1=(DHDRLM(L,M+1)/RDY(L))**2
C      PR1=(DHDR12/RDY(L))**2
C      FJ2=PR1
C      QCB0=FGB*RE(L,M)*FJ2/RAD(L,M)

C      EQ(L,M)=QCB0
C      EQ(L,M)=0.


 9181 CONTINUE

      DO 182 L=2,NZM1
      EQ(L,NR)=EQ(L,NRM1)
C      EQ(L,2)=EQ(L,3)
      EQ(L,1)=EQ(L,2)

  182 CONTINUE
      DO 183 M=1,NR
      EQ(1,M)=EQ(2,M)
      EQ(NZ,M)=EQ(NZM1,M)

  183 CONTINUE


  214 CONTINUE




C      DO 6 L=2,NZ
      DO 6 L=2,NZM1
      DO 7 M=1,NR
      PTI(M)=RAD(L,M)*RDY(L)*PLT(L,M)
      UCK(M)=(VIR(L,M)-RDZ(L,M)*VIZ(L,M))/RDY(L)
      D1(M)=0.
      D2(M)=0.
      D3(M)=0.
      C2(M)=0.

    7      P12DT(M)=PTI(M)
      VN=VNAN(L)
      UCK(NR)=VN/RDY(L)
C      UCK(NR)=UCK(NRM1)
C      IF(L.GE.LKP2) UCK(NR)=UCK(NRM1)
C      V1=VNKA(L)
C      UCK(1)=V1/RDY(L)
C      UCK(1)=UCK(2)
      IF(L.LT.LKP2) UCK(1)=UCK(2)

      IF(L.GE.LKP2) UCK(1)=-UCK(2)
C      IF(L.GE.LKP12) UCK(1)=0.
C      IF(L.GE.LKP12) UCK(1)=-UCK(2)
C      IF(L.GE.LKP12) UCK(2)=0.
      PPL=0.
      PLB=1.


      PRB=1.
      PPR=0.

      PLB=RAD(L,1)/RAD(L,2)

      PRB=RAD(L,NR)/RAD(L,NRM1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NR,DY,DT2,UCK,D1,D2,C2,D3)


      DO 8 M=1,NR
    8 PLT12(L,M)=PDTI(M)/RAD(L,M)/RDY(L)
    6      CONTINUE

      DO 46 L=2,NZM1

      DO 47 M=1,NR
      PTI(M)=RAD(L,M)*RDY(L)*PLT(L,M)*VIR(L,M)
      UCK(M)=(VIR(L,M)-RDZ(L,M)*VIZ(L,M))/RDY(L)
      DABLM=BB/2.*TEM(L,M)*PLT(L,M)*(1.+AL(L,M))
      DABLM=DAB(L,M)
      D1(M)=-RAD(L,M)*(DABLM+HFI(L,M)**2/2.)
      D2(M)=0.
      C2(M)=0.
      D3(M)=RDY(L)*(DABLM-HFI(L,M)**2/2.)

C      D1(M)=0.

C      D2(M)=DABLM
C      C2(M)=-RAD(L,M)

C      D22(M)=HRFI(L,M)**2
C      C22(M)=-1./2./RAD(L,M)

C      D22(M)=HRFI(L,M)
C      C22(M)=-HFI(L,M)

C      D3(M)=0.

   47      P12DT(M)=PTI(M)
      VN=VNAN(L)
      UCK(NR)=VN/RDY(L)
C      UCK(NR)=UCK(NRM1)
C      IF(L.GE.LKP2) UCK(NR)=UCK(NRM1)
      V1=VNKA(L)
C      UCK(1)=V1/RDY(L)
C      UCK(1)=UCK(2)
      IF(L.LT.LKP2) UCK(1)=UCK(2)

C      UCK(1)=0.
C      IF(L.GE.LKP2) UCK(1)=0.
C      IF(L.GE.LKP12) UCK(1)=-UCK(2)

C      IF(L.GE.LKP2) UCK(2)=0.
      IF(L.GE.LKP2) UCK(1)=-UCK(2)

      VR1=DR2(L)*VIZ(L,1)+V1
      VRN=DR1(L)*VIZ(L,NR)+VN


      PLB=0.
      PPL=VR1*RAD(L,1)*RDY(L)*PLT(L,1)

C      IF(L.GE.LKP2) D1(1)=0.
C      IF(L.GE.LKP2) D1(2)=0.
C      IF(L.GE.LKP2) D3(1)=0.
C      IF(L.GE.LKP2) D3(2)=0.
C      IF(L.GE.LKP2) C2(1)=0.
C      IF(L.GE.LKP2) C2(2)=0.
C      IF(L.GE.LKP2) C22(1)=0.
C      IF(L.GE.LKP2) C22(2)=0.
C      IF(L.GE.LKP2) C2(3)=0.
C      IF(L.GE.LKP2) C22(3)=0.

      IF(L.GE.LKP2) PPL=0.
      IF(L.GE.LKP2) PLB=-1.
C      IF(L.GE.LKP2) PPL=-RAD(L,2)*RDY(L)*PLT(L,2)*VIR(L,2)

      PRB=0.
      PPR=VRN*RAD(L,NR)*RDY(L)*PLT(L,NR)

C      CALL LCP2(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
C     *NR,DY,DT2,UCK,D1,D2,C2,D3,D22,C22)
      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NR,DY,DT2,UCK,D1,D2,C2,D3)

      DO 48 M=1,NR
   48      VIR12(L,M)=PDTI(M)/RAD(L,M)/RDY(L)/PLT12(L,M)
   46      CONTINUE


      DO 80 L=2,NZM1
C      DO 80 L=2,NZ
      DO 81 M=1,NR
      PTI(M)=RAD(L,M)*RDY(L)*PLT(L,M)*VIZ(L,M)
      UCK(M)=(VIR(L,M)-RDZ(L,M)*VIZ(L,M))/RDY(L)
      D1(M)=0.
      DABLM=BB/2.*TEM(L,M)*PLT(L,M)*(1.+AL(L,M))
      DABLM=DAB(L,M)

      D2(M)=DABLM+HFI(L,M)**2/2.

      C2(M)=RAD(L,M)*RDZ(L,M)
      D3(M)=0.
   81      P12DT(M)=PTI(M)
      VN=VNAN(L)
      UCK(NR)=VN/RDY(L)
C      UCK(NR)=UCK(NRM1)
C      IF(L.GE.LKP2) UCK(NR)=UCK(NRM1)
C      V1=VNKA(L)
C      UCK(1)=V1/RDY(L)
C      UCK(1)=UCK(2)
C      IF(L.GE.LKP12) UCK(1)=0.
C      IF(L.GE.LKP12) UCK(1)=-UCK(2)
      IF(L.LT.LKP2) UCK(1)=UCK(2)

C      IF(L.GE.LKP12) UCK(2)=0.
      PPL=0.
      PLB=1.
      IF(L.GE.LKP2) UCK(1)=-UCK(2)
      PRB=1.
      PPR=0.

C      IF(L.GE.LKP2) PPL=RAD(L,1)*RDY(L)*PLT(L,1)*VIZ(L,1)
C      IF(L.GE.LKP2) PLB=0.

C      IF(L.GE.LKP2) PLB=RAD(L,1)/RAD(L,2)
C      IF(L.EQ.LKP2) PLB=0.

      PLB=RAD(L,1)/RAD(L,2)

      PRB=RAD(L,NR)/RAD(L,NRM1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NR,DY,DT2,UCK,D1,D2,C2,D3)

      DO 82 M=1,NR
   82      VIZ12(L,M)=PDTI(M)/RAD(L,M)/RDY(L)/PLT12(L,M)
   80      CONTINUE

   92      CONTINUE

      IF(BBK.EQ.0.) GOTO 201

C      DO 184 L=2,NZ
      DO 184 L=2,NZM1
      DO 185 M=1,NR
      PTI(M)=RAD(L,M)*RDY(L)*PLT(L,M)*SS(L,M)
      UCK(M)=(VIR(L,M)-RDZ(L,M)*VIZ(L,M))/RDY(L)
      D1(M)=0.
      D2(M)=0.
      C2(M)=0.
      D3(M)=0.
      D3(M)=EQ(L,M)/2.
C      D3(M)=EQ(L,M)
  185      P12DT(M)=PTI(M)

      VN=VNAN(L)
      UCK(NR)=VN/RDY(L)
C      UCK(NR)=UCK(NRM1)
C      IF(L.GE.LKP2) UCK(NR)=UCK(NRM1)
C      V1=VNKA(L)
C      UCK(1)=V1/RDY(L)
C      UCK(1)=UCK(2)
C      IF(L.GE.LKP12) UCK(1)=0.
C      IF(L.GE.LKP12) UCK(1)=-UCK(2)
      IF(L.LT.LKP2) UCK(1)=UCK(2)

C      IF(L.GE.LKP12) UCK(2)=0.

      IF(L.GE.LKP2) UCK(1)=-UCK(2)
      PLB=1.
      PPL=0.
      PRB=1.
      PPR=0.
C      IF(L.GE.LKP2) PLB=RAD(L,1)/RAD(L,2)

      PLB=RAD(L,1)/RAD(L,2)

      PRB=RAD(L,NR)/RAD(L,NRM1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NR,DY,DT2,UCK,D1,D2,C2,D3)
      DO 188 M=1,NR
  188      SS12(L,M)=PDTI(M)/RAD(L,M)/RDY(L)/PLT12(L,M)
  184      CONTINUE

  201 CONTINUE


      DO 36 L=2,NZM1
      DO 37 M=1,NR
      PTI(M)=RDY(L)*HRFI(L,M)
      UCK(M)=(VIR(L,M)-RDZ(L,M)*VIZ(L,M))/RDY(L)

      D1(M)=0.

      D2(M)=0.
      C2(M)=0.


      D3(M)=0.
   37      P12DT(M)=PTI(M)


      V1=VNKA(L)
      VN=VNAN(L)
      UCK(NR)=VN/RDY(L)
C      UCK(NR)=UCK(NRM1)

       PVH1=V1*HRFI(L,1)
      PVHN=VN*HRFI(L,NR)

      IF(L.LT.LKP2) UCK(1)=V1/RDY(L)
C      UCK(1)=      UCK(2)
C      IF(L.GE.LKP2) UCK(1)=UCK(2)
C      IF(L.GE.LKP2) UCK(1)=0.
C      IF(L.GE.LKP12) UCK(1)=-UCK(2)

C      IF(L.GE.LKP2) UCK(2)=0.
      IF(L.GE.LKP2) UCK(1)=-UCK(2)

C      PLB=0.
C      PPL=0.

C      PLB=1.
C      PPL=BKA(L)

C      PLB=1.
C      PPL=BKA1(L)*RDY(L)/RE(L,1)

      PLB=0.
      PPL=RDY(L)*HRFI(1,1)

            IF(L.GE.LST) PLB=1.
      IF(L.GE.LST) PPL=BKA1(L)*RDY(L)/RE(L,1)

      IF(L.GE.LKP2) PLB=0.
      IF(L.GE.LKP2) PPL=-RDY(L)*HRFI(L,2)
C      IF(L.GE.LKP2) PPL=0.

C      PRB=1.
C      PPR=BAN(L)

      PRB=1.
      PPR=-BAN1(L)*RDY(L)/RE(L,NR)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NR,DY,DT2,UCK,D1,D2,C2,D3)


       DO 38 M=1,NR
   38      HRFI12(L,M)=PDTI(M)/RDY(L)

   36      CONTINUE

      DO 9 M=1,NR
      HRFI12(1,M)=HRFI(1,M)
      HRFI12(NZ,M)=HRFI(NZ,M)
      PLT12(NZ,M)=PLT(NZ,M)
      VIR12(NZ,M)=VIR(NZ,M)
      VIZ12(NZ,M)=VIZ(NZ,M)
      SS12(NZ,M)=SS(NZ,M)
      SS12(1,M)=SS(1,M)
      PLT12(1,M)=PLT(1,M)
      VIZ12(1,M)=VIZ(1,M)
      VIR12(1,M)=VIR(1,M)

    9 CONTINUE

          DO 22 M=1,NR
      DO 22 L=1,NZ
      HFI12(L,M)=HRFI12(L,M)/RAD(L,M)
      IF(BBK-0.) 202,202,203
  202      TEM12(L,M)=PLT12(L,M)**(GAM-1.)
      GOTO 22
  203      CONTINUE
C      TEM12(L,M)=PLT12(L,M)**GAM1*EXP(SS12(L,M))
      DAB12(L,M)=PLT12(L,M)**GAM*EXP(SS12(L,M))
C            SS(L,M)=ALOG(DAB(L,M)/(PLT(L,M)**GAM))


   22 CONTINUE


      DO 10 L=2,NZM1
C      DO 10 L=2,NZ
      DO 11 M=1,NR
      PTI(M)=RAD(L,M)*RDY(L)*PLT(L,M)
      UCK(M)=(VIR12(L,M)-RDZ(L,M)*VIZ12(L,M))/RDY(L)
      D1(M)=0.
      D2(M)=0.
      D3(M)=0.
      C2(M)=0.

   11 P12DT(M)=RAD(L,M)*RDY(L)*PLT12(L,M)
      VN=VNAN(L)
      UCK(NR)=VN/RDY(L)
C      UCK(NR)=UCK(NRM1)
C      IF(L.GE.LKP2) UCK(NR)=UCK(NRM1)
C      V1=VNKA(L)
C      UCK(1)=V1/RDY(L)
C      UCK(1)=UCK(2)
C      IF(L.GE.LKP12) UCK(1)=0.
C      IF(L.GE.LKP12) UCK(1)=-UCK(2)
      IF(L.LT.LKP2) UCK(1)=UCK(2)

C      IF(L.GE.LKP12) UCK(2)=0.
      IF(L.GE.LKP2) UCK(1)=-UCK(2)

       PPL=0.
      PLB=1.
C      IF(L.GE.LKP2) PLB=RAD(L,1)/RAD(L,2)

      PRB=1.
      PPR=0.

      PLB=RAD(L,1)/RAD(L,2)

      PRB=RAD(L,NR)/RAD(L,NRM1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NR,DY,DT,UCK,D1,D2,C2,D3)

      DO 12 M=1,NR
   12      PLT2(L,M)=PDTI(M)/RAD(L,M)/RDY(L)
   10      CONTINUE

      DO 50 L=2,NZM1

      DO 51 M=1,NR
      PTI(M)=RAD(L,M)*RDY(L)*PLT(L,M)*VIR(L,M)
      UCK(M)=(VIR12(L,M)-RDZ(L,M)*VIZ12(L,M))/RDY(L)
C      DABLM=BB/2.*TEM12(L,M)*PLT12(L,M)*(1.+ALF)
      DABLM=DAB12(L,M)
      D1(M)=-RAD(L,M)*(DABLM+HFI12(L,M)**2/2.)
      D2(M)=0.
      C2(M)=0.


      D3(M)=RDY(L)*(DABLM-HFI12(L,M)**2/2.)

C      D1(M)=0.

C      D2(M)=DABLM
C      C2(M)=-RAD(L,M)

CC      D22(M)=HRFI12(L,M)**2
CC      C22(M)=-1./2./RAD(L,M)

C      D22(M)=HRFI12(L,M)
C      C22(M)=-HFI12(L,M)

C      D3(M)=0.

   51 P12DT(M)=RAD(L,M)*RDY(L)*PLT12(L,M)*VIR12(L,M)
      VN=VNAN(L)
      UCK(NR)=VN/RDY(L)
C      UCK(NR)=UCK(NRM1)
C      IF(L.GE.LKP2) UCK(NR)=UCK(NRM1)
      V1=VNKA(L)
C      UCK(1)=V1/RDY(L)
C      UCK(1)=UCK(2)
      IF(L.LT.LKP2) UCK(1)=UCK(2)

C      IF(L.LT.LKP2) UCK(1)=UCK(2)
C      IF(L.GE.LKP2) UCK(1)=0.
      IF(L.GE.LKP12) UCK(1)=-UCK(2)

C      IF(L.GE.LKP2) UCK(2)=0.
      VR1=DR2(L)*VIZ(L,1)+V1
      VRN=DR1(L)*VIZ(L,NR)+VN

C      IF(L.GE.LKP2) UCK(1)=-UCK(2)

C      IF(L.GE.LKP2) D1(1)=0.
C      IF(L.GE.LKP2) D1(1)=0.
C      IF(L.GE.LKP2) D1(2)=0.
C      IF(L.GE.LKP2) D3(1)=0.
C      IF(L.GE.LKP2) D3(2)=0.
C      IF(L.GE.LKP2) C2(1)=0.
C      IF(L.GE.LKP2) C2(2)=0.
C      IF(L.GE.LKP2) C22(1)=0.
C      IF(L.GE.LKP2) C22(2)=0.
C      IF(L.GE.LKP2) C2(3)=0.
C      IF(L.GE.LKP2) C22(3)=0.

      PLB=0.
      PPL=VR1*RAD(L,1)*RDY(L)*PLT(L,1)
      IF(L.GE.LKP2) PPL=0.
      IF(L.GE.LKP2) PLB=-1.
C      IF(L.GE.LKP2) PPL=-RAD(L,2)*RDY(L)*PLT12(L,2)*VIR12(L,2)

      PRB=0.
      PPR=VRN*RAD(L,NR)*RDY(L)*PLT(L,NR)

C      CALL LCP2(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
C     *NR,DY,DT,UCK,D1,D2,C2,D3,D22,C22)
      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NR,DY,DT,UCK,D1,D2,C2,D3)

      DO 52 M=1,NR
   52      VIR2(L,M)=PDTI(M)/RAD(L,M)/RDY(L)/PLT2(L,M)
   50      CONTINUE

      DO 83 L=2,NZM1
C      DO 83 L=2,NZ
      DO 84 M=1,NR
      PTI(M)=RAD(L,M)*RDY(L)*PLT(L,M)*VIZ(L,M)
      UCK(M)=(VIR12(L,M)-RDZ(L,M)*VIZ12(L,M))/RDY(L)
      D1(M)=0.
C      DABLM=BB/2.*TEM12(L,M)*PLT12(L,M)*(1.+ALF)
      DABLM=DAB12(L,M)

      D2(M)=DABLM+HFI12(L,M)**2/2.

      D3(M)=0.
      C2(M)=RAD(L,M)*RDZ(L,M)
   84 P12DT(M)=RAD(L,M)*RDY(L)*PLT12(L,M)*VIZ12(L,M)
      VN=VNAN(L)
      UCK(NR)=VN/RDY(L)
C      UCK(NR)=UCK(NRM1)
C      IF(L.GE.LKP2) UCK(NR)=UCK(NRM1)
C      V1=VNKA(L)
C      UCK(1)=V1/RDY(L)
C      UCK(1)=UCK(2)
C      IF(L.GE.LKP12) UCK(1)=0.
C      IF(L.GE.LKP12) UCK(1)=-UCK(2)
      IF(L.LT.LKP2) UCK(1)=UCK(2)

C      IF(L.GE.LKP12) UCK(2)=0.
      PPL=0.
      PLB=1.
      PRB=1.
      PPR=0.
      IF(L.GE.LKP2) UCK(1)=-UCK(2)
C      IF(L.EQ.LKP2) PLB=0.

C      IF(L.GE.LKP2) PPL=RAD(L,1)*RDY(L)*PLT(L,1)*VIZ12(L,1)
c      IF(L.GE.LKP2) PLB=0.

      PLB=RAD(L,1)/RAD(L,2)

      PRB=RAD(L,NR)/RAD(L,NRM1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NR,DY,DT,UCK,D1,D2,C2,D3)

      DO 85 M=1,NR
   85      VIZ2(L,M)=PDTI(M)/RAD(L,M)/RDY(L)/PLT2(L,M)
   83      CONTINUE


      IF(BBK.EQ.0.) GOTO 204

C      DO 190 L=2,NZ
      DO 190 L=2,NZM1
      DO 191 M=1,NR
      PTI(M)=RAD(L,M)*RDY(L)*PLT(L,M)*SS(L,M)
      UCK(M)=(VIR12(L,M)-RDZ(L,M)*VIZ12(L,M))/RDY(L)
      D1(M)=0.
      D2(M)=0.
      C2(M)=0.
      D3(M)=0.

      D3(M)=EQ(L,M)/2.
C      D3(M)=EQ(L,M)
  191      P12DT(M)=RAD(L,M)*RDY(L)*PLT12(L,M)*SS12(L,M)

      VN=VNAN(L)
      UCK(NR)=VN/RDY(L)
C      UCK(NR)=UCK(NRM1)
C      IF(L.GE.LKP2) UCK(NR)=UCK(NRM1)
C      V1=VNKA(L)
C      UCK(1)=V1/RDY(L)
C      UCK(1)=UCK(2)
C      IF(L.GE.LKP12) UCK(1)=0.
C      IF(L.GE.LKP12) UCK(1)=-UCK(2)
      IF(L.LT.LKP2) UCK(1)=UCK(2)

C      IF(L.GE.LKP12) UCK(2)=0.
      IF(L.GE.LKP2) UCK(1)=-UCK(2)
      PPL=0.
      PLB=1.
      PRB=1.
      PPR=0.
       PV1=0.

      PLB=RAD(L,1)/RAD(L,2)

      PRB=RAD(L,NR)/RAD(L,NRM1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NR,DY,DT,UCK,D1,D2,C2,D3)

      DO 194 M=1,NR
  194      SS2(L,M)=PDTI(M)/RAD(L,M)/RDY(L)/PLT2(L,M)
  190      CONTINUE

  204 CONTINUE

      DO 40 L=2,NZM1
      DO 41 M=1,NR
      PTI(M)=RDY(L)*HRFI(L,M)
      UCK(M)=(VIR12(L,M)-RDZ(L,M)*VIZ12(L,M))/RDY(L)

      D1(M)=0.

      D2(M)=0.
      C2(M)=0.


      D3(M)=0.
   41 P12DT(M)=RDY(L)*HRFI12(L,M)


      V1=VNKA(L)

      VN=VNAN(L)

      UCK(NR)=VN/RDY(L)
C      UCK(NR)=UCK(NRM1)


       PVH1=V1*HRFI(L,1)
      PVHN=VN*HRFI(L,NR)

C      UCK(1)=V1/RDY(L)
      IF(L.LT.LKP2) UCK(1)=V1/RDY(L)

C      IF(L.GE.LKP2) UCK(1)=UCK(2)
C      IF(L.GE.LKP2) UCK(1)=0.
C      IF(L.GE.LKP12) UCK(1)=-UCK(2)

C      IF(L.GE.LKP2) UCK(2)=0.
      IF(L.GE.LKP2) UCK(1)=-UCK(2)


C      PLB=0.
C      PPL=0.

C      PLB=1.
C      PPL=BKA(L)

C      PLB=1.
C      PPL=BKA1(L)*RDY(L)/RE(L,1)

      PLB=0.
      PPL=RDY(L)*HRFI(1,1)

            IF(L.GE.LST) PLB=1.
      IF(L.GE.LST) PPL=BKA1(L)*RDY(L)/RE(L,1)

      IF(L.GE.LKP2) PLB=0.
      IF(L.GE.LKP2) PPL=-RDY(L)*HRFI12(L,2)
C      IF(L.GE.LKP2) PPL=0.

C      PRB=1.
C      PPR=BAN(L)

      PRB=1.
      PPR=-BAN1(L)*RDY(L)/RE(L,NR)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NR,DY,DT,UCK,D1,D2,C2,D3)


      DO 42 M=1,NR
   42      HRFI2(L,M)=PDTI(M)/RDY(L)
   40      CONTINUE

      DO 129 M=1,NR
C      HRFI2(1,M)=HRFI(1,M)
      HRFI2(NZ,M)=HRFI(NZ,M)
      SS2(NZ,M)=SS(NZM1,M)
      PLT2(NZ,M)=PLT(NZM1,M)
      VIR2(NZ,M)=VIR(NZM1,M)
      VIZ2(NZ,M)=VIZ(NZM1,M)

  129 CONTINUE




      DO 53 L=2,NZ
C      DO 53 L=2,NZM1
      DO 53 M=1,NR
      HRFI(L,M)=HRFI2(L,M)
      HFI(L,M)=HRFI(L,M)/RAD(L,M)
      VIR(L,M)=VIR2(L,M)
      VIZ(L,M)=VIZ2(L,M)
      PLT(L,M)=PLT2(L,M)
      IF(BBK-0.) 205,205,206
  205      TEM(L,M)=PLT(L,M)**(GAM-1.)
      GOTO 53
  206      SS(L,M)=SS2(L,M)

C      TEM(L,M)=PLT(L,M)**GAM1*EXP(SS2(L,M))

      DAB(L,M)=PLT(L,M)**GAM*EXP(SS(L,M))

   53      CONTINUE


      DO 15 M=2,NRM1
C      DO 15 M=1,NR
      DO 16 L=1,NZ
      PTI(L)=RAD(L,M)*RDY(L)*PLT(L,M)
      UCK(L)=VIZ(L,M)
      D1(L)=0.
      D2(L)=0.
      D3(L)=0.
      C2(L)=0.
   16 P12DT(L)=PTI(L)
      PLB=0.
      PPL=PTI(1)
      UCK(NZ)=UCK(NZM1)
      PRB=1.
      PPR=0.

      PRB=RAD(NZ,M)*RDY(NZ)/RAD(NZM1,M)/RDY(NZM1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NZ,DZ,DT2,UCK,D1,D2,C2,D3)
      DO 17 L=1,NZ
   17      PLT12(L,M)=PDTI(L)/RAD(L,M)/RDY(L)
   15      CONTINUE


      DO 55 M=2,NRM1

      DO 56 L=1,NZ
      PTI(L)=RAD(L,M)*RDY(L)*PLT(L,M)*VIR(L,M)
      UCK(L)=VIZ(L,M)
      D1(L)=0.
      D2(L)=0.
      C2(L)=0.


      D3(L)=0.
   56 P12DT(L)=PTI(L)
      PLB=0.
      PPL=RDZ(1,M)*VIZ(1,M)*RAD(1,M)*RDY(1)*PLT(1,M)
      UCK(NZ)=UCK(NZM1)

      PPR=0.
      PRB=1.

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NZ,DZ,DT2,UCK,D1,D2,C2,D3)
      DO 57 L=1,NZ
   57      VIR12(L,M)=PDTI(L)/RAD(L,M)/RDY(L)/PLT12(L,M)
   55      CONTINUE


      DO 86 M=2,NRM1
c      DO 86 M=1,NR
      DO 87 L=1,NZ
      PTI(L)=RAD(L,M)*RDY(L)*PLT(L,M)*VIZ(L,M)
      UCK(L)=VIZ(L,M)

      D1(L)=0.

C      DABLM=BB/2.*TEM(L,M)*PLT(L,M)*(1.+ALF)
      DABLM=DAB(L,M)
      D2(L)=DABLM+HFI(L,M)**2/2.
      D3(L)=PVZ(L,M)
      D3(L)=0.
      C2(L)=-RAD(L,M)*RDY(L)

   87 P12DT(L)=PTI(L)

      UCK(NZ)=UCK(NZM1)



      PPL=0.
C      PLB=RAD(1,M)/RAD(2,M)
      PLB=1.
      PRB=1.
      PPR=0.

      PRB=RAD(NZ,M)*RDY(NZ)/RAD(NZM1,M)/RDY(NZM1)

C      PLB=0.
C      PPL=PTI(1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NZ,DZ,DT2,UCK,D1,D2,C2,D3)
      DO 88 L=1,NZ
   88      VIZ12(L,M)=PDTI(L)/RAD(L,M)/RDY(L)/PLT12(L,M)
   86      CONTINUE

C      M=1

C      LNEW=NZ-LKP2+1
C      LNEWM1=LNEW-1
C      DO 787 L=1,LNEW
C      LD=L+LKP2-1
C      PTI(L)=RAD(LD,M)*RDY(LD)*PLT(LD,M)*VIZ(LD,M)
C      UCK(L)=VIZ(LD,M)

C      D1(L)=0.

C      DABLM=BB/2.*TEM(LD,M)*PLT(LD,M)*(1.+ALF)
C      D2(L)=DABLM+HFI(LD,M)**2/2.
C      D2(L)=DABLM
C      D3(L)=0.
C      C2(L)=-RAD(LD,M)*RDY(LD)
CC      C2(L)=0.
C  787 P12DT(L)=PTI(L)

C      UCK(LNEW)=UCK(LNEWM1)

CC      UCK(1)=UCK(2)
C      UCK(1)=0.
C      PPL=0.
C      PLB=0.

C      PRB=1.
C      PPR=0.

C      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
C     *LNEW,DZ,DT2,UCK,D1,D2,C2,D3)

C      DO 788 L=LKP2,NZ
C  788      VIZ12(L,M)=PDTI(L-LKP2+1)/RAD(L,M)/RDY(L)/PLT(L,M)


      IF(BBK.EQ.0.) GOTO 207

C      DO 195 M=1,NR
      DO 195 M=2,NRM1
      DO 196 L=1,NZ
      PTI(L)=RAD(L,M)*RDY(L)*PLT(L,M)*SS(L,M)
      UCK(L)=VIZ(L,M)

      D1(L)=0.

      D2(L)=0.
      C2(L)=0.
C      D3(L)=0.
      D3(L)=EQ(L,M)/2.
C      D3(L)=EQ(L,M)

  196 P12DT(L)=PTI(L)

      PLB=0.
      PPL=PTI(1)

      UCK(NZ)=UCK(NZM1)

      PRB=1.

      PPR=0.

      PRB=RAD(NZ,M)*RDY(NZ)/RAD(NZM1,M)/RDY(NZM1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NZ,DZ,DT2,UCK,D1,D2,C2,D3)


      DO 197 L=1,NZ
  197      SS12(L,M)=PDTI(L)/RAD(L,M)/RDY(L)/PLT12(L,M)
  195      CONTINUE

  207 CONTINUE

      DO 25 M=2,NRM1
      DO 26 L=1,NZ
      PTI(L)=RDY(L)*HRFI(L,M)
      UCK(L)=VIZ(L,M)
      D1(L)=0.

      HF=HFI(L,M)

      D2(L)=0.
      C2(L)=0.

      PAR1=RDY(L)*VIR(L,M)*HF

      D3(L)=PAR1
C      D3(L)=0.

   26 P12DT(L)=PTI(L)
      PLB=0.
      PPL=PTI(1)
      UCK(NZ)=UCK(NZM1)

      PPR=0.

      PRB=RDY(NZ)/RDY(NZM1)*RAD(NZ,M)/RAD(NZM1,M)
C      PRB=1.
      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NZ,DZ,DT2,UCK,D1,D2,C2,D3)
      DO 27 L=1,NZ
   27      HRFI12(L,M)=PDTI(L)/RDY(L)
   25      CONTINUE

      DO 49 L=1,NZ
      HRFI12(L,1)=HRFI(L,1)
      HRFI12(L,NR)=HRFI(L,NR)
      SS12(L,1)=SS(L,1)
      SS12(L,NR)=SS(L,NR)
      PLT12(L,1)=PLT(L,1)
      VIZ12(L,1)=VIZ(L,1)
      VIR12(L,1)=VIR(L,1)
      PLT12(L,NR)=PLT(L,NR)
      VIZ12(L,NR)=VIZ(L,NR)
      VIR12(L,NR)=VIR(L,NR)

      VIR12(L,1)=DR2(L)*VIZ12(L,1)
         VIR12(L,NR)=DR1(L)*VIZ12(L,NR)

   49 CONTINUE

C      DO 749 L=1,LKP2M1
C      VIZ12(L,1)=VIZ(L,1)

C  749 CONTINUE

      DO 253 L=1,NZ
      DO 253 M=1,NR
      HFI12(L,M)=HRFI12(L,M)/RAD(L,M)
      IF(BBK-0.) 208,208,209
  208      TEM12(L,M)=PLT12(L,M)**(GAM-1.)
      GOTO 253
  209      CONTINUE

C      TEM12(L,M)=PLT12(L,M)**GAM1*EXP(SS12(L,M))
      DAB12(L,M)=PLT12(L,M)**GAM*EXP(SS12(L,M))


  253      CONTINUE



      DO 19 M=2,NRM1
C      DO 19 M=1,NR
      DO 20 L=1,NZ
      PTI(L)=RAD(L,M)*RDY(L)*PLT(L,M)
      UCK(L)=VIZ12(L,M)
      D1(L)=0.
      D2(L)=0.
      D3(L)=0.
      C2(L)=0.
   20 P12DT(L)=RAD(L,M)*RDY(L)*PLT12(L,M)
      PLB=0.
      PPL=PTI(1)
      UCK(NZ)=UCK(NZM1)

      PPR=0.
      PRB=1.

      PRB=RAD(NZ,M)*RDY(NZ)/RAD(NZM1,M)/RDY(NZM1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NZ,DZ,DT,UCK,D1,D2,C2,D3)
      DO 21 L=1,NZ
   21      PLT2(L,M)=PDTI(L)/RAD(L,M)/RDY(L)
   19      CONTINUE



      DO 59 M=2,NRM1

      DO 60 L=1,NZ
      PTI(L)=RAD(L,M)*RDY(L)*PLT(L,M)*VIR(L,M)
      UCK(L)=VIZ12(L,M)
      D1(L)=0.
      D2(L)=0.
      C2(L)=0.

      D3(L)=0.
   60 P12DT(L)=RAD(L,M)*RDY(L)*PLT12(L,M)*VIR12(L,M)
      PLB=0.
      PPL=RDZ(1,M)*VIZ12(1,M)*RAD(1,M)*RDY(1)*PLT12(1,M)

      UCK(NZ)=UCK(NZM1)


      PPR=0.
      PRB=1.
      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NZ,DZ,DT,UCK,D1,D2,C2,D3)
      DO 61 L=1,NZ
   61      VIR2(L,M)=PDTI(L)/RAD(L,M)/RDY(L)/PLT2(L,M)
   59      CONTINUE

      DO 89 M=2,NRM1
C      DO 89 M=1,NR
      DO 90 L=1,NZ
      PTI(L)=RAD(L,M)*RDY(L)*PLT(L,M)*VIZ(L,M)
      UCK(L)=VIZ12(L,M)
      D1(L)=0.
C      DABLM=BB/2.*TEM12(L,M)*PLT12(L,M)*(1.+ALF)
      DABLM=DAB12(L,M)
      D2(L)=DABLM+HFI12(L,M)**2/2.
      D3(L)=0.
      C2(L)=-RAD(L,M)*RDY(L)

   90 P12DT(L)=RAD(L,M)*RDY(L)*PLT12(L,M)*VIZ12(L,M)

      UCK(NZ)=UCK(NZM1)


      PPL=0.
C      PLB=RAD(1,M)/RAD(2,M)
      PLB=1.
      PRB=1.

      PPR=0.

      PRB=RAD(NZ,M)*RDY(NZ)/RAD(NZM1,M)/RDY(NZM1)

C      PLB=0.
C      PPL=PTI(1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NZ,DZ,DT,UCK,D1,D2,C2,D3)
      DO 91 L=1,NZ
   91      VIZ2(L,M)=PDTI(L)/RAD(L,M)/RDY(L)/PLT2(L,M)
   89      CONTINUE

C      M=1

C      LNEW=NZ-LKP2+1
C      LNEWM1=LNEW-1
C      DO 790 L=1,LNEW
C      LD=L+LKP2-1
C      PTI(L)=RAD(LD,M)*RDY(LD)*PLT(LD,M)*VIZ(LD,M)
C      UCK(L)=VIZ12(LD,M)

C      D1(L)=0.

C      DABLM=BB/2.*TEM12(LD,M)*PLT12(LD,M)*(1.+ALF)
C      D2(L)=DABLM+HFI12(LD,M)**2/2.
C      D2(L)=DABLM
C      D3(L)=0.
C      C2(L)=-RAD(LD,M)*RDY(LD)
Cc      C2(L)=0.
C  790 P12DT(L)=RAD(LD,M)*RDY(LD)*PLT12(LD,M)*VIZ12(LD,M)

C      UCK(LNEW)=UCK(LNEWM1)

CC      UCK(1)=UCK(2)
C      UCK(1)=0.
C      PPL=0.
C      PLB=0.

C      PRB=1.
C      PPR=0.

C      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
C     *LNEW,DZ,DT,UCK,D1,D2,C2,D3)

C      DO 791 L=LKP2,NZ
C  791      VIZ2(L,M)=PDTI(L-LKP2+1)/RAD(L,M)/RDY(L)/PLT(L,M)


      IF(BBK.EQ.0.) GOTO 210

C      DO 198 M=1,NR
      DO 198 M=2,NRM1
      DO 199 L=1,NZ
      PTI(L)=RAD(L,M)*RDY(L)*PLT(L,M)*SS(L,M)
      UCK(L)=VIZ12(L,M)

      D1(L)=0.

      D2(L)=0.
      C2(L)=0.
C      D3(L)=0.
      D3(L)=EQ(L,M)/2.
C      D3(L)=EQ(L,M)

  199 P12DT(L)=RAD(L,M)*RDY(L)*PLT12(L,M)*SS12(L,M)

      PLB=0.
      PPL=PTI(1)

      UCK(NZ)=UCK(NZM1)

      PRB=1.

      PPR=0.

      PRB=RAD(NZ,M)*RDY(NZ)/RAD(NZM1,M)/RDY(NZM1)

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NZ,DZ,DT,UCK,D1,D2,C2,D3)

      DO 200 L=1,NZ
  200      SS2(L,M)=PDTI(L)/RAD(L,M)/RDY(L)/PLT2(L,M)
  198      CONTINUE

  210 CONTINUE

      DO 29 M=2,NRM1
      DO 30 L=1,NZ
      PTI(L)=RDY(L)*HRFI(L,M)
      UCK(L)=VIZ12(L,M)
      D1(L)=0.
      D2(L)=0.
      C2(L)=0.

      HF=HFI12(L,M)

      PAR1=RDY(L)*VIR12(L,M)*HF

      D3(L)=PAR1
C      D3(L)=0.

   30 P12DT(L)=RDY(L)*HRFI12(L,M)
      PLB=0.
      PPL=PTI(1)
      UCK(NZ)=UCK(NZM1)

      PPR=0.
      PRB=RDY(NZ)/RDY(NZM1)*RAD(NZ,M)/RAD(NZM1,M)
C      PRB=1.

      CALL LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NZ,DZ,DT,UCK,D1,D2,C2,D3)
      DO 31 L=1,NZ
   31      HRFI2(L,M)=PDTI(L)/RDY(L)
   29      CONTINUE

      DO 949 L=1,NZ
      HRFI2(L,1)=HRFI(L,1)
      HRFI2(L,NR)=HRFI(L,NR)
      SS2(L,1)=SS(L,1)
      SS2(L,NR)=SS(L,NR)
      VIR2(L,1)=VIR(L,1)
      VIR2(L,NR)=VIR(L,NR)
      VIZ2(L,NR)=VIZ(L,NR)
      PLT2(L,NR)=PLT(L,NR)
      VIZ2(L,1)=VIZ(L,1)
      PLT2(L,1)=PLT(L,1)

      VIR2(L,1)=DR2(L)*VIZ2(L,1)
         VIR2(L,NR)=DR1(L)*VIZ2(L,NR)

  949 CONTINUE

C      DO 7949 L=1,LKP2M1
C      VIZ2(L,1)=VIZ(L,1)

C 7949 CONTINUE

C      DO 7949 L=LKP2,NZ
C      VIR2(L,1)=-VIR2(L,2)
C      VIR2(L,1)=0.
C      VIR2(L,1)=0.
C 7949 CONTINUE

      DO 62 L=1,NZ
      DO 62 M=1,NR
      HRFI(L,M)=HRFI2(L,M)
      HFI(L,M)=HRFI(L,M)/RAD(L,M)
      VIR(L,M)=VIR2(L,M)
      VIZ(L,M)=VIZ2(L,M)
      PLT(L,M)=PLT2(L,M)
      SS(L,M)=SS2(L,M)
      IF(BBK-0.) 211,211,212
  211      TEM(L,M)=PLT(L,M)**(GAM-1.)
      GOTO 62
  212      CONTINUE
C      TEM(L,M)=PLT2(L,M)**GAM1*EXP(SS2(L,M))
      DAB(L,M)=PLT(L,M)**GAM*EXP(SS(L,M))



   62      CONTINUE

      PLT(NZ,1)=PLT(NZM1,1)
          PLT(NZ,NR)=PLT(NZM1,NR)

      HRFI(NZ,1)=HRFI(NZM1,1)
     * *RDY(NZM1)/RDY(NZ)
C      */RAD(NZM1,1)*RAD(NZ,1)

          HRFI(NZ,NR)=HRFI(NZM1,NR)
     * *RDY(NZM1)/RDY(NZ)
C      */RAD(NZM1,NR)*RAD(NZ,NR)

      SS(NZ,1)=SS(NZM1,1)
          SS(NZ,NR)=SS(NZM1,NR)

C      TEM(NZ,1)=TEM(NZM1,1)
C          TEM(NZ,NR)=TEM(NZM1,NR)

      VIZ(NZ,1)=VIZ(NZM1,1)
          VIZ(NZ,NR)=VIZ(NZM1,NR)

      VIR(NZ,1)=DR2(NZ)*VIZ(NZ,1)
         VIR(NZ,NR)=DR1(NZ)*VIZ(NZ,NR)

C      VIZ(1,1)=VIZ(1,2)
      VIZ(1,1)=PLT(2,1)/PLT(1,1)*VIZ(2,1)
     * *RDY(2)/RDY(1)
     * *RAD(2,1)/RAD(1,1)

      VIZ(1,NR)=PLT(2,NR)/PLT(1,NR)*VIZ(2,NR)
     * *RDY(2)/RDY(1)
     * *RAD(2,NR)/RAD(1,NR)

      VIR(1,1)=RDZ(1,1)*VIZ(1,1)
      VIR(1,NR)=RDZ(1,NR)*VIZ(1,NR)

C      DO 249 M=1,NR
C      IF(VIZ(1,M).LT.0.) VIZ(1,M)=0.001
C  249 CONTINUE

C      GOTO 1000

C*********PROGONKA_FOR_HRFI

      IF(RE(1,1).EQ.0.) GOTO 1110

C      DO 223 L=LKP2,NZ
C      BAN1(L)=(RE(L,NR)+RE(L,NRM1))/2.*(HRFI(L,NRM1)-HRFI(L,NR))
C      BKA1(L)=(RE(L,1)+RE(L,2))/2.*(HRFI(L,1)-HRFI(L,2))
C  223 CONTINUE
      DO 224 M=1,NR
      BIX(M)=(RE(NZ,M)+RE(NZM1,M))/2.*(HRFI(NZM1,M)-HRFI(NZ,M))
  224 CONTINUE

      DTPR=DT

      DO 132 M=2,NRM1
      MKP=M
      CALL PRO_HZ(PROZ,HRFI,RE,CNUZ,CNURZ,BIX,
     *RDY,RDZ,DY,DZ,DTPR,MKP,RXAP,NZ,NR,DHDZL)

      DO 133 L=1,NZ
  133      HRFI12(L,M)=PROZ(L)

      DO 8149 L=1,NZ
 8149      DHDZLM(L,M)=DHDZL(L)

  132      CONTINUE

      DO 134 L=1,NZ
      HRFI12(L,1)=HRFI(L,1)
      HRFI12(L,NR)=HRFI(L,NR)
  134 CONTINUE

      DO 139 L=1,NZ
      DO 139 M=1,NR
      HRFI(L,M)=HRFI12(L,M)
      HFI(L,M)=HRFI(L,M)/RAD(L,M)
  139      CONTINUE




      DO 9150 L=2,NZ
      DHDZLM(L,1)=(HRFI(L,1)-HRFI(L-1,1))/DZ
      DHDZLM(L,NR)=(HRFI(L,NR)-HRFI(L-1,NR))/DZ
 9150 CONTINUE

      DHDZLM(1,1)=0.
      DHDZLM(1,NR)=0.

C      GOTO 1000

      DO 143 L=2,NZM1
      LKPL=L
      CALL PRO_HR(PROR,HRFI,RE,CNUR,CNUZ,CNURZ,
     *BAN1,BKA1,RAD,RDY,RDZ,DY,DZ,DTPR,LKPL,RXAP,NZ,NR,DHDRM)

      DO 144 M=1,NR
  144      HRFI12(L,M)=PROR(M)

      DO 5149 M=1,NR
 5149      DHDRLM(L,M)=DHDRM(M)

  143      CONTINUE

      DO 145 M=1,NR
      HRFI12(1,M)=HRFI(1,M)
      HRFI12(NZ,M)=HRFI(NZ,M)
  145 CONTINUE


      DO 7150 M=2,NR
      DHDRLM(1,M)=(HRFI12(1,M)-HRFI12(1,M-1))/DY
      DHDRLM(NZ,M)=(HRFI12(NZ,M)-HRFI12(NZ,M-1))/DY
 7150 CONTINUE

      DHDRLM(1,1)=0.
      DHDRLM(NZ,1)=0.

      DO 146 L=1,NZ
      DO 146 M=1,NR
      HRFI(L,M)=HRFI12(L,M)
      HFI(L,M)=HRFI(L,M)/RAD(L,M)
  146      CONTINUE
C      PRINT 712
C  712 FORMAT(1X,'WE are here 1')
C      PAUSE


 1110 CONTINUE


 1000 CONTINUE

C      EPS0=1.E-2
      DO 862 L=2,NZ
      DO 862 M=1,NR

      DABLM=DAB(L,M)
      PLOT=PLT(L,M)
C      PLOTKP=0.0001
C      IF(PLOT.LT.PLOTKP) PLOT=PLOTKP
      FKT=2.E0*DABLM/BB/PLOT
C      PRINT 7122,FKT
C 7122 FORMAT(2X,F15.6)
C      PAUSE

C      ALF=AL(L,M)
C      ALF2=AL(L,M)

C      NITER=20
C      DO 762 NIT=1,NITER
C      AINT=1.E-12
C      BINT=1.+AINT

C      TEM(L,M)=2.*DABLM/BB/PLOT/(1.+AL(L,M))

      TKPFKT=TKP/FKT
      IF(TKPFKT.GT.TKPIT) TKPFKT=TKPIT
C      IF(TKPFKT.GT.TKPIT) GO TO 5858


      CAB=CAX*BB/2.E0/DABLM
      TAB=CAB*FKT**2.5*EXP(-TKPFKT)
      TASB=TAB/(1.E0+TAB)
      ATAS=SQRT(TASB)
      ALKT=ATAS

C      TAS=CAX*FKT**1.5*EXP(-TKPFKT)
C      TASP=TAS/PLOT
C      QTAS=SQRT(TASP)
C      ALKT=-TASP/2.E0+QTAS*SQRT(TASP/4.E0+1.E0)
C 5858      CONTINUE

      AINT=FKT/2.E0-0.01
      IF(AINT.LT.1.E0) AINT=1.E0
C      BINT=1.5*TEM(L,M)

C      AINT=FKT/2.
      BINT=FKT*1.1

C      AINT=0.75
C      AINT=1.E0

C      BINT=170.E0
C      AINT=FKT

C      AINT=FKT/2.-0.01
C      BINT=FKT+0.01

      EPSI=1.E-5
      EPSI=0.01
C      EPSI=0.03

      PLOTKP=0.0001

C      IF(PLOT.LT.PLOTKP) RKOP=FKT/2.E0
C      IF(PLOT.LT.PLOTKP)  GO TO 762

           IF(ALKT.LT.0.000001) RKOP=FKT
      IF(ALKT.LT.0.000001) GO TO 762



      CALL RZEKOP(AINT,BINT,RKOP,EPSI)
C      ALF2=RKOP

C           IF(NITER.EQ.2) RKOP=2.*DABLM/BB/PLOT/(1.+AL(L,M))
  762    TEM(L,M)=RKOP


C      EPS=ABS(ALF-ALF2)/ALF
C      IF(EPS.LT.EPS0) GOTO 763

C      ALF=ALF2


C      ALF=1.

C      DABLM2=BB/2.*TEMP*PLOT*(1.+ALF2)

C      SS(L,M)=ALOG(DAB(L,M)/(PLT(L,M)**GAM))




C  762      CONTINUE

C      PRINT 4171
C 4171 FORMAT(1X,'NET TOCHNOCTI')



C  763      CONTINUE



C      AL(L,M)=ALF2
C      TEM(L,M)=2.*DABLM/BB/(1.+ALF2)/PLOT

C      IF(PLOT.LT.PLOTKP) AL(L,M)=1.
C      IF(PLOT.LT.PLOTKP)  GO TO 5959

      TKPTEM=TKP/TEM(L,M)
      IF(TKPTEM.GT.TKPIT) TKPTEM=TKPIT
C      IF(TKPTEM.GT.TKPIT) GO TO 5959

      CAB=CAX*BB/2.E0/DABLM
      TAB=CAB*TEM(L,M)**2.5*EXP(-TKPTEM)
      TASB=TAB/(1.E0+TAB)
      ATAS=SQRT(TASB)
      AL(L,M)=ATAS
 5959      CONTINUE



C      TAS=CAX*TEM(L,M)**1.5*EXP(-TKP/TEM(L,M))
C      TASP=TAS/PLOT
C      QTAS=SQRT(TASP)
C      AL(L,M)=-TASP/2.E0+QTAS*SQRT(TASP/4.E0+1.E0)
C      AL(L,M)=1.

C      DAB(L,M)=PLOT**GAM*EXP(SS(L,M))
C      TEM(L,M)=DAB(L,M)*2./BB/PLOT/(1.E0+AL(L,M))


C      DAB(L,M)=BB/2.*TEM(L,M)*PLT(L,M)*(1.+AL(L,M))
C      SS(L,M)=ALOG(DAB(L,M)/(PLT(L,M)**GAM))



  862      CONTINUE


      CALL QAPA(NZ,NR,TEM,HFI,PLT,AL,QCYM)

C      PRINT 100
C  100 FORMAT('          MACCIB QCYM')
C      DO 109 L=1,NZ
C      DO 109 M=1,NR
C      FF(L,M)=QCYM(L,M)
C  109      CONTINUE
C      NZK=NZ/10
C      DO 113 M=1,NR,8
C      MJ=NR-M+1
C      PRINT 119,(FF(L,MJ),L=1,NZ,NZK)
C  113 CONTINUE
C  119 FORMAT(11(1X,F11.10))
C      PAUSE

C      GOTO 1444
C*********PROGONKA_FOR_TEM
      DTPR=DT
C      GOTO 3
      DO 8224 M=1,NR
      BIX(M)=(QCYM(NZ,M)+QCYM(NZM1,M))/2.*(TEM(NZM1,M)-TEM(NZ,M))
 8224 CONTINUE

      DO 1321 M=2,NRM1
      MKP=M
      CALL PRO_TZ(PROZ,TEM,QCYM,PLT,BIX,
     *RDY,RDZ,RAD,DY,DZ,DTPR,MKP,NZ,NR)

      DO 1331 L=1,NZ
 1331      TEMP12(L,M)=PROZ(L)
 1321      CONTINUE

      DO 1341 L=1,NZ
      TEMP12(L,1)=TEMP12(L,2)
      TEMP12(L,NR)=TEMP12(L,NRM1)
 1341 CONTINUE

      DO 1391 L=1,NZ
      DO 1391 M=1,NR
      TEM(L,M)=TEMP12(L,M)
 1391      CONTINUE



      DO 1431 L=2,NZM1
      LKPL=L
      CALL PRO_TR(PROR,TEM,QCYM,PLT,
     *RAD,RDY,RDZ,DY,DZ,DTPR,LKPL,NZ,NR)

      DO 1441 M=1,NR
 1441      TEMP12(L,M)=PROR(M)
 1431      CONTINUE

      DO 1451 M=1,NR
      TEMP12(1,M)=TEM(1,M)
      TEMP12(NZ,M)=TEMP12(NZM1,M)
 1451 CONTINUE

      DO 2391 L=1,NZ
      DO 2391 M=1,NR
      TEM(L,M)=TEMP12(L,M)
 2391      CONTINUE

C      TEM(NZ,1)=TEM(NZM1,1)
C          TEM(NZ,NR)=TEM(NZM1,NR)

C      PRINT 7121
C 7121 FORMAT(1X,'TEM(NZ)=',15X,'TEM(NZM1)=',15X,15X,'M=')

C      DO 7124 M=2,NR,2
C      PRINT 7122,TEM(NZ,M),TEM(NZM1,M),M
C 7122 FORMAT(2(5X,F15.6),5X,I4)
C 7124 CONTINUE
C      PAUSE

 1444 CONTINUE

C      DO 8105 L=1,NZ
C      DO 8105 M=1,NR
C      TAS=CAX*TEM(L,M)**1.5*EXP(-TKP/TEM(L,M))
C      TASP=TAS/PLT(L,M)
C      QTAS=SQRT(TASP)
C      AL(L,M)=-TASP/2.E0+QTAS*SQRT(TASP/4.E0+1.E0)
C      DAB(L,M)=BB/2.*TEM(L,M)*PLT(L,M)*(1.E0+AL(L,M))
C      SS(L,M)=ALOG(DAB(L,M)/(PLT(L,M)**GAM))

 8105      CONTINUE

  333 CONTINUE

      TEMPNN(NN)=TEM(71,26)
      ALNN(NN)=AL(71,26)
      BPEM(NN)=TIME

      TEMPN2(NN)=TEM(61,26)
      ALN2(NN)=AL(61,26)
      BPEM2(NN)=TIME

      write(*,"(a,i4)") "    step - ",NN
    3 CONTINUE

      PRINT 121
  121 FORMAT(1X,'TEM(1)=',15X,'AL(1)=',15X,'M=')

      DO 124 M=1,NR,10
C       ZZ=(L-1.)*DZ
      PRINT 122,TEM(1,M),AL(1,M),M
  122 FORMAT(5X,F15.6,5X,E10.4,5X,I4)
  124 CONTINUE
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)



      PRINT 118
  118 FORMAT('          MACCIB PLT')
      DO 115 L=1,NZ
      DO 115 M=1,NR
      FF(L,M)=PLT(L,M)
  115      CONTINUE
      NZK=NZ/10
      DO 116 M=1,NR,4
      MJ=NR-M+1
      PRINT 117,(FF(L,MJ),L=1,NZ,NZK)
  116 CONTINUE
  117 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 318
  318 FORMAT('          MACCIB TEM')
      DO 315 L=1,NZ
      DO 315 M=1,NR
      FF(L,M)=TEM(L,M)
  315      CONTINUE
      NZK=NZ/10
      DO 316 M=1,NR,4
      MJ=NR-M+1
      PRINT 317,(FF(L,MJ),L=1,NZ,NZK)
  316 CONTINUE
  317 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 918
  918 FORMAT('          MACCIB AL')
      DO 915 L=1,NZ
      DO 915 M=1,NR
      FF(L,M)=AL(L,M)
  915      CONTINUE
      NZK=NZ/10
      DO 916 M=1,NR,4
      MJ=NR-M+1
      PRINT 917,(FF(L,MJ),L=1,NZ,NZK)
  916 CONTINUE
  917 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 2328
 2328 FORMAT('          MACCIB SS')
      DO 2325 L=1,NZ
      DO 2325 M=1,NR
      FF(L,M)=SS(L,M)
 2325      CONTINUE
      NZK=NZ/10
      DO 2326 M=1,NR,4
      MJ=NR-M+1
      PRINT 1327,(FF(L,MJ),L=1,NZ,NZK)
 2326 CONTINUE
 2327 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 718
  718 FORMAT('          MACCIB VIZ')
      DO 715 L=1,NZ
      DO 715 M=1,NR
      FF(L,M)=VIZ(L,M)
  715      CONTINUE
      NZK=NZ/10
      DO 716 M=1,NR,4
      MJ=NR-M+1
      PRINT 717,(FF(L,MJ),L=1,NZ,NZK)
  716 CONTINUE
  717 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 618
  618 FORMAT('          MACCIB VIR')
      DO 615 L=1,NZ
      DO 615 M=1,NR
      FF(L,M)=VIR(L,M)
  615      CONTINUE
      NZK=NZ/10
      DO 616 M=1,NR,4
      MJ=NR-M+1
      PRINT 617,(FF(L,MJ),L=1,NZ,NZK)
  616 CONTINUE
  617 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

      PRINT 3218
 3218 FORMAT('          MACCIB HRFI')
      DO 3215 L=1,NZ
      DO 3215 M=1,NR
      FF(L,M)=HRFI(L,M)
 3215      CONTINUE
      NZK=NZ/10
      DO 3216 M=1,NR,4
      MJ=NR-M+1
      PRINT 3217,(FF(L,MJ),L=1,NZ,NZK)
 3216 CONTINUE
 3217 FORMAT(11(1X,F6.2))
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)


      SUM1=0.
      SUM2=0.
       DO 2712 L=1,NZM1
       DO 2712 M=1,NR
      SUM1=SUM1+(FOTEM(L,M)-TEM(L,M))**2
      SUM2=SUM2+FOTEM(L,M)**2
 2712      CONTINUE
      DTEM=SQRT(SUM1/SUM2)*100.



      SUM1=0.
      SUM2=0.
       DO 2415 L=1,NZM1
       DO 2415 M=1,NR
      SUM1=SUM1+(FOVZ(L,M)-VIZ(L,M))**2
      SUM2=SUM2+FOVZ(L,M)**2
 2415      CONTINUE
      DVZ=SQRT(SUM1/SUM2)*100.

      NZM2=NZ-2
      NZM3=NZ-3

      SUM1=0.
      SUM2=0.
       DO 2416 L=1,NZM2
       DO 2416 M=1,NR
      SUM1=SUM1+(FOVR(L,M)-VIR(L,M))**2
      SUM2=SUM2+FOVR(L,M)**2
 2416      CONTINUE
      DVR=SQRT(SUM1/SUM2)*100.


      SUM1=0.
      SUM2=0.
       DO 2418 L=1,NZM1
       DO 2418 M=1,NR
      SUM1=SUM1+(FOPL(L,M)-PLT(L,M))**2
      SUM2=SUM2+FOPL(L,M)**2
 2418      CONTINUE
      PLOT=SQRT(SUM1/SUM2)*100.

      SUM1=0.
      SUM2=0.
       DO 2419 L=1,NZM1
       DO 2419 M=1,NR
      SUM1=SUM1+(FORH(L,M)-HRFI(L,M))**2
      SUM2=SUM2+FORH(L,M)**2
 2419      CONTINUE
      DHR=SQRT(SUM1/SUM2)*100.


      PRINT 2121
 2121 FORMAT(1X,' PLT=',12X,' DVZ=',12X,'DVR=',12X,'DHR',12X,'DTEM')

      PRINT 2122,PLOT,DVZ,DVR,DHR,DTEM
 2122 FORMAT(5(2X,F11.4))



      PRINT 2125,TIME
 2125 FORMAT(1X,' TIME=',F12.3)
      !PAUSE                    ! kv (commented)
      write(*,"(a)")
      write(*,"(a)") " pause: press Enter to continue ..."
      !read(*,*)

C      OPEN(1,FILE='calul.dat')
      WRITE(1,*) PCI,HRFI,FHZ,FHR,PLT,TEM,VIF,VIR,VIZ,SS,AL

C      OPEN(2,FILE='figs_p_T.dat')
C      DO 251 M=1,NR
C      YY=(M-1.)*DY
C      WRITE(2,*) YY,PLT(141,M),TEM(141,M)
C  251      CONTINUE

C      OPEN(4,FILE='figs_Vz_Vr_Hfi.dat')
C      DO 5257 M=1,NR
C      YY=(M-1.)*DY
C      WRITE(4,*) VIZ(141,M),VIR(141,M),HFI(141,M)
C 5257      CONTINUE

C      L=121

C      OPEN(12,FILE='figs121_p_T.dat')
C      DO 9256 M=1,NR
C      YY=(M-1.)*DY
C      WRITE(12,*) YY,PLT(L,M),TEM(L,M)
C 9256      CONTINUE

C      OPEN(14,FILE='figs121_Vz_Vr_Hfi.dat')
C      DO 5251 M=1,NR
C      YY=(M-1.)*DY
C      WRITE(14,*) VIZ(L,M),VIR(L,M),HFI(L,M)
C 5251      CONTINUE

      M=26
      OPEN(9,FILE=trim(TDPDir)//'ZM26_p_T_rH.dat')
      DO 951 L=1,NZ
      ZZZ=(L-1.)*DZ
      WRITE(9,*) ZZZ,PLT(L,M),TEM(L,M),HRFI(L,M)
  951      CONTINUE

      OPEN(10,FILE=trim(TDPDir)//'ZM26_Vz_Vr_Al_Vmod.dat')
      OPEN(19,FILE=trim(TDPDir)//'ZM26_VGAZ_VZBYK.dat')
      DO 9251 L=1,NZ
      ZZZ=(L-1.)*DZ
      VMOD=SQRT(VIZ(L,M)**2+VIR(L,M)**2)
      WRITE(10,*) VIZ(L,M),VIR(L,M),AL(L,M),VMOD
      VGAZ=SQRT(BB/2.*(1.+AL(L,M))*TEM(L,M)*GAM)
      VZBYK=SQRT(VGAZ**2+HFI(L,M)**2/PLT(L,M))
      WRITE(19,*) VGAZ,VZBYK

 9251      CONTINUE

      OPEN(21,FILE=trim(TDPDir)//'TIME_L71_M26_T_Al.dat')
      DO 8951 NNN=1,NDT
      WRITE(21,*) BPEM(NNN),TEMPNN(NNN),ALNN(NNN)
 8951      CONTINUE

      OPEN(22,FILE=trim(TDPDir)//'TIME_L61_M26_T_Al.dat')
      DO 8351 NNN=1,NDT
      WRITE(22,*) BPEM2(NNN),TEMPN2(NNN),ALN2(NNN)
 8351      CONTINUE

C      M=11
C      OPEN(19,FILE='ZM_p_T_rH.dat')
C      DO 2952 L=1,NZ
C      ZZZ=(L-1.)*DZ
C      WRITE(19,*) ZZZ,PLT(L,M),TEM(L,M),HRFI(L,M)
C 2952      CONTINUE

C      OPEN(18,FILE='ZM_Vz_Vr.dat')
C      DO 4254 L=1,NZ
C      ZZZ=(L-1.)*DZ
C      WRITE(18,*) VIZ(L,M),VIR(L,M)
C 4254      CONTINUE

C      L=281
C      OPEN(12,FILE='figs_p_TeV.dat')
C      DO 3251 M=1,NR
C
C      YY=(M-1.)*DY

C      TEB=TEM(L,M)*TB/1.16/10**4

C      WRITE(12,*) YY,PLT(L,M),TEB
C 3251      CONTINUE

C      OPEN(11,FILE='figs_VzV0_VrV0_Hfi.dat')
C      DO 5221 M=1,NR
C      YY=(M-1.)*DY
C      VzV0=VIZ(L,M)*VALF0
C      VrV0=VIR(L,M)*VALF0

C      WRITE(11,*) VzV0,VrV0,HFI(L,M)
C 5221      CONTINUE

C      M=2
C      OPEN(13,FILE='ZM2_p_TeV_rH.dat')
C      DO 9951 L=1,NZ
C      ZZZ=(L-1.)*DZ
C      TEB=TEM(L,M)*TB/1.16/10**4

C      WRITE(13,*) ZZZ,PLT(L,M),TEB,HRFI(L,M)
C 9951      CONTINUE

C      OPEN(14,FILE='ZM2_VzV0_VrV0.dat')
C      DO 9551 L=1,NZ
C      ZZZ=(L-1.)*DZ

C      VzV0=VIZ(L,M)*VALF0
C      VrV0=VIR(L,M)*VALF0

C      WRITE(14,*) VzV0,VrV0
C 9551      CONTINUE


C      M=41
C      OPEN(15,FILE='ZM41_p_TeV_rH.dat')
C      DO 2951 L=1,NZ
C      ZZZ=(L-1.)*DZ
C      TEB=TEM(L,M)*TB/1.16/10**4

C      WRITE(15,*) ZZZ,PLT(L,M),TEB,HRFI(L,M)
C 2951      CONTINUE

C      OPEN(16,FILE='ZM41_VzV0_VrV0_VmV0.dat')
C      DO 2551 L=1,NZ
C      ZZZ=(L-1.)*DZ

C      VzV0=VIZ(L,M)*VALF0
C      VrV0=VIR(L,M)*VALF0
C      VmV0=SQRT(VzV0**2+VrV0**2)
C      WRITE(16,*) VzV0,VrV0,VmV0
C 2551      CONTINUE


      M=2
C      M=161

C      OPEN(5,FILE='M161_Z_WETE_QELN5_QPOL5.dat')
C      OPEN(6,FILE='M161_Z_QLIN_QPEK_QTOP.dat')
C      OPEN(8,FILE='ZM161_QFOT_RETOK2.dat')

      OPEN(5,FILE=trim(TDPDir)//'M2_Z_WETE_QELN5_QPOL5.dat')
      OPEN(6,FILE=trim(TDPDir)//'M2_Z_QLIN_QPEK_QTOP.dat')
      OPEN(8,FILE=trim(TDPDir)//'ZM2_QFOT_RETOK2.dat')
       DO 851 L=1,NZ
      ZZZ=(L-1.)*DZ

      TBD1000=TB/1000.

      TEMPR=TEM(L,M)
      PLOT=PLT(L,M)
      POLE=ABS(HFI(L,M))

      FWETE=0.28/1000.*TBD1000**2/SQRT(BB*PB16)
      WETE=FWETE*POLE*SQRT(TEMPR**3)/PLOT
      GWETE=(11.92+4.66*WETE**2)/(3.77+14.79*WETE**2+WETE**4)
      FEL1=1.62/100000.*SQRT(PM*BB**3)*TBD1000**2/PB16/DLCM/ZION
      QELN=FEL1*SQRT(TEMPR**5)*GWETE
      QELN5=QELN*10**5
      GWETE2=WETE*(21.67+2.5*WETE**2)/(3.77+14.79*WETE**2+WETE**4)
      QPOL=FEL1*SQRT(TEMPR**5)*GWETE2
      QPOL5=QPOL*10**5

      FGB=(GAM-1.)/BB*RDY(L)/TEM(L,M)

      TEB=TEM(L,M)*TB/1.16/10**4
      FQ=QB2*PB16**2*PLT(L,M)**2/10**5
      QLIN=6.25*FQ*ZION**4
      QLIN=80.*FQ*ZION**7/SQRT(TEB**3)

      QPEK=4.4*FQ*ZION**5/SQRT(TEB)
      QTOP=0.154*FQ*ZION**3*SQRT(TEB)
      QFOT=QLIN+QPEK+QTOP
C      QFOT=QPEK+QTOP

      QCB1=QFOT*FGB*RAD(L,M)
      QCB0=EQ(L,M)+QCB1
      RETOK2=QCB0/FGB/RAD(L,M)


      WRITE(5,*) ZZZ,WETE,QELN5,QPOL5
      WRITE(6,*) ZZZ,QLIN,QPEK,QTOP
      WRITE(8,*) QFOT,RETOK2
  851 CONTINUE


      DO 7251 M=1,NR
      PEX(M)=PLT(66,M)
      TEX(M)=TEM(66,M)
      VZEX(M)=VIZ(66,M)
      VREX(M)=VIR(66,M)
      HREX(M)=HRFI(66,M)

 7251      CONTINUE

      OPEN(7,FILE=trim(TDPDir)//'R_PTVzVrHr.dat')
      WRITE(7,*) PEX,TEX,VZEX,VREX,HREX

      L=273
      L=69

      SUM4=0.
      DO 8251 M=1,NRM1
C       L=NZ
      PVNEX=PLT(L,M)*VIZ(L,M)*RAD(L,M)*(RAD(L,M+1)-RAD(L,M))
      SUM4=SUM4+PVNEX

 8251      CONTINUE

C      L=NZ
           SUM5=0.
      DO 8451 M=1,NRM1
      PVNEX=PLT(L,M)*VIZ(L,M)**2*RAD(L,M)*(RAD(L,M+1)-RAD(L,M))
      SUM5=SUM5+PVNEX

 8451      CONTINUE

C      L=NZ
           SUM6=0.
      DO 8651 M=1,NRM1
      VVD2=(VIZ(L,M)**2+VIR(L,M)**2+VIF(L,M)**2)/2.
      PVNEX=PLT(L,M)*VIZ(L,M)*VVD2*RAD(L,M)*(RAD(L,M+1)-RAD(L,M))
      SUM6=SUM6+PVNEX

 8651      CONTINUE

C      FMAN=2.*PI*DZ*RXAP*SUM1
C      FMKA=2.*PI*SUM3
C      FMIN=2.*PI*SUM2
      FMEX=2.*PI*SUM4
      FPEX=2.*PI*SUM5
      FNEX=2.*PI*SUM6
C      DFM=FMIN/ABS(FMKA)
C      DJM=2.*PI*E1*RXAP/ABS(FMEX)
C      PRINT 7126,FMAN,FMKA,FMIN,FMEX,DJM
C      PRINT 7126,FMAN,FMKA,FMIN,FMEX
C 7126 FORMAT(' FMAN=',F8.3,1X,'FMKA=',F8.3,1X,'FMIN=',F8.3,1X,
C     *'FMEX=',F8.3)
C     *'FMEX=',F8.3,1X,'DJM=',F8.3)

      FLNJP=ALOG(TOKBKA)
C      FLNPAC=ALOG(TOKBKA/DJM)
      FMEXGC=1.67/10**8*PM*PB16*VALF*DLCM**2*FMEX
      FLNPAC=ALOG(FMEXGC*1.6*100./PM/1.67)

      PEXCI=1.67/10**8*PM*PB16*VALF**2*DLCM**2*FPEX/10**5
      FNEXCI=1.67/10**8*PM*PB16*VALF**3*DLCM**2*FNEX/10**7/10**6
      PRINT 3123,TOKBKA,FMEXGC,FLNJP,FLNPAC
 3123 FORMAT(10X,' TOKBKA=',F8.3,1X,'FMEXGC=',F8.3,1X,
     *'FLNJP=',F8.3,1X,'FLNPAC=',F8.3)
C     *'FMEXGC=',F8.3,20X,'PEXCI(H)=',F10.3,1X,'FNEXCI(MBT)=',F11.3)
      PRINT 7127,TOKBKA,FMEXGC,PEXCI,FNEXCI
 7127 FORMAT(5X,' TOKBKA=',F8.3,10X,
     *'FMEXGC=',F8.3,10X,'PEXCI(H)=',F10.3,10X,'FNEXCI(MBT)=',F11.3)

      DEALLOCATE(TEMPNN,ALNN,BPEM)
      DEALLOCATE(TEMPN2,ALN2,BPEM2)

      STOP
      END


      SUBROUTINE LCP(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NH,DX,DT,VFI,D1,D2,C2,D3)
      DIMENSION PTI(NH),PDTI(NH),P12DT(NH),
     *VFI(NH),D1(NH),D2(NH),C2(NH),D3(NH),D(NH)
      DIMENSION FMU(NH),FTMD(NH),FDNU(NH),FTNU(NH),FDNUT(NH)
      DTDX=DT/DX
      NHM1=NH-1
      DO 1 N=1,NHM1
      VF12=(VFI(N+1)+VFI(N))/2.
      P12=(P12DT(N+1)+P12DT(N))/2.
C      E12=DTDX*VF12
C*******2-oi
      E121=DTDX*VF12
      E12=E121-DTDX*E121/2.*(VFI(N+1)-VFI(N))
C*******2-oi
      FKON=E12*P12
      FNU=(0.5+E12**2)/3.
      FMU(N)=(1.-E12**2)/6.
C      FD1=-DTDX/2.*(D1(N+1)+D1(N))
C*******2-oi
      A12=(1.+E121)/2.
      B12=(1.-E121)/2.
      FD1=-DTDX*(A12*D1(N)+B12*D1(N+1))
C*******2-oi
      FTMD(N)=FKON+FD1
      FDNU(N)=-FNU*(PTI(N+1)-PTI(N))
C*******2-oi
C      FNUT=E12**2/2.
C      FDNUT(N)=-FNUT*(PTI(N+1)-PTI(N))
C*******2-oi
    1 CONTINUE
      DO 2 N=2,NHM1
C      QI=DT*D3(N)
C      DFD2=DTDX/2.*C2(N)*(D2(N+1)-D2(N-1))
C*******2-oi
      QI=DT*(D3(N)-DTDX/4.*(D3(N)*(VFI(N+1)-VFI(N-1))+
     *VFI(N)*(D3(N+1)-D3(N-1))))
      VCP12=(C2(N+1)+C2(N))*(VFI(N+1)+VFI(N))/4.
      VCM12=(C2(N-1)+C2(N))*(VFI(N-1)+VFI(N))/4.
      DFD2=DTDX/2.*(C2(N)*(D2(N+1)-D2(N-1))-
     *DTDX*(VCP12*(D2(N+1)-D2(N))-VCM12*(D2(N)-D2(N-1))))
C*******2-oi
      PTMD=PTI(N)-FTMD(N)+FTMD(N-1)+DFD2+QI
      D(N)=PTMD-FDNU(N)+FDNU(N-1)
C      PDTI(N)=PTMD
C*******2-oi
C      PDTI(N)=PTMD-FDNUT(N)+FDNUT(N-1)
C*******2-oi
    2 CONTINUE
C      PDTI(1)=PLB*PDTI(2)+PPL
C      PDTI(NH)=PRB*PDTI(NHM1)+PPR
      D(1)=PLB*D(2)+PPL
      D(NH)=PRB*D(NHM1)+PPR
C      PDTI(1)=D1(1)
C      PDTI(NH)=D1(NH)

      DO 3 N=1,NHM1
C      FTNU(N)=FMU(N)*(PDTI(N+1)-PDTI(N))
      FTNU(N)=FMU(N)*(D(N+1)-D(N))

    3 CONTINUE

      NHM2=NH-2
      F1=0.
      DO 4 N=2,NHM2
      S=SGN(D(N+1)-D(N))
      F2=ABS(FTNU(N))
      F3=S*(D(N+2)-D(N+1))
      F4=S*(D(N)-D(N-1))
      F21=AMIN1(F2,F3,F4)
      FDNU(N)=S*AMAX1(F1,F21)
    4 CONTINUE
      S=SGN(D(2)-D(1))
      F2=ABS(FTNU(1))
      F3=S*(D(3)-D(2))
      F21=AMIN1(F2,F3)
      FDNU(1)=S*AMAX1(F1,F21)
      S=SGN(D(NH)-D(NHM1))
      F2=ABS(FTNU(NHM1))
      F4=S*(D(NHM1)-D(NHM2))
      F21=AMIN1(F2,F4)
      FDNU(NHM1)=S*AMAX1(F1,F21)
      DO 5 N=2,NHM1
      PDTI(N)=D(N)-FDNU(N)+FDNU(N-1)
    5 CONTINUE
      PDTI(1)=PLB*PDTI(2)+PPL
      PDTI(NH)=PRB*PDTI(NHM1)+PPR
      RETURN
      END

      SUBROUTINE LCP2(PTI,PDTI,P12DT,PLB,PPL,PRB,PPR,
     *NH,DX,DT,VFI,D1,D2,C2,D3,D22,C22)
      DIMENSION PTI(NH),PDTI(NH),P12DT(NH),D(NH),
     *VFI(NH),D1(NH),D2(NH),C2(NH),D3(NH),D22(NH),C22(NH)
      DIMENSION FMU(NH),FTMD(NH),FDNU(NH),FTNU(NH),FDNUT(NH)
      DTDX=DT/DX
      NHM1=NH-1
      DO 1 N=1,NHM1
      VF12=(VFI(N+1)+VFI(N))/2.
      P12=(P12DT(N+1)+P12DT(N))/2.
C      E12=DTDX*VF12
C*******2-oi
      E121=DTDX*VF12
      E12=E121-DTDX*E121/2.*(VFI(N+1)-VFI(N))
C*******2-oi
      FKON=E12*P12
      FNU=(0.5+E12**2)/3.
      FMU(N)=(1.-E12**2)/6.
C      FD1=-DTDX/2.*(D1(N+1)+D1(N))
C*******2-oi
      A12=(1.+E121)/2.
      B12=(1.-E121)/2.
      FD1=-DTDX*(A12*D1(N)+B12*D1(N+1))
C*******2-oi
      FTMD(N)=FKON+FD1
      FDNU(N)=-FNU*(PTI(N+1)-PTI(N))
C*******2-oi
C      FNUT=E12**2/2.
C      FDNUT(N)=-FNUT*(PTI(N+1)-PTI(N))
C*******2-oi
    1 CONTINUE
      DO 2 N=2,NHM1
C      QI=DT*D3(N)
C      DFD2=DTDX/2.*C2(N)*(D2(N+1)-D2(N-1))
C      DFD22=DTDX/2.*C22(N)*(D22(N+1)-D22(N-1))

C*******2-oi
      QI=DT*(D3(N)-DTDX/4.*(D3(N)*(VFI(N+1)-VFI(N-1))+
     *VFI(N)*(D3(N+1)-D3(N-1))))
      VCP12=(C2(N+1)+C2(N))*(VFI(N+1)+VFI(N))/4.
      VCM12=(C2(N-1)+C2(N))*(VFI(N-1)+VFI(N))/4.
      DFD2=DTDX/2.*(C2(N)*(D2(N+1)-D2(N-1))-
     *DTDX*(VCP12*(D2(N+1)-D2(N))-VCM12*(D2(N)-D2(N-1))))
      VCP22=(C22(N+1)+C22(N))*(VFI(N+1)+VFI(N))/4.
      VCM22=(C22(N-1)+C22(N))*(VFI(N-1)+VFI(N))/4.
      DFD22=DTDX/2.*(C22(N)*(D22(N+1)-D22(N-1))-
     *DTDX*(VCP22*(D22(N+1)-D22(N))-VCM22*(D22(N)-D22(N-1))))

C*******2-oi
      PTMD=PTI(N)-FTMD(N)+FTMD(N-1)+DFD2+DFD22+QI
      D(N)=PTMD-FDNU(N)+FDNU(N-1)
C      PDTI(N)=PTMD
C*******2-oi
C      PDTI(N)=PTMD-FDNUT(N)+FDNUT(N-1)
C*******2-oi
    2 CONTINUE
C      PDTI(1)=PLB*PDTI(2)+PPL
C      PDTI(NH)=PRB*PDTI(NHM1)+PPR
C      D(1)=PLB*D(2)+PPL
C      D(NH)=PRB*D(NHM1)+PPR
C      PDTI(1)=D(1)
C      PDTI(NH)=D(NH)
C      D(1)=PDTI(1)
C      D(NH)=PDTI(NH)

      D(1)=PLB*D(2)+PPL
      D(NH)=PRB*D(NHM1)+PPR

      DO 3 N=1,NHM1
C      FTNU(N)=FMU(N)*(PDTI(N+1)-PDTI(N))
      FTNU(N)=FMU(N)*(D(N+1)-D(N))

    3 CONTINUE

      NHM2=NH-2
      F1=0.
      DO 4 N=2,NHM2
      S=SGN(D(N+1)-D(N))
      F2=ABS(FTNU(N))
      F3=S*(D(N+2)-D(N+1))
      F4=S*(D(N)-D(N-1))
      F21=AMIN1(F2,F3,F4)
      FDNU(N)=S*AMAX1(F1,F21)
    4 CONTINUE
      S=SGN(D(2)-D(1))
      F2=ABS(FTNU(1))
      F3=S*(D(3)-D(2))
      F21=AMIN1(F2,F3)
      FDNU(1)=S*AMAX1(F1,F21)
      S=SGN(D(NH)-D(NHM1))
      F2=ABS(FTNU(NHM1))
      F4=S*(D(NHM1)-D(NHM2))
      F21=AMIN1(F2,F4)
      FDNU(NHM1)=S*AMAX1(F1,F21)
      DO 5 N=2,NHM1
      PDTI(N)=D(N)-FDNU(N)+FDNU(N-1)
    5 CONTINUE
      PDTI(1)=PLB*PDTI(2)+PPL
      PDTI(NH)=PRB*PDTI(NHM1)+PPR
      RETURN
      END

      SUBROUTINE LCPR(PTI,PDTI,P12DT,PV1,PVN,
     *NH,DX,DT,VFI,D1,D2,C2,D3)
      DIMENSION PTI(NH),PDTI(NH),P12DT(NH),
     *VFI(NH),D1(NH),D2(NH),C2(NH),D3(NH),D(NH)
      DIMENSION FMU(NH),FTMD(NH),FDNU(NH),FTNU(NH),FDNUT(NH)

      DTDX=DT/DX
      NHM1=NH-1
      DO 1 N=1,NHM1
      VF12=(VFI(N+1)+VFI(N))/2.
      P12=(P12DT(N+1)+P12DT(N))/2.
C      E12=DTDX*VF12
C*******2-oi
      E121=DTDX*VF12
      E12=E121-DTDX*E121/2.*(VFI(N+1)-VFI(N))
C*******2-oi
      FKON=E12*P12
      FNU=(0.5+E12**2)/3.
      FMU(N)=(1.-E12**2)/6.
C      FD1=-DTDX/2.*(D1(N+1)+D1(N))
C*******2-oi
      A12=(1.+E121)/2.
      B12=(1.-E121)/2.
      FD1=-DTDX*(A12*D1(N)+B12*D1(N+1))
C*******2-oi
      FTMD(N)=FKON+FD1
      FDNU(N)=-FNU*(PTI(N+1)-PTI(N))
C*******2-oi
      FNUT=E12**2/2.
      FDNUT(N)=-FNUT*(PTI(N+1)-PTI(N))
C*******2-oi
    1 CONTINUE
      DO 2 N=2,NHM1
C      QI=DT*D3(N)
C      DFD2=DTDX/2.*C2(N)*(D2(N+1)-D2(N-1))
C*******2-oi
      QI=DT*(D3(N)-DTDX/4.*(D3(N)*(VFI(N+1)-VFI(N-1))+
     *VFI(N)*(D3(N+1)-D3(N-1))))
      VCP12=(C2(N+1)+C2(N))*(VFI(N+1)+VFI(N))/4.
      VCM12=(C2(N-1)+C2(N))*(VFI(N-1)+VFI(N))/4.
      DFD2=DTDX/2.*(C2(N)*(D2(N+1)-D2(N-1))-
     *DTDX*(VCP12*(D2(N+1)-D2(N))-VCM12*(D2(N)-D2(N-1))))
C*******2-oi
      PTMD=PTI(N)-FTMD(N)+FTMD(N-1)+DFD2+QI
      D1(N)=PTMD-FDNU(N)+FDNU(N-1)
C      PDTI(N)=PTMD
C*******2-oi
      PDTI(N)=PTMD-FDNUT(N)+FDNUT(N-1)
C*******2-oi
    2 CONTINUE
      F1=DTDX*(PV1-D1(1))
      FN=DTDX*(PVN-D1(NH))
      FNM12=FTMD(NHM1)+FDNU(NHM1)
      F32=FTMD(1)+FDNU(1)
      D1(1)=PTI(1)-2.*(F32-F1)+DT*D3(1)
     *+DTDX*C2(1)*(D2(2)-D2(1))
      D1(NH)=PTI(NH)-2.*(FN-FNM12)+DT*D3(NH)
     *+DTDX*C2(NH)*(D2(NH)-D2(NHM1))
      PDTI(1)=D1(1)
      PDTI(NH)=D1(NH)
      DO 3 N=1,NHM1
      FTNU(N)=FMU(N)*(PDTI(N+1)-PDTI(N))
    3 CONTINUE
      NHM2=NH-2
      F1=0.
      DO 4 N=2,NHM2
      S=SGN(D1(N+1)-D1(N))
      F2=ABS(FTNU(N))
      F3=S*(D1(N+2)-D1(N+1))
      F4=S*(D1(N)-D1(N-1))
      F21=AMIN1(F2,F3,F4)
      FDNU(N)=S*AMAX1(F1,F21)
    4 CONTINUE
      S=SGN(D1(2)-D1(1))
      F2=ABS(FTNU(1))
      F3=S*(D1(3)-D1(2))
      F21=AMIN1(F2,F3)
      FDNU(1)=S*AMAX1(F1,F21)
      S=SGN(D1(NH)-D1(NHM1))
      F2=ABS(FTNU(NHM1))
      F4=S*(D1(NHM1)-D1(NHM2))
      F21=AMIN1(F2,F4)
      FDNU(NHM1)=S*AMAX1(F1,F21)
      DO 5 N=2,NHM1
      PDTI(N)=D1(N)-FDNU(N)+FDNU(N-1)
    5 CONTINUE
C      PDTI(1)=PLB*PDTI(2)+PPL
C      PDTI(NH)=PRB*PDTI(NHM1)+PPR
      RETURN
      END



      FUNCTION SGN(X)
      IF(X) 5,8,13
   5  SGN=-1.
      GO TO 15
   8  SGN=0.
      GO TO 15
  13  SGN=1.
  15  CONTINUE
      RETURN
      END

      SUBROUTINE PRO_HZ(PROZ,HRFI,RE,CNUZ,CNURZ,BIX,
     *RDY,RDZ,DY,DZ,DTPR,M,RXAP,NZ,NR,DHDZL)
C      PARAMETER (NZ=161,NR=81)
      DIMENSION RE(NZ,NR),CNUZ(NZ,NR),CNURZ(NZ,NR)
      DIMENSION HRFI(NZ,NR),PROZ(NZ),RDY(NZ),RDZ(NZ,NR)
      DIMENSION BIX(NR),DHDZL(NZ)
      DIMENSION A(NZ),B(NZ),C(NZ),F(NZ),BET(NZ),
     *GAMM(NZ),DELT(NZ),WL(NZ)
      NZM1=NZ-1
      DO 1 L=2,NZ
      C(L)=RDY(L)/DTPR
      A(L)=(RDY(L)*CNUZ(L,M)+RDY(L-1)*CNUZ(L-1,M))/2./DZ**2

    1 CONTINUE

      DO 2 L=2,NZM1
      B(L)=(RDY(L)*CNUZ(L,M)+RDY(L+1)*CNUZ(L+1,M))/2./DZ**2
      FLZYD2=-(RDZ(L,M)*CNUZ(L,M)+CNURZ(L,M))*RE(L,M)/4./DZ/DY*
     *(HRFI(L+1,M+1)-HRFI(L+1,M-1)-HRFI(L-1,M+1)+HRFI(L-1,M-1))
      FL1=-(RE(L,M+1)*RDZ(L,M+1)*CNUZ(L,M+1)-
     *RE(L,M-1)*RDZ(L,M-1)*CNUZ(L,M-1))/4./DZ/DY*
     *(HRFI(L+1,M)-HRFI(L-1,M))
      F(L)=RDY(L)*HRFI(L,M)/DTPR+FLZYD2+FL1

    2 CONTINUE

      GAP1=1.
C      DAMD1=0.
      DAMD1=2./(RE(1,M)+RE(2,M))
      XUP1=-RXAP

      CNAM=A(2)*GAP1+C(2)*DAMD1
      BET(2)=A(2)*DAMD1/CNAM
      GAMM(2)=A(2)*XUP1/CNAM

      DO 3 L=2,NZM1
      RE12=(RE(L+1,M)+RE(L,M))/2.
C      CIGM=1./RE12
      DELT(L)=A(L)*(A(L+1)*RE12+C(L+1))+RE12*B(L)*C(L+1)*BET(L)
C      DELT(L)=A(L)*(A(L+1)+C(L+1)*CIGM)+B(L)*C(L+1)*BET(L)

      BET(L+1)=A(L+1)/DELT(L)*(RE12*B(L)*BET(L)+A(L))
C      BET(L+1)=A(L+1)/DELT(L)*(B(L)*BET(L)+A(L)*CIGM)
      GAMM(L+1)=A(L+1)*RE12/DELT(L)*(A(L)*GAMM(L)+BET(L)*F(L))
C      GAMM(L+1)=A(L+1)/DELT(L)*(A(L)*GAMM(L)+BET(L)*F(L))

    3 CONTINUE



      GAP2=0.
      DAMD2=1.
      XUP2=0.

      REL12=(RE(NZ,M)+RE(NZM1,M))/2.

      GAP2=0.
      DAMD2=1.
      XUP2=-BIX(M)

C      GAP2=0.
C      DAMD2=1.
C      XUP2=(RDY(NZM1)/RDY(NZ)-1.)*HRFI(NZM1,M)*REL12

      WL(NZ)=(GAMM(NZ)*GAP2-REL12*XUP2)/
     *(REL12*DAMD2+BET(NZ)*GAP2)


      DO 4 N=2,NZM1
      L=NZ-N+1
      WL(L)=((1.-C(L)/A(L)*BET(L))*(B(L)*WL(L+1)-F(L))+
     *C(L)*GAMM(L))/A(L)
      PROZ(L)=(A(L)*GAMM(L)+BET(L)*(F(L)-B(L)*WL(L+1)))/A(L)
    4 CONTINUE

C      PROZ(NZ)=(BET(NZ)*XUP2+GAMM(NZ)*DAMD2)/
C     *(REL12*DAMD2+BET(NZ)*GAP2)

      PROZ(NZ)=PROZ(NZM1)-WL(NZ)/REL12

      PROZ(1)=-RXAP

      DO 14 L=2,NZ
      RE12=(RE(L-1,M)+RE(L,M))/2.
      DHDZL(L)=-WL(L)/RE12/DZ
   14 CONTINUE
      DHDZL(1)=0.

      RETURN
      END


      SUBROUTINE PRO_HR(PROR,HRFI,RE,CNUR,CNUZ,CNURZ,
     *BAN1,BKA1,RAD,RDY,RDZ,DY,DZ,DTPR,L,RXAP,NZ,NR,DHDRM)
C      PARAMETER (NZ=161,NR=81)
      DIMENSION RE(NZ,NR),CNUR(NZ,NR),CNUZ(NZ,NR),CNURZ(NZ,NR)
      DIMENSION HRFI(NZ,NR),PROR(NR),RDY(NZ),RDZ(NZ,NR),RAD(NZ,NR)
      DIMENSION BAN1(NZ),BKA1(NZ),DHDRM(NR)
      DIMENSION A(NR),B(NR),C(NR),F(NR),BET(NR),
     *GAMM(NR),DELT(NR),WL(NR)
      COMMON /LKP/ LKP2
      COMMON /LS/ LST

      NRM1=NR-1

C            LKP2=42


      DO 1 M=2,NR
      C(M)=RDY(L)/DTPR
C      C(M)=RDY(L)**2/DTPR
      PAR1=RAD(L,M)*(CNUR(L,M)+CNUR(L,M-1))/(RAD(L,M)+RAD(L,M-1))
C      PAR1=0.
C      PAR1=RAD(L,M)*(CNUR(L,M)/RAD(L,M)+CNUR(L,M-1)/RAD(L,M-1))/2.
      PAR2=(RDZ(L,M)**2*CNUZ(L,M)+RDZ(L,M-1)**2*CNUZ(L,M-1))/2.
      PAR3=2.*CNURZ(L,M)*(RDZ(L,M)+RDZ(L,M-1))/2.

      A(M)=(PAR1+PAR2+PAR3)/DY**2/RDY(L)
C      A(M)=(PAR1+PAR2+PAR3)/DY**2
    1 CONTINUE


      DO 2 M=2,NRM1
      PAR1=RAD(L,M)*(CNUR(L,M)+CNUR(L,M+1))/(RAD(L,M)+RAD(L,M+1))
C      PAR1=RAD(L,M)*(CNUR(L,M)/RAD(L,M)+CNUR(L,M+1)/RAD(L,M+1))/2.
C      PAR1=0.
      PAR2=(RDZ(L,M)**2*CNUZ(L,M)+RDZ(L,M+1)**2*CNUZ(L,M+1))/2.
      PAR3=2.*CNURZ(L,M)*(RDZ(L,M)+RDZ(L,M+1))/2.
      B(M)=(PAR1+PAR2+PAR3)/DY**2/RDY(L)
C      B(M)=(PAR1+PAR2+PAR3)/DY**2

      FLYZD2=-(RDZ(L,M)*CNUZ(L,M)+CNURZ(L,M))*RE(L,M)/4./DZ/DY*
     *(HRFI(L+1,M+1)-HRFI(L+1,M-1)-HRFI(L-1,M+1)+HRFI(L-1,M-1))
      FL2=-(RE(L+1,M)*RDZ(L+1,M)*CNUZ(L+1,M)-
     *RE(L-1,M)*RDZ(L-1,M)*CNUZ(L-1,M))/4./DZ/DY*
     *(HRFI(L,M+1)-HRFI(L,M-1))
      F(M)=RDY(L)*HRFI(L,M)/DTPR+FLYZD2+FL2
C      F(M)=RDY(L)**2*HRFI(L,M)/DTPR+(FLYZD2+FL2)*RDY(L)

    2 CONTINUE

      IF(L.GE.LKP2) GOTO 715

      IF(L.GE.LST) GOTO 717

      GAP1=1.
C      DAMD1=0.
      DAMD1=2./(RE(L,1)+RE(L,2))
      XUP1=-RXAP

      GOTO 716


  717 CONTINUE

      GAP1=0.
C      DAMD1=2./(RE(L,1)+RE(L,2))
      DAMD1=1.
C      *RDY(L)
      XUP1=BKA1(L)
C      XUP1=HRFI(L,1)
      GOTO 716
  715 CONTINUE

      GAP1=1.
      DAMD1=2./(RE(L,1)+RE(L,2))

C      GAP1=1.
      DAMD1=0.
      XUP1=0.

      GAP1=0.
      DAMD1=1.
C      XUP1=(RE(L,1)+RE(L,2))/2.*(HRFI(L,1)-HRFI(L,2))

      GAP1=1.
      DAMD1=0.
      XUP1=-HRFI(L,2)
      XUP1=0.

  716 CONTINUE

      CNAM=A(2)*GAP1+C(2)*DAMD1
      BET(2)=A(2)*DAMD1/CNAM
      GAMM(2)=A(2)*XUP1/CNAM

      DO 3 M=2,NRM1
      RE12=(RE(L,M+1)+RE(L,M))/2.
      DELT(M)=A(M)*(A(M+1)*RE12+C(M+1))+RE12*B(M)*C(M+1)*BET(M)
      BET(M+1)=A(M+1)/DELT(M)*(RE12*B(M)*BET(M)+A(M))
      GAMM(M+1)=A(M+1)*RE12/DELT(M)*(A(M)*GAMM(M)+BET(M)*F(M))

    3 CONTINUE

      REM12=(RE(L,NR)+RE(L,NRM1))/2.

      GAP2=0.
      DAMD2=1.
      XUP2=-BAN1(L)
C      XUP2=HRFI(L,NR)

      WL(NR)=(GAMM(NR)*GAP2-REM12*XUP2)/
     *(REM12*DAMD2+BET(NR)*GAP2)

C      WL(NR)=BAN1(L)
C      WL(NR)=WUP


      DO 4 N=2,NRM1
      M=NR-N+1
      WL(M)=((1.-C(M)/A(M)*BET(M))*(B(M)*WL(M+1)-F(M))+
     *C(M)*GAMM(M))/A(M)
      PROR(M)=(A(M)*GAMM(M)+BET(M)*(F(M)-B(M)*WL(M+1)))/A(M)
    4 CONTINUE


      RE32=(RE(L,1)+RE(L,2))/2.

      IF(L.GE.LKP2) GOTO 725

      IF(L.GE.LST) GOTO 727

      PROR(1)=-RXAP
      GOTO 726

  727 CONTINUE

      PROR(1)=PROR(2)+WL(2)/RE32

      GOTO 726
  725 CONTINUE
      PROR(1)=0.

  726 CONTINUE


C      PROR(1)=PROR(2)+WL(2)/RE32

      PROR(NR)=PROR(NRM1)-WL(NR)/REM12

      DO 14 M=2,NR
      RE12=(RE(L,M-1)+RE(L,M))/2.
      DHDRM(M)=-WL(M)/RE12/DY
   14 CONTINUE
      DHDRM(1)=0.

      RETURN
      END

      SUBROUTINE QAPA(NZ,NR,TEM,HFI,PLT,AL,QCYM)
      DIMENSION TEM(NZ,NR),HFI(NZ,NR),PLT(NZ,NR),QCYM(NZ,NR),AL(NZ,NR)
      COMMON /BBB/ BB,TB,PB16,DLCM,GAM1,PM,ZION

      TBD1000=TB/1000.

      DO 1 L=1,NZ
      DO 1 M=1,NR

      TEMP=TEM(L,M)
      PLOT=PLT(L,M)
      POLE=ABS(HFI(L,M))
      CKOP=SQRT(BB*TEMP)
      ALFA=AL(L,M)

      FLOC=0.6/10000.*SQRT(TB)/PB16
C      DLOC(N)=FLOC*SQRT(TEMP)/PLOT/(1.-ALFA)
C      DLOC=FLOC*SQRT(TEMP)/PLOT/(1.-ALFA)
      DLOC=FLOC*SQRT(TEM(1,1))/PLT(1,1)/(1.-AL(1,1))
C      DLOCER=DLOC(1)
      FLYCH=8.52*SQRT(PM*BB**3)*TBD1000**3/SQRT(TB)/PB16/DLCM
C      QLY=FLYCH*DLOCER*TEMP**3
      QLY=FLYCH*DLOC*TEMP**3
C      QLY(N)=FLYCH*DLOCER*TEM**3
C      QLY=0.

      FWETE=0.28/1000.*TBD1000**2/SQRT(BB*PB16)
      WETE=FWETE*POLE*SQRT(TEMP**3)/PLOT
      GWETE=(11.92+4.66*WETE**2)/(3.77+14.79*WETE**2+WETE**4)
      FEL=1.62/100000.*SQRT(PM*BB**3)*TBD1000**2*GWETE/PB16/DLCM/ZION
      QELN=FEL*SQRT(TEMP**5)
C      QEL(N)=QELN

      FAT=6.76/10000.*SQRT(PM*BB**3)*TB**0.25/PB16/DLCM
      QATN=FAT*TEMP**0.75
C      QAT(N)=QATN
      QELDAT=QELN/QATN
      S1=1.156/100000.
      S2=1.44*(TEMP*TB)**0.16
      PART1=QELN*ALFA/(ALFA+S1*QELDAT*(1.-ALFA))
      PART2=QATN*(1.-ALFA)/(1.-ALFA+S2*ALFA)
C      QATEL(N)=PART1+PART2


      QCYM(L,M)=PART1+PART2+QLY

    1      CONTINUE



      RETURN
      END

      SUBROUTINE PRO_TZ(PROZ,TEM,QCYM,PLT,BIX,
     *RDY,RDZ,RAD,DY,DZ,DTPR,M,NZ,NR)

      DIMENSION QCYM(NZ,NR),PLT(NZ,NR),BIX(NR)
      DIMENSION TEM(NZ,NR),PROZ(NZ),RDY(NZ),RDZ(NZ,NR),RAD(NZ,NR)
      DIMENSION A(NZ),B(NZ),C(NZ),G(NZ),BET(NZ),
     *GAMM(NZ),DELT(NZ),WL(NZ)
      COMMON /BBB/ BB,TB,PB16,DLCM,GAM1,PM,ZION

      NZM1=NZ-1

      DO 1 L=2,NZ
      C(L)=RAD(L,M)*RDY(L)**2*PLT(L,M)*BB/GAM1/DTPR
      A(L)=(RAD(L,M)*RDY(L)**2)/DZ**2
    1 CONTINUE

      DO 2 L=2,NZM1
      B(L)=(RAD(L,M)*RDY(L)**2)/DZ**2
      FLZY=RAD(L,M)*RDY(L)**2/8./DZ/DY*(QCYM(L+1,M)*RDZ(L+1,M)/
     *RDY(L+1)*(TEM(L+1,M+1)-TEM(L+1,M-1))-
     *QCYM(L-1,M)*RDZ(L-1,M)/RDY(L-1)*(TEM(L-1,M+1)-TEM(L-1,M-1)))
      FLYZ=RAD(L,M)*RDY(L)*RDZ(L,M)/8./DY/DZ*(QCYM(L,M+1)*
     *(TEM(L+1,M+1)-TEM(L-1,M+1))-QCYM(L,M-1)*
     *(TEM(L+1,M-1)-TEM(L-1,M-1)))

      G(L)=C(L)*TEM(L,M)-FLZY-FLYZ
    2 CONTINUE

      GAP1=1.
      DAMD1=2./(QCYM(1,M)+QCYM(2,M))
      XUP1=1.
CCCCCCCCCCCCCCCCCCCCCCCCCc
      XUP1=TEM(1,M)

      CNAM=A(2)*GAP1+C(2)*DAMD1
      BET(2)=A(2)*DAMD1/CNAM
      GAMM(2)=A(2)*XUP1/CNAM

      DO 3 L=2,NZM1
      QCYM12=(QCYM(L+1,M)+QCYM(L,M))/2.
      DELT(L)=A(L)*(A(L+1)*QCYM12+C(L+1))+QCYM12*B(L)*C(L+1)*BET(L)

      BET(L+1)=A(L+1)/DELT(L)*(QCYM12*B(L)*BET(L)+A(L))
      GAMM(L+1)=A(L+1)*QCYM12/DELT(L)*(A(L)*GAMM(L)+BET(L)*G(L))
    3 CONTINUE

      GAP2=0.
      DAMD2=1.
      XUP2=0.

      GAP2=0.
      DAMD2=1.
      XUP2=-BIX(M)

      QCYML12=(QCYM(NZ,M)+QCYM(NZM1,M))/2.

C      GAP2=0.
C      DAMD2=1.
C      XUP2=(RDY(NZM1)/RDY(NZ)-1.)*TEM(NZM1,M)*QCYML12

      WL(NZ)=(GAMM(NZ)*GAP2-QCYML12*XUP2)/
     *(QCYML12*DAMD2+BET(NZ)*GAP2)
C      PROZ(NZ)=(BET(NZ)*XUP2+GAMM(NZ)*DAMD2)/
C      *(QCYML12*DAMD2+BET(NZ)*GAP2)


      DO 4 N=2,NZM1
      L=NZ-N+1
      WL(L)=((1.-C(L)/A(L)*BET(L))*(B(L)*WL(L+1)-G(L))+
     *C(L)*GAMM(L))/A(L)
      PROZ(L)=(A(L)*GAMM(L)+BET(L)*(G(L)-B(L)*WL(L+1)))/A(L)
    4 CONTINUE

      PROZ(NZ)=PROZ(NZM1)-WL(NZ)/QCYML12
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      PROZ(NZ)=PROZ(NZM1)

      PROZ(1)=1.
CCCCCCCCCCCCCCCCCCCCCCCCC
      PROZ(1)=TEM(1,M)

      RETURN
      END

      SUBROUTINE PRO_TR(PROR,TEM,QCYM,PLT,
     *RAD,RDY,RDZ,DY,DZ,DTPR,L,NZ,NR)
      DIMENSION QCYM(NZ,NR),PLT(NZ,NR)
      DIMENSION TEM(NZ,NR),PROR(NR),RDY(NZ),RDZ(NZ,NR),RAD(NZ,NR)
      DIMENSION A(NR),B(NR),C(NR),G(NR),BET(NR),
     *GAMM(NR),DELT(NR),WL(NR)
       COMMON /BBB/ BB,TB,PB16,DLCM,GAM1,PM,ZION
      COMMON /LKP/ LKP2


      NRM1=NR-1


      DO 1 M=2,NR
      C(M)=RAD(L,M)*RDY(L)**2*PLT(L,M)*BB/GAM1/DTPR
      A(M)=(RAD(L,M)*RDZ(L,M)*(RDZ(L,M-1)+RDZ(L,M))/2.+
     *(RAD(L,M-1)+RAD(L,M))/2.)/DY**2
    1 CONTINUE

      DO 2 M=2,NRM1
      B(M)=(RAD(L,M)*RDZ(L,M)*(RDZ(L,M+1)+RDZ(L,M))/2.+
     *(RAD(L,M+1)+RAD(L,M))/2.)/DY**2

      FLZY=RAD(L,M)*RDY(L)**2/8./DZ/DY*(QCYM(L+1,M)*RDZ(L+1,M)/
     *RDY(L+1)*(TEM(L+1,M+1)-TEM(L+1,M-1))-
     *QCYM(L-1,M)*RDZ(L-1,M)/RDY(L-1)*(TEM(L-1,M+1)-TEM(L-1,M-1)))
      FLYZ=RAD(L,M)*RDY(L)*RDZ(L,M)/8./DY/DZ*(QCYM(L,M+1)*
     *(TEM(L+1,M+1)-TEM(L-1,M+1))-QCYM(L,M-1)*
     *(TEM(L+1,M-1)-TEM(L-1,M-1)))

      G(M)=C(M)*TEM(L,M)-FLZY-FLYZ
    2 CONTINUE

      IF(L.GE.LKP2) GOTO 7177
      GAP1=1.
      DAMD1=0.
      XUP1=TEM(L,2)
CCCCCCCCCCCCCCCCCCCCCCC
       GAP1=0.
      DAMD1=1.
      XUP1=0.

      GOTO 7178
 7177 CONTINUE

       GAP1=0.
      DAMD1=1.
      XUP1=0.


 7178 CONTINUE

      CNAM=A(2)*GAP1+C(2)*DAMD1
      BET(2)=A(2)*DAMD1/CNAM
      GAMM(2)=A(2)*XUP1/CNAM


      DO 3 M=2,NRM1
      QCYM12=(QCYM(L,M+1)+QCYM(L,M))/2.
      DELT(M)=A(M)*(A(M+1)*QCYM12+C(M+1))+QCYM12*B(M)*C(M+1)*BET(M)
      BET(M+1)=A(M+1)/DELT(M)*(QCYM12*B(M)*BET(M)+A(M))
      GAMM(M+1)=A(M+1)*QCYM12/DELT(M)*(A(M)*GAMM(M)+BET(M)*G(M))
    3 CONTINUE

      QCYMM12=(QCYM(L,NR)+QCYM(L,NRM1))/2.

      GAP2=0.
      DAMD2=1.
      XUP2=0.


      WL(NR)=(GAMM(NR)*GAP2-QCYMM12*XUP2)/
     *(QCYMM12*DAMD2+BET(NR)*GAP2)
C      PROR(NR)=(BET(NR)*XUP2+GAMM(NR)*DAMD2)/
C      *(QCYMM12*DAMD2+BET(NR)*GAP2)

      DO 4 N=2,NRM1
      M=NR-N+1
      WL(M)=((1.-C(M)/A(M)*BET(M))*(B(M)*WL(M+1)-G(M))+
     *C(M)*GAMM(M))/A(M)
      PROR(M)=(A(M)*GAMM(M)+BET(M)*(G(M)-B(M)*WL(M+1)))/A(M)
    4 CONTINUE

      QCYM32=(QCYM(L,1)+QCYM(L,2))/2.
      PROR(1)=PROR(2)+WL(2)/QCYM32
      PROR(NR)=PROR(NRM1)-WL(NR)/QCYMM12

      PROR(1)=PROR(2)
      PROR(NR)=PROR(NRM1)

      RETURN
      END


      SUBROUTINE RZEKOP(AI,BI,XKOP,EPS)
      COMMON /NUMITER/ NITER
      NITER=1
      NKPIT=100
      X1=AI
      X2=BI
      FA=FCN(X1)
      FB=FCN(X2)
C  174 DIF=ABS(X2-X1)/ABS(X1)
  174 DIF=ABS(X2-X1)
      NITER=NITER+1
      IF(NITER.GE.NKPIT) GOTO 183
      IF(DIF-EPS) 181,181,182
  181 XKOP=X1+(X2-X1)/2.
      GOTO 180
  182 CONTINUE
      IF(FA*FB-0.) 171,172,277
  172 IF(FA-0.) 176,175,176
  176 XKOP=X2
      GOTO 180
  175 XKOP=X1
      GOTO 180
  171 X3=X1+(X2-X1)/2.
      FX3=FCN(X3)
      IF(FA*FX3-0.) 177,178,179
  178 XKOP=X3
      GOTO 180
  177 X2=X3
      FB=FX3
      GOTO 174
  179 X1=X3
      FA=FX3
      GOTO 174
  180 CONTINUE
      GOTO 184

  277      PRINT 7171
 7171 FORMAT(1X,'NO SOLUTION')
C      XKOP=X1+(X2-X1)/2.

      PRINT 112,NITER
  112 FORMAT(3X,'NITER',I4)

      GOTO 184

  183      PRINT 8171
 8171 FORMAT(1X,'in RZEKOP NITER=NKPIT')


  184 CONTINUE

      RETURN
      END

      FUNCTION FCN(X)

      COMMON /BBTKP/ BB,TKP,DABLM,PLOT,CAX

C      TEMP=2.*DABLM/BB/(1.+X)/PLOT
      TEMP=X
      TKPIT=70.
C      X1=-TKP/TEMP
C      X2=X1*X1/2.E0
C      X3=X2*X1/3.E0
C      X4=X3*X1/4.E0
C      X5=X4*X1/5.E0
C      X6=X5*X1/6.E0
C      X7=X6*X1/7.E0
C      X8=X7*X1/8.E0
C      TAS=CAX*TEMP**1.5*(1.E0+X1+X2+X3+X4+X5+X6+X7+X8)
C       TAS=CAX*TEMP**1.5*EXP(-TKP/TEMP)
C       TASP=TAS/PLOT
C      QTAS=SQRT(TASP)
C      ALF2=-TASP/2.E0+QTAS*SQRT(TASP/4.E0+1.E0)
      TKPFKT=TKP/TEMP
      IF(TKPFKT.GT.TKPIT) TKPFKT=TKPIT

      CAB=CAX*BB/2.E0/DABLM
      TAB=CAB*TEMP**2.5*EXP(-TKPFKT)
      TASB=TAB/(1.E0+TAB)
      ATAS=SQRT(TASB)
      ALF2=ATAS

C      FCN=X-ALF2
C      FCN=DABLM-BB/2.E0*(1.E0+ALF2)*PLOT*X
      FCN=1.E0-BB/2.E0*(1.E0+ALF2)*PLOT*TEMP/DABLM

C      BP=2.E0*DABLM/BB/PLOT
C      FCN=(1.E0-2.E0*TASP)*X**2-BP*X*(2.E0-TASP)+BP**2
C      TB=TAS*X*BB/2.E0/DABLM
C      FCN=X**2-2.E0*BP*X*(1.E0+TB)+(1.E0+TB)*BP**2

C      FCN=DABLM*2.E0/BB/PLOT-(1.E0+ALF2)*X

      RETURN
      END

      SUBROUTINE FEXT (NEQ, T, Y, YDOT)
      PARAMETER (NH=101)
      DOUBLE PRECISION T, Y, YDOT,TKPD,DUHMD
      DOUBLE PRECISION RE1,RE2,PART2,T3,FK1,FK2,FK3,FK4,ZIOND,RELMD
      DOUBLE PRECISION UPLD,PLELD,UHMD,BIND,APEKD,ALD,QFOTD
      DOUBLE PRECISION SRE1D,SRE2D,HXD,BBD,GAMD,DZETAD,PART3
      DIMENSION Y(NH), YDOT(NH),DUHMD(NH)
      COMMON /PAR/ SRE1D,SRE2D,HXD,BBD
      COMMON /PART/ GAMD,DZETAD,TKPD
      COMMON /BBB/ BB,TB,PB16,DLCM,GAM1,PM,ZION
C      DIMENSION UPLD(NH),PLELD(NH),TEMPD(NH),UCKD(NH),ECED(NH)
      DIMENSION UPLD(NH),PLELD(NH),FOTIND(NH),FOTPED(NH),QFOTD(NH)
      DOUBLE PRECISION FOTIND,FOTPED                                    ! <!> KV
C      DIMENSION UHMD(NH),BIND(NH),APEKD(NH),DEDXD(NH),RED(NH),ALD(NH)
      DIMENSION UHMD(NH),BIND(NH),APEKD(NH),ALD(NH)
C      COMMON /MAC/ UPLD,PLELD,TEMPD,UHMD,BIND,APEKD
      COMMON /MAC/ UPLD,PLELD,UHMD,BIND,APEKD
      COMMON /MACF/ FOTIND,FOTPED
C      COMMON /MACT/ UCKD,DEDXD,ECED
      COMMON /VIC/ ALD
      COMMON /ABCYM/ EX,ACYM,BCYM,E1,VF
      COMMON /DHD/ DUHMD
      COMMON /QB/ QB2

C      NH=NEQ
      NHM1=NH-1
      DO 768 N=1,NH
      TEMPN=SNGL(Y(N))
      TEB=TEMPN/EX
      PLOTN=SNGL(UPLD(N))

      BPR0=DLCM/VF*PB16*5.6/SQRT(10.)/10**8
      EA1=2.334733
      EA2=0.250621
      EB1=3.330657
      EB2=1.681534
      APR0=0.332462*BPR0/10**5

      XE=E1/TEB
      E12X=(XE**2+EA1*XE+EA2)/
     *(XE**2+EB1*XE+EB2)/XE
      E1X=EXP(-XE)*E12X
      BPRF=SQRT(TEB)*E1X
      BPR=BPRF*BPR0
      APR=APR0*E12X/TEB


      KCYM=20
      CYM=0.
      DO 1517 K=1,KCYM

      FLK=0.3*ALOG(1.017+0.462*TEB/E1*K**2*(K+1.)**2/(2.*K+1.))

      X=TEB*K**2*(K+1)**2/E1/(2*K+1)
      CYM=CYM+2.*(2*K+1)*EXP(-E1/TEB/(K+1)**2)/K**3/(K+1)**4/FLK
C      /S(X)
 1517 CONTINUE

      APEKN=ACYM/CYM/TEB**2+APR
      BINN=BCYM*EXP(-E1/TEB)/SQRT(TEB)/CYM+BPR
      BIND(N)=DBLE(BINN)
      APEKD(N)=DBLE(APEKN)
      FTINN=0.312*DLCM/VF*E1**2*EXP(-E1/TEB)*TEB
      FTPEN=6.457*DLCM*PB16/VF*E1**2*SQRT(TEB)/10**7
      FOTIND(N)=DBLE(FTINN)
      FOTPED(N)=DBLE(FTPEN)

C      TEB=TEM(L,M)*TB/1.16/10**4
      FQ=QB2*PB16**2*PLOTN**2/10**5
C      uzing the encyclopedia for T < 4 eV
      QLIN=6.25*FQ*ZION**4
C      uzing the paper by Kogan
      QLIN=80.*FQ*ZION**7/SQRT(TEB**3)

      QPEK=4.4*FQ*ZION**5/SQRT(TEB)
      QTOP=0.154*FQ*ZION**3*SQRT(TEB)

      QFOT=QLIN+QPEK+QTOP
      QFOTD(N)=DBLE(QFOT)
C      QFOT=QPEK+QTOP
C      FGB=2.*(GAM-1.)/BB*RDY(L)/TEM(L,M)/(1.+AL(L,M))

C      QCB1=QFOT*FGB*RAD(L,M)
C      QCB1=0.

  768 CONTINUE

      ZIOND = DBLE(ZION)                      ! <!>
      DO 1 N=2,NHM1
      T3=Y(N)
      RE1=(UPLD(N)-PLELD(N))/PLELD(N)/SRE1D
      RE2=1.D0/SRE2D/DSQRT(T3*T3*T3)
      RELMD=ZIOND*(RE1+RE2)
      FK1=(UPLD(N)-PLELD(N))*PLELD(N)*BIND(N)
      FK2=PLELD(N)*PLELD(N)*PLELD(N)*APEKD(N)
      FK3=(UPLD(N)-PLELD(N))*FOTIND(N)
      FK4=PLELD(N)*PLELD(N)*FOTPED(N)
C      PART2=2.D0/BBD/UPLD(N)/(1.D0+ALD(N))*(RE1+RE2)*(GAMD-1.D0)*
C     *DUHMD(N+1)/HXD*
C     *DUHMD(N+1)/HXD
C     *(UHMD(N+1)-UHMD(N-1))/HXD/2.D0*
C     *(UHMD(N+1)-UHMD(N-1))/HXD/2.D0
C     *DUHMD(N)/HXD*
C     *DUHMD(N)/HXD
C     *(DUHMD(N)+DUHMD(N+1))/2.D0/HXD*
C     *(DUHMD(N)+DUHMD(N+1))/2.D0/HXD

      PART2=0.D0
      PART2=2.D0/BBD/UPLD(N)/(1.D0+ALD(N))*(GAMD-1.D0)*QFOTD(N)
      PART3=2.D0/BBD/UPLD(N)/(1.D0+ALD(N))*(GAMD-1.D0)*RELMD*DUHMD(N)
      PART2=0.D0
      PART3=0.D0

      YDOT(N)=PART3-(FK1-FK2+FK3-FK4)*(GAMD-1.D0)*DZETAD*TKPD/UPLD(N)/
     *(1.D0+ALD(N))-PART2
C       YDOT(N)=PART3
C       YDOT(N)=0.D0
    1 CONTINUE
      YDOT(1)=0.D0
      YDOT(NH)=YDOT(NHM1)

      RETURN
      END

      SUBROUTINE JEXT (NEQ, TD, Y, ML, MU, PD, NRPD)
      DOUBLE PRECISION PD, TD, Y
      DIMENSION Y(NEQ), PD(NRPD,NEQ)
      RETURN
      END

*      ! ------------- get array from scalar grid function
*      subroutine get_array_from_scalar_gf(SCL,NZ,NR,ARRAY)
*      use module_GridFunct
*      implicit none
*      ! input:
*      type(TGridNodalScalar) :: SCL          !
*      integer                :: NZ,NR
*      ! output:
*      real, dimension(NZ,NR) :: ARRAY
*      ! local:
*      integer :: i,j,inode
*
*      inode = 0
*      do i = 1,NZ
*          do j = 1,NR
*              inode = inode + 1
*              ARRAY(i,j) = real(SCL%fval(inode),4)
*          enddo
*      enddo
*
*      return
*      end subroutine get_array_from_scalar_gf
*
*      ! ------------- get array from vector grid function component
*      subroutine get_array_from_vector_gf(VEC,iCmp,NZ,NR,ARRAY)
*      use module_GridFunct
*      implicit none
*      ! input:
*      type(TGridNodalVector) :: VEC          !
*      integer                :: iCmp
*      integer                :: NZ,NR
*      ! output:
*      real, dimension(NZ,NR) :: ARRAY
*      ! local:
*      integer :: i,j,inode
*
*      inode = 0
*      do i = 1,NZ
*          do j = 1,NR
*              inode = inode + 1
*              ARRAY(i,j) = real(VEC%fvec(inode)%cmp(iCmp),4)
*          enddo
*      enddo
*
*      return
*      end subroutine get_array_from_vector_gf
