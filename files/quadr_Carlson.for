      SUBROUTINE Carlsonquadr(Lmax,Um,Smu,WG)
c Программа расчета угловой сетки Карлсона на четверти сферы
c Сетка полностью определяется параметром Lmax
c и состоит из:
c сетки по гамма - косинусу полярного угла theta (Lmax узлов)
c и сеток по мю - косинусу азимутального угла fi
c Для каждого слоя по гамма используется своя сетка по мю
c Все веса квадратуры здесь одинаковые и их сумма равна pi
      parameter (pi=3.1415926)
      INTEGER Lmax
      REAL(4) :: WG
      DIMENSION UM(LMAX),CFI(LMAX*(LMAX+1)/2),WL(Lmax),UMP(Lmax),
     *UM2(Lmax),Smu(Lmax*(Lmax+1))
      call SN (LMAX,UM,CFI,WL,UMP,UM2,Smu,WG,pi)
c      write (*,*) 'N_th =',LMAX
c      write (*,*) 'N_angle_mesh =',LMAX*(LMAX+1)
c      write (*,*)
c      write (*,*) 'Weight of the quadrature =',wg
c      write (*,*) 'Points summ of the quadrature =',pi
c      write (*,*)
c      write (*,*) 'cos(th)'
c      write (*,'(6F11.8)') Um
c      write (*,*)
c      write (*,*) 'cos(fi)'
      i1=1
      Do 12 k=1,Lmax
      i2=2*(Lmax-k+1)+i1-1
c      write (*,*) 'Number of the theta layer =',k
c      write (*,'(6F11.7)') (smu(i),i=i1,i2)
12    i1=i2+1
      return
      end


      SUBROUTINE SN (LMAX,UM,CFI,WL,UMP,UM2,Smu,WG,pi)
      DIMENSION UM(LMAX),CFI(LMAX*(LMAX+1)/2),WL(Lmax),UMP(Lmax),
     *UM2(Lmax),Smu(Lmax*(Lmax+1))
c  Расчет узлов и весов ESn-квадратуры
      EPS=0.01
c  сумма весов равна pi
      WG=pi/float(LMAX*(LMAX+1))
c  окончание вычисления весов
c  расчет узлов по гамма
      E0=0.0
      E1=1.0
      CALL MUMESH(UM,WL,LMAX,UMP,UM2,E0,E1)
c  расчет узлов по фи
      CALL FIMESH(CFI,LMAX,UM,WL,EPS,pi)
      K=0
      DO 2 I=1,LMAX
      JMAX=LMAX-I+1
      DO 2 J=1,JMAX
      K=K+1
      RAB=COS(CFI(K))
    2 CFI(K)=RAB

      N=1
      i1=1
      Do 17 k=1,Lmax
      i2=Lmax-k+i1
      j=i2
      Do 4 m=1,Lmax-k+1
      Smu(N)=-CFI(j)
      j=j-1
4     N=N+1
      j=i1
      Do 8 m=1,Lmax-k+1
      Smu(N)=CFI(j)
      j=j+1
8     N=N+1
17    i1=i2+1


      RETURN
      END


      SUBROUTINE MUMESH(UM,W,LMAX,UMP,UM2,Gm0,Gp0)
      DIMENSION UM(LMAX),UMP(LMAX),UM2(LMAX),W(LMAX)
c  расчет узлов по гамма
      N=LMAX*2
      A=0.
      B=0.
      C=0.
      Gp=amax1(abs(Gm0),abs(Gp0))
      Gm=amin1(abs(Gm0),abs(Gp0))
      DN=float(N*(N+2))/abs(Gp-Gm)
      DO 1 L=1,LMAX
      DL=float(N-2*L+2)
      UMP(L)=Gp-DL*DL/DN
      UM2(L)=Gp-DL*(DL+2)/DN
      W(L)=4.0*DL/DN
      DX=W(L)*UM2(L)
      A=A+DX*UM2(L)
      B=B+DX*UMP(L)
    1 C=C+W(L)*UMP(L)*UMP(L)
c  параметр f для счета центров
      F0=(Gp*Gp*Gp-Gm*Gm*Gm)/3.0
      F=(-B+SQRT(B*B-A*(C-F0)))/A
      DO 2 L=1,LMAX
    2 UM(L)=UMP(L)+F*UM2(L)
C*****

      S=0.
      DO 3 I=1,LMAX
    3 S=S+W(I)*UM(I)*UM(I)
c      S=0.0
c      Do L=1,Lmax
c      S=s+w(L)
c      write (111,'(I3,E12.4)') L,S
c      end Do
      RETURN
      END

      SUBROUTINE FIMESH(FI,LMAX,UM,W,EPS,pi)
      DIMENSION FI(LMAX*(LMAX+1)/2),W(LMAX),UM(LMAX)
c  Вычисление узлов по фи
      N=LMAX*2
      C=0.
      CK=N*(N+2)*SQRT(2.)/8.
      DO 1 L=1,LMAX
    1 C=C+W(L)*UM(L)
      C=C*CK
      A0=1.
      A1=1.+0.89/N
c  вычисление An
    4 CALL FUN(A0,LMAX,UM,F0)
      F0=F0-C
c  вычисление An
      CALL FUN(A1,LMAX,UM,F1)
      F1=F1-C
      A2=A1-F1*(A1-A0)/(F1-F0)
c  вычисление An
      CALL FUN(A2,LMAX,UM,F2)
      IF(ABS(F2-C)-EPS)3,3,2
    2 A0=A1
      A1=A2
      GO TO 4
    3 FK=(1.-A2)/2.
      K=0
      DO 5 L=1,LMAX
      MMAX=LMAX-L+1
      NL=N-2*L+2
      DO 5 M=MMAX,1,-1
      FF=A2*(2*M-1)/NL+FK
      K=K+1
    5 FI(K)=PI*FF/2.

      RETURN
      END

      SUBROUTINE FUN(A,LMAX,UM,F)
      DIMENSION UM(LMAX)
c  решение ур-ия  F(An)-C=0
      N=LMAX*2
      F=0.
      B=3.1415926*A/4.
      DO 1 L=1,LMAX
      C=B*2./(N-2*L+2)
      C=SIN(C)
      C=SQRT(1-UM(L)*UM(L))/C
    1 F=F+C
      F=F*SIN(B)
      RETURN
      END


