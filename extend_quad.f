!-----CALCULATES BOUND STATE OF PS-
!
!-----AUTHOR: WILLIAM MITCHELL (UNIVERSITY OF NORTH TEXAS)
!
      IMPLICIT REAL * 16 (A-H,O-Z)
      REAL * 16 FIHFI,FIFI,RK,RL,RM,TP,SP,R3P,XSQ,CSQ,ALPHA2,GAMMA2,PI
      REAL * 16 C,EPS3,CJB4
      COMMON /CJB4/FIHFI(364,364),FIFI(364,364),FIHFIN(364,364)
      COMMON RK(364),RL(364),RM(364)
      COMMON SP(12),TP(12),R3P(12)
      COMMON XSQ(15),CSQ(15)
      COMMON ALPHA2,GAMMA2,PI
      COMMON K(364),L(364),M(364)
      COMMON N1F,N2F,IHDPTF,N1PN2F
      COMMON/EPS3/C(203)

      ! YOUR INPUT FILE IS PS.DATA, OUTPUT IS PS.TXT
      OPEN(UNIT=5,FILE='Ps.data',STATUS='old')
      OPEN(UNIT=6,FILE='Ps.txt',STATUS='unknown')

      ! DESCRIPTIONS OF FUNCTIONS ARE WITHIN THEM
      CALL CONST
      CALL SQUARE

      READ(5,4) IPRINT
    4 FORMAT(I2)
      WRITE(6,43)IPRINT
   43 FORMAT(3X,"IPRINT=",I5)
      IF(IPRINT.NE.0) GO TO 5
      WRITE(6,30)
   30 FORMAT(" FIHFI")
      DO 32 I=1,N1PN2F
      WRITE(6,33) I,(FIHFI(I,J),J=1,N1PN2F)
   33 FORMAT(I5/(10D12.4))
   32 CONTINUE 
      WRITE(6,31)
   31 FORMAT(" FIFI")
      DO 34 I=1,N1PN2F
      WRITE(6,33) I,(FIFI(I,J),J=1,N1PN2F)
   34 CONTINUE
    5 CONTINUE

      CALL BOUND

      STOP

      END

      SUBROUTINE CONST
              !THIS SUBROUTINE DOES A LOT OF NUMERICAL WORK BEFORE
              !BEGINNING THE CALCULATION OF THE VARIOUS MATRIX ELEMENTS

              IMPLICIT REAL * 16 (a-h,o-z)
              REAL * 16 FIHFI,FIFI,RK,RL,RM,SP,R3P,TP,XSQ,CSQ,ALPHA2
              REAL * 16 C,EPS3,CJB4,GAMMA2,PI
              COMMON /CJB4/FIHFI(364,364),FIFI(364,364),FIHFIN(364,364)
              COMMON RK(364),RL(364),RM(364)
              COMMON SP(12),TP(12),R3P(12)
              COMMON XSQ(15),CSQ(15)
              COMMON ALPHA2,GAMMA2,PI
              COMMON K(364),L(364),M(364)
              COMMON N1F,N2F,IHDPTF,N1PN2F
              COMMON/EPS3/C(203)

              ! THE PRIMARY PURPOSE OF THIS SUBROUTINE IS TO GENERATE
              ! THE POWERS OF THE SHORT-RANGE EXPANSION TERMS, AND THEN
              ! ORDER THEM
              DO 28 I=1,364
              K(I)=0.
              L(I)=0.
              M(I)=0.
              RK(I)=0.0Q0
              RL(I)=0.0Q0
              RM(I)=0.0Q0
   28         CONTINUE

              PI=4.Q0*ATAN(1.Q0)

    1         FORMAT(2F8.5)
              READ(5,*) ALPHA2,GAMMA2
              READ(5,2) IHDPTF
    2         FORMAT(I3)
              WRITE (6,10) ALPHA2,GAMMA2
   10         FORMAT (" A2=",F9.5,1x," G2=",F9.5)
              WRITE(6,11) IHDPTF
   11         FORMAT(" DEGREE OF POLYNOMIAL=",I3)
              IHDPP1=IHDPTF+1 
              WRITE(6,15)IHDPP1
   15         FORMAT(3X,"IHDPP1=",I3)

              IND=0

              DO 65 I=1,IHDPP1
              DO 66 IK=1,I
              ILMAX=I-IK+1.
              DO 67 JK=1,ILMAX
              IND=IND+1.
              I1=IK-1.
              I2=JK-1.
              I12=I-1.-I1-I2
              K(IND)=I1
              L(IND)=I2
              M(IND)=I12
              RK(IND)=K(IND)
              RL(IND)=L(IND)
              RM(IND)=M(IND)
   67         CONTINUE
   66         CONTINUE
   65         CONTINUE

              NT=IND
              SP(1)=1.0Q0
              TP(1)=1.0Q0
              R3P(1)=1.0Q0

              PRINT *, "K-L-M"
              WRITE(6,72)(K(I),L(I),M(I),I,I=1,NT)
   72         FORMAT(4I5)

              KC=0.

              DO 100 I=1,NT
              IRES=MOD(L(I),2)
              If (IRES.EQ.0) Then
                      KC=KC+1.
                      K(KC)=K(I)
                      L(KC)=L(I)
                      M(KC)=M(I)
                      RK(KC)=RK(I)
                      RL(KC)=RL(I)
                      RM(KC)=RM(I)
              END IF
  100         CONTINUE

              KK=KC
              N1F=KK
              N2F=0.
              N1PN2F=N1F+N2F

              WRITE(6,12)N1F,N2F,N1PN2F
   12         FORMAT(3X,"N1F=",I3,"N2F=",I3,3X,"N1PN2F=",I3)
              WRITE(6,72)(K(KC),L(KC),M(KC),KC,KC=1,KK)
              RETURN
              END

      SUBROUTINE SQUARE
              !THIS SUBROUTINE CALCULATES SHORT-RANGE--SHORT-RANGE
              !MATRIX ELEMENTS
              IMPLICIT REAL * 16 (a-h,o-z)
              REAL * 16 FIHFI,FIFI,RK,RL,RM,TP,SP,R3P,XSQ,CSQ,ALPHA2
              REAL * 16 C,EPS3,CJB4,GAMMA2,PI
              REAL * 16 XA,XB,XC,EXFAC,CX,CY,CZ,CXINV,CYINV,CZINV,X,Y
              REAL * 16 Z,R1,R2,R3,S,T,RECR1,RECR2,RECR3,RECS,COEF 
              REAL * 16 ST,RECST,R3T,SR3,RECS2,RECR32,STR31,STR32,SR38
              REAL * 16 FIJCF,TERM27,TERM67,TERM11,TERM1,TERM2,TERM3
              REAL * 16 TERM4,TERM55,TERM5,TERM6,TERM7,TERM8,POTTOL
              REAL * 16 TERM
              COMMON /CJB4/FIHFI(364,364),FIFI(364,364),FIHFIN(364,364)
              COMMON RK(364),RL(364),RM(364)
              COMMON SP(12),TP(12),R3P(12)
              COMMON XSQ(15),CSQ(15)
              COMMON ALPHA2,GAMMA2,PI
              COMMON K(364),L(364),M(364)
              COMMON N1F,N2F,IHDPTF,N1PN2F
              COMMON/EPS3/C(203)
              Dimension FI(364)

              DO 14 I=1,364
              DO 16 J=1,364
              FIFI(I,J)=0.0Q0
              FIHFI(I,J)=0.0Q0
   16         CONTINUE
   14         CONTINUE

              NS=IHDPTF+2

              CALL LWTS(NS,XSQ,CSQ,1.0Q0,0.0Q0)

              PRINT *, "Abscissas and Weights:"

              WRITE(6,13)( XSQ(I),CSQ(I),I,I=1,NS)
   13         FORMAT(3X,2D20.10,2X,I2)

              !BELOW IS THE INTEGRATION CODE
              XA=ALPHA2+ALPHA2
              XB=XA
              XC=GAMMA2+GAMMA2
              EXFAC=0.25Q0*PI
    4         CX=0.5Q0*(XA+XB)
              CY=0.5Q0*(XB+XC)
              CZ=0.5Q0*(XC+XA)
              CXINV=1.0Q0/CX
              CYINV=1.0Q0/CY
              CZINV=1.0Q0/CZ
              EXFAC=EXFAC*CXINV*CYINV*CZINV

              DO 5 IX=1,NS
              X=XSQ(IX)*CXINV
              DO 6 IY=1,NS
              Y=XSQ(IY)*CYINV
              DO 7 IZ=1,NS
              Z=XSQ(IZ)*CZINV
              R1=0.5Q0*(X+Z)
              R2=0.5Q0*(X+Y)
              R3=0.5Q0*(Y+Z)
              S=R1+R2
              T=R1-R2
              RECR1=1.0Q0/R1
              RECR2=1.0Q0/R2
              RECR3=1.0Q0/R3
              RECS=1.0Q0/S

              DO 8 I=2,12
              SP(I)=SP(I-1)*S
              TP(I)=TP(I-1)*T
              R3P(I)=R3P(I-1)*R3
    8         CONTINUE

              COEF=EXFAC*CSQ(IX)*CSQ(IY)*CSQ(IZ)*R1*R2*R3

              DO 11 IA=1,N1F
              FI(IA)=SP(K(IA)+1)*TP(L(IA)+1)*R3P(M(IA)+1)
   11         CONTINUE

              ST=S*S-T*T
              RECST=1.0Q0/ST
              R3T=R3*R3-T*T
              SR3=S*S-R3*R3
              RECS2=RECS*RECS
              RECR32=RECR3*RECR3
              STR31=3.0Q0*S*S-T*T-2.0Q0*R3*R3
              STR32=S*S-3.0Q0*T*T+2.0Q0*R3*R3
              SR38=8.0Q0*S*R3

              DO 9 JA=1,N1F
              FIJCF=SP(K(JA)+1)*TP(L(JA)+1)*R3P(M(JA)+1)*COEF
              TERM27=RK(JA)*RECS-ALPHA2
              TERM67=RM(JA)*RECR3-GAMMA2
              TERM11=ALPHA2*ALPHA2-2.0Q0*ALPHA2*RK(JA)*RECS
     1              +K(JA)*(K(JA)-1.0Q0)*RECS2
              TERM1=STR31*R3*TERM11
              TERM2=SR38*TERM27

              IF(L(JA).EQ.0.OR.L(JA).EQ.1)THEN
                      TERM3=0.0Q0
              ELSE
               TERM3=-4.0Q0/3.0Q0*RECST*RECR3*STR32*R3*
     1               RL(JA)*(RL(JA)-1.0Q0)*SP(K(JA)+1)*TP(L(JA)-1)*
     1               R3P(M(JA)+1)*COEF
              END IF

              TERM4=-8.0Q0*R3*RL(JA)
              TERM55=GAMMA2*GAMMA2-2.0Q0*GAMMA2*RM(JA)*RECR3+RM(JA)*
     1               (RM(JA)-1.0Q0)*RECR32
              TERM5=ST*R3*TERM55
              TERM6=2.0Q0*ST*TERM67
              TERM7=R3T*2.0Q0*S*TERM67*TERM27
              TERM8=SR3*2.0Q0*L(JA)*TERM67
              POTTOL=-(RECR1+RECR2-RECR3)
              TERM=4.0Q0/3.0Q0*(-RECST*RECR3*(TERM1+TERM2+TERM4+
     1        TERM5+TERM6+TERM7+TERM8)+POTTOL+0.25Q0)

              DO 10 IA=1,JA
              FIHFI(IA,JA)=FIHFI(IA,JA)+FI(IA)*(TERM*FIJCF+TERM3)
              FIFI(IA,JA)=FIFI(IA,JA)+FI(IA)*FIJCF
   10         CONTINUE

    9         CONTINUE
    7         CONTINUE
    6         CONTINUE
    5         CONTINUE

              DO 20 JA=1,N1F
              DO 21 IA=1,JA
              FIHFI(IA,JA)=0.75Q0*FIHFI(IA,JA)
   21         CONTINUE
   20         CONTINUE

              RETURN

              END

      SUBROUTINE LWTS(N,Y,W,A,B)
              !THIS IS A HELPER-SUBROUTINE THAT HELPS WITH THE
              !GAUSS-LAGUERRE QUADRATURE DONE IN OTHER SUBROUTINES
              IMPLICIT REAL * 16 (a-h,o-z)
              REAL * 16 FIHFI,FIFI,RK,RL,RM,TP,SP,R3P,XSQ,CSQ,ALPHA2
              REAL * 16 C,EPS3,CJB4,GAMMA2,PI
              REAL * 16 Y,W,A,B,C1,C2,X,DX,F,D,CAB
              Dimension Y(N), W(N)

     0DATA D1, D3, D15, D19, D24, D25, D255/1.0Q0, 3.0Q0, 
     1         15.0Q0, 1.9Q0, 2.4Q0, 2.5Q0, 2.55Q0/ 

              FN=N
              C1=D1/D19
              C2=D255/D19

              DO 1 I=1,N
              IF(I.GT.2) GO TO 5
              GO TO (3,4),I 
    3         X=D3/(D1+D24*FN)
              GO TO 6

    4         X=D15/(D25*FN+D1)+Y(I-1)
              GO TO 6

    5         X=(C1/I+D1+C2)*Y(I-1)-(C1/I+C2)*Y(I-2)

    6         CALL LN(F,D,X,N)

              DX=F/D
              X=X-DX
              
              IF (ABS(DX).LE.1.0Q-12) GO TO 7
              GO TO 6

    7         W(I)=D1/(X*D*D)
              Y(I)=X
    1         CONTINUE

              CAB=EXP(-A*B)/A

              DO 8 I=1,N
              Y(I)=Y(I)/A+B
              W(I)=CAB*W(I)
    8         CONTINUE

              RETURN

              END

      SUBROUTINE LN(F, D, X, N)
              !THIS IS A SECOND HELPER SUBROUTINE THAT HELPS WITH THE
              !GAUSS-LAGUERRE QUADRATURE
              IMPLICIT REAL * 16 (a-h,o-z)
              REAL * 16 FIHFI,FIFI,RK,RL,RM,TP,SP,R3P,XSQ,CSQ,ALPHA2
              REAL * 16 C,EPS3,CJB4,GAMMA2,PI
              REAL * 16 L,F,D,X
              Dimension L(3)

              L(1)=0.0Q0
              L(2)=1.0Q0

              DO 1 I=1,N
              L(3)=((2*I-1-X)*L(2)-(I-1)*L(1))/I
              L(1)=L(2)
              L(2)=L(3)
    1         CONTINUE

              F=L(2)
              D=N/X*(F-L(1))

              RETURN

              END

      SUBROUTINE BOUND
              !THIS IS THE SUBROUTINE THAT CALCULATES THE BINDING ENERGY
              IMPLICIT REAL * 16 (a-h,o-z)
              REAL * 16 FIHFI,FIFI,RK,RL,RM,TP,SP,R3P,XSQ,CSQ,ALPHA2
              REAL * 16 C,EPS3,CJB4,GAMMA2,PI
              REAL * 8 FIFI8,FIHFI8,WORK1,WORK,R
              INTEGER IA,IB,N,LWORK,INFO
              COMMON /CJB4/FIHFI(364,364),FIFI(364,364),FIHFIN(364,364)
              COMMON RK(364),RL(364),RM(364)
              COMMON SP(12),TP(12),R3P(12)
              COMMON XSQ(15),CSQ(15)
              COMMON ALPHA2,GAMMA2,PI
              COMMON K(364),L(364),M(364)
              COMMON N1F,N2F,IHDPTF,N1PN2F
              COMMON/EPS3/C(203)
              Dimension R(203),V(203,203),DL(203),E(203)
              Dimension FIFI8(364,364),FIHFI8(364,364)
              Dimension WORK1(1)
              Allocatable WORK(:)

              DO 100 I=1,203
              R(I)=0.0
              DL(I)=0.0
              E(I)=0.0

              DO 200 J=1,203
              V(I,J)=0.0
  200         CONTINUE
  100         CONTINUE

              WRITE(6,50)N1F
   50         FORMAT(3X,"The Square Matrix is of Order N1F=",I4)

              IA=INT(364)
              IB=INT(364)
              N=INT(N1F)
              IV=203
              IFAIL=0
              INFO=1

              ! THE ROUTINE WE USE IS ONLY DOUBLE-PRECISION, WHICH
              ! NECESSITATES THIS STEP
              FIHFI8 = FIHFI
              FIFI8 = FIFI

              !WE USE THE LAPACK DSYGV SUBROUTINE FROM LAPACK TO SOLVE
              !THE EIGNEVALUE PROBLEM
              LWORK = -1
              PRINT *, "LWORK pre CALL 1=", LWORK
              PRINT *, "SHAPE of WORK1 pre CALL 1=", SHAPE(WORK1)
              CALL DSYGV(INT(1),'V','U',N,FIHFI8,IA,FIFI8,IB,R,WORK1,
     1        LWORK,INFO)
              PRINT *, "Optimal LWORK from WORK1(1)=", WORK1(1)
              LWORK = WORK1(1)
              ALLOCATE (WORK(LWORK))
              PRINT *, "LWORK post CALL 1=", LWORK
              PRINT *, "SHAPE of WORK post CALL 1=", SHAPE(WORK)

              CALL DSYGV(INT(1),'V','U',N,FIHFI8,IA,FIFI8,IB,R,WORK,
     1        LWORK,INFO)

              WRITE(6,400)INFO
  400         FORMAT(3X,"INFO=",I5)
              WRITE(6,15)
   15         FORMAT(3X,"Eigenvalues in Ascending Order")

              DO 20 I=1,N1F
              WRITE(6,25)I,R(I)
   25         FORMAT(3X,I5,3X,D24.16)
   20         CONTINUE

              WRITE(6,30)
   30         FORMAT(3X,"Eigenvectors"/,"The Eigenvectors 
     1         Corresponding to the eigenvalue R(I) is in FIHFI8(J,I)")

              DO 40 I=1,N1F
              WRITE(6,35)(FIHFI8(J,I),J=1,N1F)
   35         FORMAT(10D12.4)
   40         CONTINUE

              WRITE(6,41)
   41         FORMAT(3X,"Linear Parameters for Bound State Wavefunction",/,
     1          "Stored in C(I)")

              DO 42 I=1,N1F
              C(I)=FIHFI8(I,1)
              WRITE(6,43)I,C(I)
   42         CONTINUE

   43         FORMAT(3X,"I=",I3,"C(I)=",D15.7)

              PRINT *, C(1)

              RETURN

              END

!Sample Input for executable below:
!0.40000.050
!07
!1
!EOJ
