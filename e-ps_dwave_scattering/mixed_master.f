! THIS FILE COMPUTES THE SINGLET AND TRIPLET DWAVE SCATTERING PHASE
! SHIFTS FOR ELASTIC E-PS SCATTERING. THIS FILE INCLUDES THE MIXED
! SYMMETRY TERM IN THE SHORT-RANGE EXPANSION OF THE TRIAL WAVE FUNCTION
!
! AUTHOR: WILLIAM MITCHELL (UNIVERSITY OF NORTH TEXAS) UNDER DR. SANDRA
! WARD QUINTANILLA

      IMPLICIT REAL * 8 (A-H,O-Z)

      COMMON /CJB4/FILFI(2000,2000),FIFI(2000,2000),
     1FILS(2000,2,10),FILC(2000,2,10)
      COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
      COMMON AK(10),AKAP(10)
      COMMON RK(2000),RL(2000),RM(2000)
      COMMON SP(17),TP(17),R3P(17)
      COMMON XSIN(999),CCSIN(999),XCOL(999),CCOL(999),XSQ(30)
      COMMON CSQ(30)
      COMMON ALPHA2,GAMMA2,ZETA,PI
      COMMON K(2000),L(2000),M(2000)
      COMMON NK,NS,NC,N1F,N2F,N3F,IHDPTF,N1FP1,N1PN2F
      COMMON/EPS2/SPINN,TB,ISPIN,NB

      OPEN(unit=5,file='Pd.data',status='old')
      OPEN(unit=6,file='Pd.txt',status='unknown')

      CALL CONST
      WRITE(6,12)
   12 FORMAT(" GAUSSIAN POINTS") 
      WRITE(6,13)( XSIN(I),CCSIN(I),I,I=1,NS)
      WRITE(6,13)( XCOL(I),CCOL(I),I,I=1,NC)
   13 FORMAT(2D20.10,2X,I3)
      WRITE(6,39)NS,NC
   39 FORMAT(3X,"NS=",I3,2X,"NC=",I3)

      CALL SINGLE
      WRITE(6,14)
   14 FORMAT(" SLS")
      WRITE(6,1) SLS
    1 FORMAT(4D20.8)
      WRITE(6,15)
   15 FORMAT(" SLC")
      WRITE(6,1) SLC
      WRITE(6,16)
   16 FORMAT(" CLS")
      WRITE(6,1) CLS
      WRITE(6,17)
   17 FORMAT(" CLC")
      WRITE(6,1) CLC

      CALL COLUMN
      WRITE(6,20)
   20 FORMAT(18X,"FILS",36X,"FILC")
      DO 21 IK=1,NK
      DO 22 I=1,N3F
      WRITE(6,23) FILS(I,1,IK),FILS(I,2,IK),FILC(I,1,IK),FILC(I,2,IK)
   23 FORMAT(4D25.15)
   22 CONTINUE
   21 CONTINUE

      CALL SQUARE
      READ(5,4) IPRINT
    4 FORMAT(I2)
      WRITE(6,43)IPRINT
   43 FORMAT(3X,"IPRINT=",I5)
      IF(IPRINT.NE.0) GO TO 5
      WRITE(6,30)
   30 FORMAT(" FILFI")
      DO 32 I=1,N3F
      WRITE(6,33) I,(FILFI(I,J),J=1,N3F)
   33 FORMAT(I5/(10D12.4))
   32 CONTINUE
      WRITE(6,31)
   31 FORMAT(" FIFI")
      DO 34 I=1,N3F
      WRITE(6,33) I,(FIFI(I,J),J=1,N3F)
   34 CONTINUE 
    5 CONTINUE 

      CALL RMAT

      STOP
      END

      SUBROUTINE CONST
              ! THIS SUBROUTINE DOES PRELIMINARY NUMERICAL WORK AND
              ! ORDERING FOR THE LATER SUBROUTINES

              IMPLICIT REAL * 8 (A-H,O-Z)

              COMMON /CJB4/FILFI(2000,2000),FIFI(2000,2000),
     1        FILS(2000,2,10),FILC(2000,2,10)
              COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
              COMMON AK(10),AKAP(10)
              COMMON RK(2000),RL(2000),RM(2000)
              COMMON SP(17),TP(17),R3P(17)
              COMMON XSIN(999),CCSIN(999),XCOL(999),CCOL(999),XSQ(30)
              COMMON CSQ(30)
              COMMON ALPHA2,GAMMA2,ZETA,PI
              COMMON K(2000),L(2000),M(2000)
              COMMON NK,NS,NC,N1F,N2F,N3F,IHDPTF,N1FP1,N1PN2F
              COMMON/EPS2/SPINN,TB,ISPIN,NB

              DIMENSION KF(2000),LF(2000),MF(2000)

              N1FP1=0

              DO 28 I=1,2000
              K(I) = 0
              L(I) = 0
              M(I) = 0
              RK(I) = 0.0D0
              RL(I) = 0.0D0
              RM(I) = 0.0D0
              KF(I) = 0
              LF(I) = 0
              MF(I) = 0
   28         CONTINUE

              PI = 4.0D0*DATAN(1.0D0)

    1         FORMAT(3F7.4)
              READ(5,1)ALPHA2,GAMMA2,ZETA
              READ(5,3) NS,NC
    3         FORMAT(2I3)
              WRITE(6,44)NS,NC
   44         FORMAT(3X,"NS=",I3,3X,"NC=",I3)

              N = NS

              CALL LWTS(N,XSIN,CCSIN,1.0D0,0.0D0)
              CALL LWTS(NC,XCOL,CCOL,1.0D0,0.0D0)

              READ(5,2) IHDPTF
    2         FORMAT(5I3)
              WRITE(6,10) ALPHA2,GAMMA2,ZETA
   10         FORMAT(" A2=",F6.3," G2=",F6.3," ZT=",F6.3)
              WRITE(6,11) IHDPTF
   11         FORMAT(" DEGREE OF POLYNOMIAL=",I3)
              READ(5,2)ISPIN

              SPINN = ((-1.0)**ISPIN)

              WRITE(6,88)ISPIN,SPINN
   88         FORMAT(2X,"ISPIN=",I1,3X,"SPINN=",D12.4)
              READ(5,89)NB
   89         FORMAT(I3)

              TB = REAL(NB)

              WRITE(6,91)
              WRITE(6,90)NB,TB
   91         FORMAT(3X,"NB IS THE POWER TO WHICH THE SHIELDING FUNCTION
     1 IS RAISED TO")
   90         FORMAT(3X,"NB=",I3,3X,"TB=",D12.4)

              IHDPP1 = IHDPTF+1

              WRITE(6,15)IHDPP1
   15         FORMAT(3X,"IHDPP1=",I3)

              IND = 0

              DO 65 I=1,IHDPP1
              DO 66 IK=1,I
              ILMAX = I-IK+1
              DO 67 JK=1,ILMAX
              IND = IND+1
              I1 = IK-1
              I2 = JK-1
              I12 = I-1-I1-I2
              KF(IND) = I1
              LF(IND) = I2
              MF(IND) = I12
   67         CONTINUE
   66         CONTINUE
   65         CONTINUE

              NT = IND
              SP(1) = 1.0D0
              TP(1) = 1.0D0
              R3P(1) = 1.0D0

              WRITE(6,72)(KF(I),LF(I),MF(I),I,I=1,NT)
   72         FORMAT(4I5)

              K0UNT = 0 

              DO 50 I=1,NT
              IRES = MOD(LF(I),2)
              IF(IRES.EQ.0) THEN
                      K0UNT = K0UNT+1
              END IF
   50         CONTINUE

              WRITE(6,51)IHDPTF,K0UNT
   51         FORMAT(3X,"FOR W=",I3," THERE ARE",I3," EVEN POWERS OF T")

              N1F = K0UNT
              N2F = NT-N1F
              N1PN2F = NT

              IF (ISPIN.EQ.0) THEN
                      N3F = N1PN2F+N1F
              ELSE 
                      N3F = N1PN2F+N2F
              END IF

              WRITE(6,52)IHDPTF,N1F,N2F,N1PN2F,N3F
   52         FORMAT(3X,"IHDPTF=W=",I3," N1F=",I3," N2F=",I3," N1PN2F=",
     1        I3," N3F=",I3)
              WRITE(6,53)
   53         FORMAT(3X,"CHECK THAT N1F+N2F=N1PN2F:")

              KE = 0
              KO = N1F

              DO 100 I=1,NT
              IRES = MOD(LF(I),2)
              IF(IRES.EQ.0) THEN
                      KE = KE+1
                      K(KE) = KF(I)
                      L(KE) = LF(I)
                      M(KE) = MF(I)
                      RK(KE) = K(KE)
                      RL(KE) = L(KE)
                      RM(KE) = M(KE)
              ELSE
                      KO = KO+1
                      K(KO) = KF(I)
                      L(KO) = LF(I)
                      M(KO) = MF(I)
                      RK(KO) = K(KO)
                      RL(KO) = L(KO)
                      RM(KO) = M(KO)
              END IF
  100         CONTINUE

              WRITE(6,72)(K(I),L(I),M(I),I,I=1,NT)
              READ(5,2)NK
              WRITE(6,83)NK
   83         FORMAT(1X,"NUMBER OF ENERGIES=",I3)

              DO 73 I=1,NK
              READ(5,1)AK(I)
              WRITE(6,84)NK,I,AK(I)
   84         FORMAT(1X,"NK=",I3,3X,"I=",I2,3X,"AK(I)=",D14.7)
   73         CONTINUE

              IF (ISPIN.EQ.0) THEN
                      DO 444 I=N1PN2F+1,N3F
                      K(I) = K(I-N1PN2F)
                      RK(I) = RK(I-N1PN2F)
                      L(I) = L(I-N1PN2F)
                      RL(I) = RL(I-N1PN2F)
                      M(I) = M(I-N1PN2F)
                      RM(I) = RM(I-N1PN2F)
  444                 CONTINUE

              ELSE IF (ISPIN.EQ.1) THEN
                      DO 445 I=N1PN2F+1,N3F
                      K(I) = K(I-N2F)
                      RK(I) = RK(I-N2F)
                      L(I) = L(I-N2F)
                      RL(I) = RL(I-N2F)
                      M(I) = M(I-N2F)
                      RM(I) = RM(I-N2F)
  445                 CONTINUE
              END IF

              RETURN
              END

      SUBROUTINE SINGLE
              ! THIS SUBROUTINE COMPUTES THE LONG-RANGE--LONG-RANGE
              ! MATRIX ELEMENTS
              IMPLICIT REAL * 8 (A-H,O-Z)

              COMMON /CJB4/FILFI(2000,2000),FIFI(2000,2000),
     1        FILS(2000,2,10),FILC(2000,2,10)
              COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
              COMMON AK(10),AKAP(10)
              COMMON RK(2000),RL(2000),RM(2000)
              COMMON SP(17),TP(17),R3P(17)
              COMMON XSIN(999),CCSIN(999),XCOL(999),CCOL(999),XSQ(30)
              COMMON CSQ(30)
              COMMON ALPHA2,GAMMA2,ZETA,PI
              COMMON K(2000),L(2000),M(2000)
              COMMON NK,NS,NC,N1F,N2F,N3F,IHDPTF,N1FP1,N1PN2F
              COMMON/EPS2/SPINN,TB,ISPIN,NB

              COMPLEX*16 CPS

              DIMENSION SLC11(10),CLS11(10),SLS11(10),CLC11(10)
              DIMENSION SLC21(10),SLS21(10),CLS21(10),CLC21(10)

              ALIMIT = 1.0D-27

              DO 23 II=1,2
              DO 24 JJ=1,2
              DO 27 IK=1,NK
              SLS(II,JJ,IK) = 0.0D0
              SLC(II,JJ,IK) = 0.0D0
              CLS(II,JJ,IK) = 0.0D0
              CLC(II,JJ,IK) = 0.0D0
   27         CONTINUE
   24         CONTINUE
   23         CONTINUE

              DO 10 IK=1,NK
              SLS11(IK) = 0.0D0
              CLS11(IK) = 0.0D0
              SLC11(IK) = 1.0D0
              CLC11(IK) = 0.0D0
              SLC21(IK) = 0.0D0
              SLS21(IK) = 0.0D0
              CLS21(IK) = 0.0D0
              CLC21(IK) = 0.0D0
   10         CONTINUE

              ! THE FOLLOWING COMPUTES THE CLC DIRECT-DIRECT ELEMENT
              RECZ = 1.0D0/ZETA
              NN = NS
              CNST1 = -1.0D0

              DO 20 IRHO=1,NN
              RHOO = XSIN(IRHO)*RECZ
              EXZ = EXP(-ZETA*RHOO)
              SHF = (1.0D0-EXP(-ZETA*RHOO))
              SHFN = SHF**TB
              SHFN1 = SHF**(TB-1.0D0)
              SHFN2 = SHF**(TB-2.0D0)
              ZMSHF = (TB*EXZ-1.0D0)
              RH1 = 1.0D0/RHOO
              RH2 = RH1*RH1
              RH3 = RH1**3

              DO 25 IK=1,NK
              AKK = AK(IK)
              AK1 = 1.0D0/AKK
              AK2 = AK1*AK1
              AK3 = AK1**3
              COS1 = COS(AKK*RHOO)
              SIN1 = SIN(AKK*RHOO)
              CY11 = (COS1*AK1*RH1)-(3.0D0*SIN1*AK2*RH2)-
     1               (3.0D0*COS1*AK3*RH3)

              BRA1 = -TB*SHFN1*2.0D0*RHOO*(SIN1+COS1*AK1*
     1               RH1+2.0D0*CY11)
              BRA2 = TB*ZETA*(RHOO**2)*SHFN2*ZMSHF*CY11
              FUNCC = SHFN*AKK*CNST1*CY11*(BRA1+BRA2)
              CLC11(IK) = CLC11(IK)+FUNCC*CCSIN(IRHO)
   25         CONTINUE
   20         CONTINUE

              WRITE(6,100)
  100         FORMAT(3X,"CLC11(IK)")
              DO 130 IK=1,NK
              WRITE(6,120)IK,CLC11(IK)
  120         FORMAT(3X,I2,3X,D12.4)
  130         CONTINUE

              ! THE FOLLOWING COMPUTES ALL OF THE DIRECT-EXCHANGE
              ! ELEMENTS
              CX = 0.5D0
              CY = 0.25D0
              CZ = 0.25D0
              CXINV = 1.0D0/CX
              CYINV = 1.0D0/CY
              CZINV = 1.0D0/CZ
              EXFAC = (5.0D0/(48.0D0*PI))*CXINV*CYINV*CZINV*
     1                (1.0D0/4.0D0)

              DO 1 IX=1,NS
              X = XSIN(IX)*CXINV
              DO 2 IY=1,NS
              Y = XSIN(IY)*CYINV
              DO 3 IZ=1,NS
              Z = XSIN(IZ)*CZINV
              R1 = 0.5D0*(X+Z)
              R2 = 0.5D0*(X+Y)
              R3 = 0.5D0*(Y+Z)
              RHO = 0.5D0*SQRT(2.0D0*(R3*R3+R1*R1)-R2*R2)
              RHOEX = 0.5D0*SQRT(2.0D0*(R2*R2+R3*R3)-R1*R1)

              FAC1D = 1.0D0-EXP(-ZETA*RHO)
              FAC1EX = 1.0D0-EXP(-ZETA*RHOEX)

              FACN = FAC1D**NB
              FACNEX = FAC1EX**NB
              FACN1 = FAC1D**(TB-1)
              FACN2 = FAC1D**(TB-2)
              EX1D = EXP(-ZETA*RHO)
              FACNN = (TB*EX1D-1)

              RECR1 = 1.0D0/R1
              RECR2 = 1.0D0/R2
              RECR3 = 1.0D0/R3

              RECRHO = 1.0D0/RHO
              RRHOEX = 1.0D0/RHOEX
              R12 = RECRHO*RRHOEX
              RREC2 = RECRHO*RECRHO
              RREC3 = RECRHO**3

              RECEX2 = RRHOEX*RRHOEX
              RECEX3 = RRHOEX**3

              POT13 = (RECR3-RECR1)

              COSA=((RHO**2)+(0.25D0*R2**2)-(R3**2))/(RHO*R2)
              RTERM=9.0D0*RECEX2*(((9.0D0*PI*R2**2)/20.0D0)*
     1              ((1.0D0/3.0D0)+(2.0D0*COSA**2)/3.0D0)+
     1              ((RHO*PI)/5.0D0)*(RHO-3.0D0*COSA*R2))-4.0D0*PI

              COEF = R1*R2*R3*CCSIN(IX)*CCSIN(IY)*CCSIN(IZ)

              IF(COEF.LT.ALIMIT)GO TO 2

              DO 4 IK=1,NK
              ARG1 = AK(IK)*RHO
              ARG2 = AK(IK)*RHOEX
              AREC1 = 1.0D0/AK(IK)
              AREC2 = AREC1*AREC1
              AREC3 = AREC1**3
              SIN1 = SIN(ARG1)
              SIN2 = SIN(ARG2)
              COS1 = COS(ARG1)
              COS2 = COS(ARG2)

              SY = (3.0D0*SIN1*RREC3*AREC3)-(SIN1*RECRHO*AREC1)-
     1             (3.0D0*COS1*AREC2*RREC2)              
              SYEX = (3.0D0*SIN2*RECEX3*AREC3)-(SIN2*RRHOEX*AREC1)-
     1               (3.0D0*COS2*AREC2*RECEX2)

              CY1 = (COS1*AREC1*RECRHO)-(3.0D0*SIN1*AREC2*RREC2)-
     1              (3.0D0*COS1*AREC3*RREC3)
              CY = CY1*FACN

              CY1EX =(COS2*RRHOEX*AREC1)-(3.0D0*SIN2*RECEX2*AREC2)-
     1               (3.0D0*COS2*RECEX3*AREC3)
              CYEX = CY1EX*FACNEX

              SYL = POT13*SY

              FINT1 = -TB*ZETA*FACN1*EX1D*(2.0D0/RHO)*(SIN1+COS1*AREC1*
     1                RECRHO+2.0D0*CY1)
              FINT2 = EX1D*TB*(ZETA**2)*FACN2*FACNN*CY1
              FINT3 = -4.0D0/3.0D0*POT13*CY1*FACN

              CYL = FINT1+FINT2+FINT3

              SS21 = RTERM*SYEX*SYL
              SLS21(IK) = SLS21(IK)+SS21*COEF

              CS21 = RTERM*CYEX*(-SYL)
              CLS21(IK) = CLS21(IK)+CS21*COEF

              CC21 = RTERM*CYEX*CYL
              CLC21(IK) = CLC21(IK)+CC21*COEF
    4         CONTINUE
    3         CONTINUE
    2         CONTINUE
    1         CONTINUE

              DO 5 IK=1,NK
              FAC = AK(IK)*EXFAC
              SLS21(IK) = SLS21(IK)*FAC
              CLS21(IK) = CLS21(IK)*FAC
              CLC21(IK) = CLC21(IK)*FAC*(-3.0D0/4.0D0)
              SLC21(IK) = CLS21(IK)
    5         CONTINUE

              DO 7 IK=1,NK
              SLC(1,1,IK) = SLC11(IK)+SPINN*SLC21(IK)
              SLS(1,1,IK) = SLS11(IK)+SPINN*SLS21(IK)
              CLS(1,1,IK) = CLS11(IK)+SPINN*CLS21(IK)
              CLC(1,1,IK) = CLC11(IK)+SPINN*CLC21(IK)
    7         CONTINUE

              RETURN
              END

      SUBROUTINE COLUMN
              ! THIS SUBROUTINE COMPUTES ALL OF THE MIXED-RANGE
              ! ELEMENTS
              IMPLICIT REAL * 8 (A-H,O-Z)

              COMMON /CJB4/FILFI(2000,2000),FIFI(2000,2000),
     1        FILS(2000,2,10),FILC(2000,2,10)
              COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
              COMMON AK(10),AKAP(10)
              COMMON RK(2000),RL(2000),RM(2000)
              COMMON SP(17),TP(17),R3P(17)
              COMMON XSIN(999),CCSIN(999),XCOL(999),CCOL(999),XSQ(30)
              COMMON CSQ(30)
              COMMON ALPHA2,GAMMA2,ZETA,PI
              COMMON K(2000),L(2000),M(2000)
              COMMON NK,NS,NC,N1F,N2F,N3F,IHDPTF,N1FP1,N1PN2F
              COMMON/EPS2/SPINN,TB,ISPIN,NB

              DIMENSION FILS1(15504),FILC1(15504),FILS2(15504)
              DIMENSION FILC2(15504)
              DIMENSION FICOEF(2000)
              EQUIVALENCE(FILS1(1),FILS2(1)),(FILC1(1),FILC2(1))

              ALIMIT = 1.0D-40

              DO 20 I=1,15504
              FILS1(I) = 0.0D0
              FILC1(I) = 0.0D0
   20         CONTINUE

    1         XA = ALPHA2
              XB = ALPHA2+0.5D0
              XC = GAMMA2
    4         CX = 0.5D0*(XA+XB)
              CY = 0.5D0*(XB+XC)
              CZ = 0.5D0*(XC+XA)
              CXINV = 1.0D0/CX
              CYINV = 1.0D0/CY
              CZINV = 1.0D0/CZ
              EXFAC = CXINV*CYINV*CZINV

              DO 5 IX=1,NC
              X = XCOL(IX)*CXINV
              DO 6 IY=1,NC
              Y = XCOL(IY)*CYINV
              DO 7 IZ=1,NC
              Z = XCOL(IZ)*CZINV
              R1 = 0.5D0*(X+Z)
              R2 = 0.5D0*(X+Y)
              R3 = 0.5D0*(Y+Z)
              S = R1+R2
              T = R1-R2
              RECR1 = 1.0D0/R1
              RECR2 = 1.0D0/R2
              RECR3 = 1.0D0/R3
              R1S = R1*R1
              R2S = R2*R2
              R3S = R3*R3
              RHO = 0.5D0*SQRT(2.0D0*(R3*R3+R1*R1)-R2*R2)
              RECRHO = 1.0D0/RHO
              RREC2 = RECRHO*RECRHO
              RREC3 = RECRHO**3
              RHOS = RHO*RHO
              COS12 = 0.5D0*RECR1*RECR2*(R1S+R2S-R3S)

              ANGM = 7.0D0/8.0D0*R1**4-1.0D0/8.0D0*R2**4-0.75D0*R1**2*
     1               R2**2-9.0D0/8.0D0*R3**4+17.0D0/4.0D0*R1**2*R3**2+
     1               0.25D0*R2**2*R3**2

              ANGP = 43.0D0/24.0D0*R1**4+19.0D0/24.0D0*R2**4+91.0D0/
     1               24.0D0*R3**4-23.0D0/12.0D0*R1**2*R2**2+13.0D0/
     1               12.0D0*R1**2*R3**2-35.0D0/12.0D0*R2**2*R3**2

              ANGMIX = 0.4D0*R1**4-0.2D0*R2**4-2.0D0*R3**4+0.2D0*
     1                 R1**2*R2**2-1.6D0*R1**2*R3**2+1.4D0*R2**2*R3**2

              POT13 = (RECR3-RECR1)

              COEF = R1*R2*R3*CCOL(IX)*CCOL(IY)*CCOL(IZ)*EXFAC
              IF(COEF.LT.ALIMIT) GO TO 6

              DO 8 I=2,IHDPTF+1
              SP(I) = SP(I-1)*S
              TP(I) = TP(I-1)*T
              R3P(I) = R3P(I-1)*R3
    8         CONTINUE

              DO 9 I=1,N1PN2F
              FICOEF(I) = COEF*SP(K(I)+1)*TP(L(I)+1)*R3P(M(I)+1)
    9         CONTINUE

              IF (ISPIN.EQ.0) THEN
                      DO 99 I=N1PN2F+1,N3F
                      FICOEF(I) = FICOEF(I-N1PN2F)
   99                 CONTINUE
              ELSE IF (ISPIN.EQ.1) THEN
                      DO 98 I=N1PN2F+1,N3F
                      FICOEF(I) = FICOEF(I-N2F)
   98                 CONTINUE
              END IF

              IF (ISPIN.EQ.0) THEN
                      SING = 2.0D0
                      TRIP = 0.0D0
              ELSE IF (ISPIN.EQ.1) THEN
                      SING = 0.0D0
                      TRIP = 2.0D0
              END IF

              SHFN = (1.0D0-EXP(-ZETA*RHO))
              EXZ = EXP(-ZETA*RHO)
              SHFN1 = SHFN**(TB-1)
              SHFN2 = SHFN**(TB-2)
              SHFNN = SHFN**TB
              SHF = ZETA*(1.0D0-TB*EXZ)
              II = 0

              DO 15 IK=1,NK
              AKK = AK(IK)
              AKNORM = SQRT(AKK)
              AREC1 = 1.0D0/AK(IK)
              AREC2 = AREC1*AREC1
              AREC3 = AREC1**3
              ARG2 = AKK*RHO
              SIN1 = SIN(ARG2)
              COS1 = COS(ARG2)

              SY = (3.0D0*SIN1*RREC3*AREC3)-(SIN1*RECRHO*AREC1)-
     1             (3.0D0*COS1*AREC2*RREC2)

              CONSTS = (PI/32.0D0)*SQRT(AKK/(2.0D0*PI))
              CONSTSM = 0.25D0/32.0D0*SQRT(5.0D0*AKK/3.0D0)
              S2 = RREC2*POT13*SY

              CY1 = (COS1*AREC1*RECRHO)-(3.0D0*SIN1*AREC2*RREC2)-
     1              (3.0D0*COS1*AREC3*RREC3)

              FINT = 0.25D0*EXZ*(-NB*ZETA*SHFN1*RECRHO*(2.0D0*SIN1+
     1               2.0D0*COS1/ARG2+4.0D0*CY1)+NB*ZETA**2*SHFN2*
     1               (NB*EXZ-1.0D0)*CY1)-(1.0D0/3.0D0)*POT13*CY1*SHFNN

              CONSTC = (3.0D0*PI/32.0D0)*SQRT(AKK/(2.0D0*PI))
              CONSTCM = 3.0D0*0.25D0/32.0D0*SQRT(5.0D0*AKK/3.0D0)
              C2 = RREC2*FINT

              DO 16 I=1,N3F
              II = II+1

              IF (I.LE.N1F) THEN
                      ANGS = CONSTS*SING*ANGP + CONSTS*TRIP*ANGM
                      ANGC = CONSTC*SING*ANGP + CONSTC*TRIP*ANGM
              ELSE IF (I.GT.N1F.AND.I.LE.N1PN2F) THEN
                      ANGS = CONSTS*SING*ANGM + CONSTS*TRIP*ANGP
                      ANGC = CONSTC*SING*ANGM + CONSTC*TRIP*ANGP
              ELSE IF (I.GT.N1PN2F) THEN
                      ANGS = CONSTSM*ANGMIX*2.0D0
                      ANGC = CONSTCM*ANGMIX*2.0D0
              END IF

              FICF = FICOEF(I)
              FILS2(II) = FILS2(II)+FICF*S2*ANGS
              FILC2(II) = FILC2(II)+FICF*C2*ANGC
   16         CONTINUE
   15         CONTINUE
    7         CONTINUE
    6         CONTINUE
    5         CONTINUE

              II = 0

              DO 17 IK=1,NK
              DO 18 I=1,N3F
              II = II+1
              FILS(I,1,IK) = FILS1(II)
              FILC(I,1,IK) = FILC1(II)
!-------------FILS1 AND FILS2 ARE IN EQUIVALENCE
              !SAME FOR FILC1 AND FILC2
   18         CONTINUE
   17         CONTINUE

              WRITE(6,111)
  111         FORMAT(3X,"AKNORM")

              RETURN
              END

      SUBROUTINE SQUARE
              ! THIS SUBROUTINE COMPUTES THE SHORT-RANGE--SHORT-RANGE
              ! MATRIX ELEMENTS
              IMPLICIT REAL * 8 (A-H,O-Z)

              COMMON /CJB4/FILFI(2000,2000),FIFI(2000,2000),
     1        FILS(2000,2,10),FILC(2000,2,10)
              COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
              COMMON AK(10),AKAP(10)
              COMMON RK(2000),RL(2000),RM(2000)
              COMMON SP(17),TP(17),R3P(17)
              COMMON XSIN(999),CCSIN(999),XCOL(999),CCOL(999),XSQ(30)
              COMMON CSQ(30)
              COMMON ALPHA2,GAMMA2,ZETA,PI
              COMMON K(2000),L(2000),M(2000)
              COMMON NK,NS,NC,N1F,N2F,N3F,IHDPTF,N1FP1,N1PN2F
              COMMON/EPS2/SPINN,TB,ISPIN,NB

              DIMENSION FI(2000),FIJCF(2000),FLFSW(2000)

              DO 15 I=1,2000
              DO 16 J=1,2000
              FILFI(I,J) = 0.0D0
              FIFI(I,J) = 0.0D0
   16         CONTINUE
   15         CONTINUE

              DO 17 I=1,2000
              FI(I) = 0.0D0
              FIJCF(I) = 0.0D0
              FLFSW(I) = 0.0D0
   17         CONTINUE

              IF(ISPIN.EQ.0)THEN
                      WRITE(6,60)ISPIN
                      I1 = 1
                      I2 = N1F
                      I3 = N1F+1
                      I4 = N1PN2F
                      I5 = 1
                      I6 = N1F+1
                      WRITE(6,61)I1,I2,I3,I4,I5,I6
              ELSE IF (ISPIN.EQ.1)THEN
                      WRITE(6,60)ISPIN
                      I1 = N1F+1
                      I2 = N1PN2F
                      I3 = 1
                      I4 = N1F
                      I5 = N1F+1
                      I6 = 1
                      WRITE(6,61)I1,I2,I3,I4,I5,I6
              END IF

   60         FORMAT(3X,"SUBROUTINE SQUARE ISPIN=",I4)
   61         FORMAT(3X,"I1=",I3," I2=",I3," I3=",I3," I4=",I3," I5=",
     1               I3," I6=",I3)

              NS = IHDPTF+6

              CALL LWTS(NS,XSQ,CSQ,1.0D0,0.0D0)

              WRITE(6,13)( XSQ(I),CSQ(I),I,I=1,NS)
   13         FORMAT(3X,2D20.10,2X,I2)

              XA = ALPHA2+ALPHA2
              XB = XA
              XC = GAMMA2+GAMMA2
              EXFAC = 0.25D0*PI
    4         CX = 0.5D0*(XA+XB)
              CY = 0.5D0*(XB+XC)
              CZ = 0.5D0*(XC+XA)
              CXINV = 1.0D0/CX
              CYINV = 1.0D0/CY
              CZINV = 1.0D0/CZ
              EXFAC = EXFAC*CXINV*CYINV*CZINV

              DO 5 IX=1,NS
              X = XSQ(IX)*CXINV
              DO 6 IY=1,NS
              Y = XSQ(IY)*CYINV
              DO 7 IZ=1,NS
              Z = XSQ(IZ)*CZINV
              R1 = 0.5D0*(X+Z)
              R2 = 0.5D0*(X+Y)
              R3 = 0.5D0*(Y+Z)
              S = R1+R2
              T = R1-R2
              RECR1 = 1.0D0/R1
              RECR2 = 1.0D0/R2
              RECR3 = 1.0D0/R3
              RECS = 1.0D0/S
              R12 = R1*R1
              R22 = R2*R2
              R32 = R3*R3
              COS12 = 0.5D0*RECR1*RECR2*(R12+R22-R32)
              R3FRE = 0.25D0*(2.0D0*(R12+R22)-R32)
              R3294 = 9.0D0/4.0D0*R32
              AR3 = GAMMA2*R3
              ST1 = S*T
              BS = ALPHA2*S
              BT = ALPHA2*T
              ST34 = 0.75D0*ST1

              DO 8 I=2,IHDPTF+1
              SP(I) = SP(I-1)*S
              TP(I) = TP(I-1)*T
              R3P(I) = R3P(I-1)*R3
    8         CONTINUE

              COEF = EXFAC*CSQ(IX)*CSQ(IY)*CSQ(IZ)*R1*R2*R3

              DO 11 IA=1,N1PN2F
              FI(IA) = SP(K(IA)+1)*TP(L(IA)+1)*R3P(M(IA)+1)
   11         CONTINUE

              ST = S*S-T*T
              RECST = 1.0D0/ST
              R3T = R3*R3-T*T
              SR3 = S*S-R3*R3
              RECS2 = RECS*RECS
              RECR32 = RECR3*RECR3
              STR31 = 3.0D0*S*S-T*T-2.0D0*R3*R3
              STR32 = S*S-3.0D0*T*T+2.0D0*R3*R3
              SR38 = 8.0*S*R3

              DO 9 JA=1,N1PN2F
              FIJCF(JA) = SP(K(JA)+1)*TP(L(JA)+1)*R3P(M(JA)+1)*COEF
              TERM27 = RK(JA)*RECS-ALPHA2
              TERM67 = RM(JA)*RECR3-GAMMA2
              TERM11 = ALPHA2*ALPHA2-2.0D0*ALPHA2*RK(JA)*RECS+
     1                 RK(JA)*(RK(JA)-1.0D0)*RECS2
              TERM1 = STR31*R3*TERM11
              TERM2 = SR38*TERM27

              IF(L(JA).EQ.0.OR.L(JA).EQ.1)THEN
                      TERM3 = 0.0D0
              ELSE
                      TERM3 = -4.0D0/3.0D0*RECST*RECR3*STR32*R3*RL(JA)*
     1               (RL(JA)-1.0D0)*SP(K(JA)+1)*TP(L(JA)-1)*R3P(M(JA)+1)
     1                *COEF
              END IF

              TERM4 = -8.0D0*R3*RL(JA)
              TERM55 = GAMMA2*GAMMA2-2.0D0*GAMMA2*RM(JA)*RECR3+RM(JA)*
     1                 (RM(JA)-1.0D0)*RECR32
              TERM5 = ST*R3*TERM55
              TERM6 = 2.0D0*ST*TERM67
              TERM7 = R3T*2.0D0*S*TERM67*TERM27
              TERM8 = SR3*2.0D0*RL(JA)*TERM67
              POTTOL = -(RECR1+RECR2-RECR3)
              TERM = 4.0D0/3.0D0*(-RECST*RECR3*(TERM1+TERM2+TERM4+TERM5+
     1               TERM6+TERM7+TERM8)+POTTOL+0.25D0)
              FLFSW(JA) = TERM*FIJCF(JA)+TERM3
    9         CONTINUE

              IF (ISPIN.EQ.0) THEN
                      DO 111 IA=N1PN2F+1,N3F
                      FI(IA) = FI(IA-N1PN2F)
                      FIJCF(IA) = FIJCF(IA-N1PN2F)
                      FLFSW(IA) = FLFSW(IA-N1PN2F)
  111                 CONTINUE
              ELSE IF (ISPIN.EQ.1) THEN
                      DO 112 IA=N1PN2F+1,N3F
                      FI(IA) = FI(IA-N2F)
                      FIJCF(IA) = FIJCF(IA-N2F)
                      FLFSW(IA) = FLFSW(IA-N2F)
  112                 CONTINUE
              END IF

              ! THIS IS THE +L+ TERM
              DO 100 JA=I1,I2
              CONST = 5.0D0/12.0D0

              ! SF TERM
              PW1 = CONST*FLFSW(JA)*0.75D0*(31.0D0/20.0D0*R1**4+31.0D0/
     1              20.0D0*R2**4+91.0D0/20.0D0*R3**4-2.3D0*R1**2*R2**2-
     1              1.1D0*R1**2*R3**2-1.1D0*R2**2*R3**2)

              ! ALPHA TERM
              PW2 = CONST*(-ALPHA2)*S*(-12.0D0*(R1**2-R1*R2+R2**2)-1.2D0
     1              *R1*R2+6.0D0*COS12*(R1**2-R1*R2+R2**2)+25.2D0*R1*R2*
     1              COS12-20.4D0*R1*R2*COS12**2)
              
              ! K-TERM
              PW3 = CONST*RK(JA)*(-12.0D0*(R1**2-R1*R2+R2**2)-1.2D0
     1              *R1*R2+6.0D0*COS12*(R1**2-R1*R2+R2**2)+25.2D0*R1*R2*
     1              COS12-20.4D0*R1*R2*COS12**2)

              ! L-TERM
              PW4 = CONST*RL(JA)*(-12.0D0*(R1**2+R1*R2+R2**2)+1.2D0
     1              *R1*R2-6.0D0*COS12*(R1**2+R1*R2+R2**2)+25.2D0*R1*R2*
     1              COS12+20.4D0*R1*R2*COS12**2)

              ! M-TERM
              PW5 = CONST*(RM(JA)/R3-GAMMA2)/R3*(-2.7D0*R1**4-2.7D0*
     1              R2**4-17.1D0*R3**4+5.4D0*R1**2*R2**2+1.8D0*R1**2*
     1              R3**2+1.8D0*R2**2*R3**2)

              ! ELECTRON ENERGY TERM
              AKK = AK(1)
              AKS = AKK**2
              PVE = CONST*(-0.75D0)*AKS*(31.0D0/20.0D0*R1**4+31.0D0/
     1              20.0D0*R2**4+91.0D0/20.0D0*R3**4-2.3D0*R1**2*R2**2-
     1              1.1D0*R1**2*R3**2-1.1D0*R2**2*R3**2)

              PT = PW1+FIJCF(JA)*(PW2+PW3+PW4+PW5+PVE)

              DO 150 IA=I5,JA
              FILFI(IA,JA) = FILFI(IA,JA)+PT*FI(IA)
  150         CONTINUE
  100         CONTINUE

              !THIS IS THE -L- TERM
              DO 300 JA=I3,I4
              CONST = 15.0D0/64.00

              ! SF TERM
              PW1 = CONST*FLFSW(JA)*0.75D0*(0.8D0*R1**4+0.8D0*R2**4-
     1              2.4D0*R3**4-1.6D0*R1**2*R2**2+4.8D0*R1**2*R3**2+
     1              4.8D0*R2**2*R3**2)

              ! ALPHA TERM
              PW2 = CONST*(-ALPHA2)*S*(-4.8D0*(R1**2-R1*R2+R2**2)+1.6D0*
     1              R1*R2-8.0D0*R3**2-9.6D0*R3**2*COS12+3.2D0*COS12*
     1              (R1**2-R1*R2+R2**2))

              ! K-TERM
              PW3 = CONST*RK(JA)*(-4.8D0*(R1**2-R1*R2+R2**2)+1.6D0*
     1              R1*R2-8.0D0*R3**2-9.6D0*R3**2*COS12+3.2D0*COS12*
     1              (R1**2-R1*R2+R2**2))

              ! L-TERM
              PW4 = CONST*RL(JA)*(-4.8D0*(R1**2+R1*R2+R2**2)-1.6D0*
     1              R1*R2-8.0D0*R3**2+9.6D0*R3**2*COS12-3.2D0*COS12*
     1              (R1**2+R1*R2+R2**2))

              ! M-TERM
              PW5 = CONST*(RM(JA)/R3-GAMMA2)/R3*(-1.6D0*R1**4-1.6D0*
     1              R2**4+4.8D0*R3**4+3.2D0*R1**2*R2**2-9.6D0*R1**2*
     1              R3**2-9.6D0*R2**2*R3**2)

              ! ELECTRON ENERGY TERM
              AKK = AK(1)
              AKS = AKK**2
              PVE = CONST*(-0.75D0)*AKS*(0.8D0*R1**4+0.8D0*R2**4-
     1              2.4D0*R3**4-1.6D0*R1**2*R2**2+4.8D0*R1**2*R3**2+
     1              4.8D0*R2**2*R3**2)

              PT = PW1+FIJCF(JA)*(PW2+PW3+PW4+PW5+PVE)

              DO 350 IA=I6,JA
              FILFI(IA,JA) = FILFI(IA,JA)+PT*FI(IA)
  350         CONTINUE
  300         CONTINUE

              ! THIS IS THE MIX-L-MIX TERM
              DO 301 JA=N1PN2F+1,N3F
              CONST = 3.0D0/96.0D0*(1.0D0/PI)

              ! SF TERM
              PW1 = CONST*FLFSW(JA)*0.75D0*(-23.0D0/20.0D0*R1**4-23.0D0/
     1              20.0D0*R2**4+73.0D0/20.0D0*R3**4+3.1D0*R1**2*R2**2+
     1              0.7D0*R1**2*R3**2+0.7D0*R2**2*R3**2)

              ! ALPHA TERM
              PW2 = CONST*(-ALPHA2)*S*(-24.0D0/5.0D0*(R1**2-R1*R2+R2**2)
     1              -78.0D0/5.0D0*R1*R2+48.0D0/5.0D0*COS12*(R1**2-R1*R2+
     1              R2**2)+108.0D0/5.0D0*R1*R2*COS12-66.0D0/5.0D0*R1*R2*
     1              COS12**2)

              ! K-TERM
              PW3 = CONST*RK(JA)*(-24.0D0/5.0D0*(R1**2-R1*R2+R2**2)
     1              -78.0D0/5.0D0*R1*R2+48.0D0/5.0D0*COS12*(R1**2-R1*R2+
     1              R2**2)+108.0D0/5.0D0*R1*R2*COS12-66.0D0/5.0D0*R1*R2*
     1              COS12**2)

              ! L-TERM
              PW4 = CONST*RL(JA)*(-24.0D0/5.0D0*(R1**2+R1*R2+R2**2)
     1              +78.0D0/5.0D0*R1*R2-48.0D0/5.0D0*COS12*(R1**2+R1*R2+
     1              R2**2)+108.0D0/5.0D0*R1*R2*COS12+66.0D0/5.0D0*R1*R2*
     1              COS12**2)

              ! M-TERM
              PW5 = CONST*(RM(JA)/R3-GAMMA2)/R3*(2.7D0*R1**4+2.7D0*
     1              R2**4-15.3D0*R3**4-27.0D0/5.0D0*R1**2*R2**2-9.0D0/
     1              5.0D0*R1**2*R3**2-9.0D0/5.0D0*R2**2*R3**2)

              ! ELECTRON ENERGY TERM
              AKK = AK(1)
              AKS = AKK**2
              PVE = CONST*(-0.75D0*AKS)*(-23.0D0/20.0D0*R1**4-23.0D0/
     1              20.0D0*R2**4+73.0D0/20.0D0*R3**4+3.1D0*R1**2*R2**2+
     1              0.7D0*R1**2*R3**2+0.7D0*R2**2*R3**2)

              PT = PW1+FIJCF(JA)*(PW2+PW3+PW4+PW5+PVE)

              DO 351 IA=N1PN2F+1,JA
              FILFI(IA,JA) = FILFI(IA,JA)+PT*FI(IA)
  351         CONTINUE
  301         CONTINUE

              !FIRST IS +L- 
              IF (ISPIN.EQ.0) THEN
                      DO 200 JA=N1F+1,N1PN2F
                      CONST = 5.0D0/16.0D0

                      ! SF TERM
                      PW1 = CONST*FLFSW(JA)*0.75D0*(0.8D0*R1**4-0.8D0*
     1                      R2**4+3.2D0*R1**2*R3**2-3.2D0*R2**2*R3**2)

                      ! ALPHA TERM
                      PW2 = CONST*(-ALPHA2)*T*(-16.0D0*(R1**2+R1*R2+
     1                      R2**2)+1.6D0*R1*R2-8.0D0*COS12*(R1**2+R1*
     1                      R2+R2**2)+33.6D0*R1*R2*COS12+27.2D0*R1*R2*
     1                      COS12**2)

                      ! K-TERM
                      PW3 = CONST*RK(JA)/S*T*(-16.0D0*(R1**2+R1*R2+
     1                      R2**2)+1.6D0*R1*R2-8.0D0*COS12*(R1**2+R1*
     1                      R2+R2**2)+33.6D0*R1*R2*COS12+27.2D0*R1*R2*
     1                      COS12**2)

                      ! L-TERM
                      IF (L(JA).EQ.0) THEN
                              PW4 = 0.0D0
                      ELSE
                              PW4 = SP(K(JA)+1)*TP(L(JA))*R3P(M(JA)+1)
     1                              *COEF*RL(JA)*S*CONST*(-16.0D0*(R1**2
     1                              -R1*R2+R2**2)-1.6D0*R1*R2+8.0D0*
     1                              COS12*(R1**2-R1*R2+R2**2)+33.6D0*R1*
     1                              R2*COS12-27.2D0*R1*R2*COS12**2)
                      END IF

                      ! M-TERM
                      PW5 = CONST*(RM(JA)/R3-GAMMA2)/R3*(-1.6D0*R1**4+
     1                      1.6D0*R2**4-6.4D0*R1**2*R3**2+6.4D0*R2**2*
     1                      R3**2)

                      ! ELECTRON ENERGY TERM
                      AKK = AK(1)
                      AKS = AKK**2
                      PVE = CONST*(-0.75D0)*AKS*(0.8D0*R1**4-0.8D0*
     1                      R2**4+3.2D0*R1**2*R3**2-3.2D0*R2**2*R3**2)

                      PT = PW1+FIJCF(JA)*(PW2+PW3+PW5+PVE)+PW4

                      DO 250 IA=1,N1F
                      FILFI(IA,JA) = FILFI(IA,JA)+PT*FI(IA)
  250                 CONTINUE
  200                 CONTINUE

                      ! THIS IS THE MIX/L/- TERM
                      DO 201 JA=N1F+1,N1PN2F
                      CONST = 3.0D0/32.0D0*SQRT(5.0D0/(6.0D0*PI))

                      ! SF TERM
                      PW1 = CONST*FLFSW(JA)*0.75D0*(0.8D0*R1**4-0.8D0*
     1                      R2**4-4.0D0*R1**2*R3**2+4.0D0*R2**2*R3**2)

                      ! ALPHA TERM
                      PW2 = CONST*(-ALPHA2)*T*(64.0D0/5.0D0*(R1**2+R1*
     1                      R2+R2**2)-28.0D0/5.0D0*R1*R2+32.0D0/5.0D0*
     1                      COS12*(R1**2+R1*R2+R2**2)-192.0D0/5.0D0*R1*
     1                      R2*COS12-116.0D0/5.0D0*R1*R2*COS12**2)

                      ! K-TERM
                      PW3 = CONST*RK(JA)/S*T*(64.0D0/5.0D0*(R1**2+R1*
     1                      R2+R2**2)-28.0D0/5.0D0*R1*R2+32.0D0/5.0D0*
     1                      COS12*(R1**2+R1*R2+R2**2)-192.0D0/5.0D0*R1*
     1                      R2*COS12-116.0D0/5.0D0*R1*R2*COS12**2)

                      ! L-TERM
                      IF (L(JA).EQ.0) THEN
                              PW4 = 0.0D0
                      ELSE
                              PW4 = SP(K(JA)+1)*TP(L(JA))*R3P(M(JA)+1)
     1                              *COEF*RL(JA)*S*CONST*(64.0D0/5.0D0*
     1                              (R1**2-R1*R2+R2**2)+28.0D0/5.0D0*R1*
     1                              R2-32.0D0/5.0D0*COS12*(R1**2-R1*R2+
     1                              R2**2)-192.0D0/5.0D0*R1*R2*COS12+
     1                              116.0D0/5.0D0*R1*R2*COS12**2)
                      END IF

                      ! M-TERM
                      PW5 = CONST*(RM(JA)/R3-GAMMA2)/R3*(-8.0D0/5.0D0*
     1                      R1**4+8.0D0/5.0D0*R2**4+8.0D0*R1**2*R3**2-
     1                      8.0D0*R2**2*R3**2)

                      ! ELECTRON ENERGY TERM
                      AKK = AK(1)
                      AKS = AKK**2
                      PVE = CONST*(-0.75D0)*AKS*(0.8D0*R1**4-0.8D0*
     1                      R2**4-4.0D0*R1**2*R3**2+4.0D0*R2**2*R3**2)

                      PT = PW1+FIJCF(JA)*(PW2+PW3+PW5+PVE)+PW4

                      DO 251 IA=N1PN2F+1,N1PN2F+N1F
                      FILFI(JA,IA) = FILFI(JA,IA)+PT*FI(IA)
  251                 CONTINUE
  201                 CONTINUE

                      ! THIS IS THE MIX/L/+ TERM
                      DO 202 JA=1,N1F
                      CONST = 1.0D0/8.0D0*SQRT(5.0D0/(6.0D0*PI))

                      ! SF TERM
                      PW1 = CONST*FLFSW(JA)*0.75D0*(0.2D0*R1**4+0.2D0*
     1                      R2**4-4.0D0*R3**4+0.4D0*R1**2*R2**2-0.2D0*
     1                      R1**2*R3**2-0.2D0*R2**2*R3**2)

                      ! ALPHA TERM
                      PW2 = CONST*(-ALPHA2)*S*(48.0D0/5.0D0*(R1**2-R1*R2
     1                      +R2**2)+21.0D0/5.0D0*R1*R2-24.0D0/5.0D0*
     1                      COS12*(R1**2-R1*R2+R2**2)-144.0D0/5.0D0*R1*
     1                      R2*COS12+87.0D0/5.0D0*R1*R2*COS12**2)

                      ! K-TERM
                      PW3 = CONST*RK(JA)*(48.0D0/5.0D0*(R1**2-R1*R2
     1                      +R2**2)+21.0D0/5.0D0*R1*R2-24.0D0/5.0D0*
     1                      COS12*(R1**2-R1*R2+R2**2)-144.0D0/5.0D0*R1*
     1                      R2*COS12+87.0D0/5.0D0*R1*R2*COS12**2)

                      ! L-TERM
                      PW4 = CONST*RL(JA)*(48.0D0/5.0D0*(R1**2+R1*R2
     1                      +R2**2)-21.0D0/5.0D0*R1*R2+24.0D0/5.0D0*
     1                      COS12*(R1**2+R1*R2+R2**2)-144.0D0/5.0D0*R1*
     1                      R2*COS12-87.0D0/5.0D0*R1*R2*COS12**2)

                      ! M-TERM
                      PW5 = CONST*(RM(JA)/R3-GAMMA2)/R3*(-2.7D0*R1**4-
     1                      2.7D0*R2**4+15.3D0*R3**4+27.0D0/5.0D0*R1**2*
     1                      R2**2+9.0D0/5.0D0*R1**2*R3**2+9.0D0/5.0D0*
     1                      R2**2*R3**2)

                      ! ELECTRON ENERGY TERM
                      AKK = AK(1)
                      AKS = AKK**2
                      PVE = CONST*(-0.75D0)*AKS*(0.2D0*R1**4+0.2D0*
     1                      R2**4-4.0D0*R3**4+0.4D0*R1**2*R2**2-0.2D0*
     1                      R1**2*R3**2-0.2D0*R2**2*R3**2)

                      PT = PW1+FIJCF(JA)*(PW2+PW3+PW4+PW5+PVE)

                      DO 252 IA=N1PN2F+1,N1PN2F+N1F
                      FILFI(JA,IA) = FILFI(JA,IA)+PT*FI(IA)
  252                 CONTINUE
  202                 CONTINUE

              !THIS IS THE -L+ TERM
              ELSE IF(ISPIN.EQ.1)THEN
                      DO 400 JA=N1F+1,N1PN2F
                      CONST = 5.0D0/16.0D0

                      ! SF TERM
                      PW1 = CONST*FLFSW(JA)*0.75D0*(0.8D0*R1**4-0.8D0*
     1                      R2**4+3.2D0*R1**2*R3**2-3.2D0*R2**2*R3**2)

                      ! ALPHA TERM
                      PW2 = CONST*(-ALPHA2)*T*(-9.6D0*(R1**2+R1*R2+
     1                      R2**2)+4.8D0*R1*R2+4.8D0*COS12*(R1**2+R1*R2+
     1                      R2**2)+4.8D0*R1*R2*COS12-14.4D0*R1*R2*
     1                      COS12**2)

                      ! K-TERM
                      PW3 = CONST*RK(JA)/S*T*(-9.6D0*(R1**2+R1*R2+
     1                      R2**2)+4.8D0*R1*R2+4.8D0*COS12*(R1**2+R1*R2+
     1                      R2**2)+4.8D0*R1*R2*COS12-14.4D0*R1*R2*
     1                      COS12**2)

                      ! L-TERM
                      IF (L(JA).EQ.0) THEN
                              PW4 = 0.0D0
                      ELSE
                              PW4 = CONST*SP(K(JA)+1)*TP(L(JA))*
     1                              R3P(M(JA)+1)*COEF*RL(JA)*S*
     1                              (-9.6D0*(R1**2-R1*R2+R2**2)-4.8D0*R1
     1                              *R2-4.8D0*COS12*(R1**2-R1*R2+R2**2)+
     1                              4.8D0*R1*R2*COS12+14.4D0*R1*R2*
     1                              COS12**2)
                      END IF

                      ! M-TERM
                      PW5 = CONST*(RM(JA)/R3-GAMMA2)/R3*(-14.4D0*R1**4+
     1                      14.4D0*R2**4+28.8D0*R1**3*R2*COS12-28.8D0*
     1                      R1*R2**3*COS12)

                      ! ELECTRON ENERGY TERM
                      AKK = AK(1)
                      AKS = AKK**2
                      PVE = CONST*(-0.75D0)*AKS*(0.8D0*R1**4-0.8D0*
     1                      R2**4+3.2D0*R1**2*R3**2-3.2D0*R2**2*R3**2)

                      PT = PW1+FIJCF(JA)*(PW2+PW3+PW5+PVE)+PW4

                      DO 450 IA=1,N1F
                      FILFI(IA,JA) = FILFI(IA,JA)+PT*FI(IA)
  450                 CONTINUE
  400                 CONTINUE

                      ! THIS IS THE MIX/L/+ TERM
                      DO 203 JA=N1F+1,N1PN2F
                      CONST = 1.0D0/8.0D0*SQRT(5.0D0/(6.0D0*PI))

                      ! SF TERM
                      PW1 = CONST*FLFSW(JA)*0.75D0*(0.2D0*R1**4+0.2D0*
     1                      R2**4-4.0D0*R3**4+0.4D0*R1**2*R2**2-0.2D0*
     1                      R1**2*R3**2-0.2D0*R2**2*R3**2)

                      ! ALPHA TERM
                      PW2 = CONST*(-ALPHA2)*S*(48.0D0/5.0D0*(R1**2-R1*R2
     1                      +R2**2)+21.0D0/5.0D0*R1*R2-24.0D0/5.0D0*
     1                      COS12*(R1**2-R1*R2+R2**2)-144.0D0/5.0D0*R1*
     1                      R2*COS12+87.0D0/5.0D0*R1*R2*COS12**2)

                      ! K-TERM
                      PW3 = CONST*RK(JA)*(48.0D0/5.0D0*(R1**2-R1*R2
     1                      +R2**2)+21.0D0/5.0D0*R1*R2-24.0D0/5.0D0*
     1                      COS12*(R1**2-R1*R2+R2**2)-144.0D0/5.0D0*R1*
     1                      R2*COS12+87.0D0/5.0D0*R1*R2*COS12**2)

                      ! L-TERM
                      PW4 = CONST*RL(JA)*(48.0D0/5.0D0*(R1**2+R1*R2
     1                      +R2**2)-21.0D0/5.0D0*R1*R2+24.0D0/5.0D0*
     1                      COS12*(R1**2+R1*R2+R2**2)-144.0D0/5.0D0*R1*
     1                      R2*COS12-87.0D0/5.0D0*R1*R2*COS12**2)

                      ! M-TERM
                      PW5 = CONST*(RM(JA)/R3-GAMMA2)/R3*(-2.7D0*R1**4-
     1                      2.7D0*R2**4+15.3D0*R3**4+27.0D0/5.0D0*R1**2*
     1                      R2**2+9.0D0/5.0D0*R1**2*R3**2+9.0D0/5.0D0*
     1                      R2**2*R3**2)

                      ! ELECTRON ENERGY TERM
                      AKK = AK(1)
                      AKS = AKK**2
                      PVE = CONST*(-0.75D0)*AKS*(0.2D0*R1**4+0.2D0*
     1                      R2**4-4.0D0*R3**4+0.4D0*R1**2*R2**2-0.2D0*
     1                      R1**2*R3**2-0.2D0*R2**2*R3**2)

                      PT = PW1+FIJCF(JA)*(PW2+PW3+PW4+PW5+PVE)

                      DO 253 IA=N1PN2F+1,N1PN2F+N2F
                      FILFI(JA,IA) = FILFI(JA,IA)+PT*FI(IA)
  253                 CONTINUE
  203                 CONTINUE

                      ! THIS IS THE MIX/L/- TERM
                      DO 204 JA=1,N1F
                      CONST = 3.0D0/32.0D0*SQRT(5.0D0/(6.0D0*PI))

                      ! SF TERM
                      PW1 = CONST*FLFSW(JA)*0.75D0*(0.8D0*R1**4-0.8D0*
     1                      R2**4-4.0D0*R1**2*R3**2+4.0D0*R2**2*R3**2)

                      ! ALPHA TERM
                      PW2 = CONST*(-ALPHA2)*T*(64.0D0/5.0D0*(R1**2+R1*
     1                      R2+R2**2)-28.0D0/5.0D0*R1*R2+32.0D0/5.0D0*
     1                      COS12*(R1**2+R1*R2+R2**2)-192.0D0/5.0D0*R1*
     1                      R2*COS12-116.0D0/5.0D0*R1*R2*COS12**2)

                      ! K-TERM
                      PW3 = CONST*RK(JA)/S*T*(64.0D0/5.0D0*(R1**2+R1*
     1                      R2+R2**2)-28.0D0/5.0D0*R1*R2+32.0D0/5.0D0*
     1                      COS12*(R1**2+R1*R2+R2**2)-192.0D0/5.0D0*R1*
     1                      R2*COS12-116.0D0/5.0D0*R1*R2*COS12**2)

                      ! L-TERM
                      IF (L(JA).EQ.0) THEN
                              PW4 = 0.0D0
                      ELSE
                              PW4 = SP(K(JA)+1)*TP(L(JA))*R3P(M(JA)+1)
     1                              *COEF*RL(JA)*S*CONST*(64.0D0/5.0D0*
     1                              (R1**2-R1*R2+R2**2)+28.0D0/5.0D0*R1*
     1                              R2-32.0D0/5.0D0*COS12*(R1**2-R1*R2+
     1                              R2**2)-192.0D0/5.0D0*R1*R2*COS12+
     1                              116.0D0/5.0D0*R1*R2*COS12**2)
                      END IF

                      ! M-TERM
                      PW5 = CONST*(RM(JA)/R3-GAMMA2)/R3*(-8.0D0/5.0D0*
     1                      R1**4+8.0D0/5.0D0*R2**4+8.0D0*R1**2*R3**2-
     1                      8.0D0*R2**2*R3**2)

                      ! ELECTRON ENERGY TERM
                      AKK = AK(1)
                      AKS = AKK**2
                      PVE = CONST*(-0.75D0)*AKS*(0.8D0*R1**4-0.8D0*
     1                      R2**4-4.0D0*R1**2*R3**2+4.0D0*R2**2*R3**2)

                      PT = PW1+FIJCF(JA)*(PW2+PW3+PW5+PVE)+PW4

                      DO 254 IA=N1PN2F+1,N1PN2F+N2F
                      FILFI(JA,IA) = FILFI(JA,IA)+PT*FI(IA)
  254                 CONTINUE
  204                 CONTINUE
              END IF
    7         CONTINUE
    6         CONTINUE
    5         CONTINUE

              RETURN
              END

      SUBROUTINE RMAT
              ! THIS SUBROUTINE SOLVES THE MATRIX EQUATION AND COMPUTES
              ! OUR PHASE SHIFTS
              IMPLICIT REAL * 8 (A-H,O-Z)

              COMMON /CJB4/FILFI(2000,2000),FIFI(2000,2000),
     1        FILS(2000,2,10),FILC(2000,2,10)
              COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
              COMMON AK(10),AKAP(10)
              COMMON RK(2000),RL(2000),RM(2000)
              COMMON SP(17),TP(17),R3P(17)
              COMMON XSIN(999),CCSIN(999),XCOL(999),CCOL(999),XSQ(30)
              COMMON CSQ(30)
              COMMON ALPHA2,GAMMA2,ZETA,PI
              COMMON K(2000),L(2000),M(2000)
              COMMON NK,NS,NC,N1F,N2F,N3F,IHDPTF,N1FP1,N1PN2F
              COMMON/EPS2/SPINN,TB,ISPIN,NB

              DIMENSION A(2000,2000),X(2000,2),B(2000,2),RBAR(2,2)
              DIMENSION DELR(2,2)
              DIMENSION DELU(2,2),UBAR(2,2)

              DIMENSION CLCTIL(2,2,10), CLSTIL(2,2,10)
              DIMENSION FILCTIL(2000,2,10), FILSTIL(2000,2,10)
              DIMENSION SLSTIL(2,2,10), SLCTIL(2,2,10)
              DIMENSION UMAT(4)

              ALLOCATABLE BL(:,:), AL(:,:), BTR(:,:), IPIV(:)
              ALLOCATABLE BLCOM(:,:), ALCOM(:,:), BTRCOM(:,:)
              COMPLEX*16 BLCOM, ALCOM, BTRCOM

              COMPLEX*16 CLCCOM(2,2,10), CLSCOM(2,2,10)
              COMPLEX*16 FILCCOM(2000,2,10), FILSCOM(2000,2,10)
              COMPLEX*16 SLSCOM(2,2,10), SLCCOM(2,2,10)
              COMPLEX*16 ACOM(2000,2000), XCOM(2000,2), BCOM(2000,2)
              COMPLEX*16 UMATCOM(4), ELLTCOM, ELLVCOM, PSTCOM, PSSCOM
              COMPLEX*16 DETUCOM

              DIMENSION WORK1(1)
              COMPLEX * 16 WORK1COM(1), WORKCOM
              ALLOCATABLE WORK(:)
              ALLOCATABLE WORKCOM(:)

              IF (ISPIN.EQ.0) THEN
                      NT = N1F+N2F+N1F
              ELSE IF (ISPIN.EQ.1) THEN
                      NT = N1F+N2F+N2F
              END IF

              NTP1 = NT+1
              NTP2 = NT+2
              NTP3 = NT+3
              NTP4 = NT+4

              ALLOCATE (AL(NTP1, NTP1))
              ALLOCATE (BL(NTP1, 1))
              ALLOCATE (BTR(1, NTP1))
              ALLOCATE (IPIV(INT(NTP1)))
              ALLOCATE (ALCOM(NTP1, NTP1))
              ALLOCATE (BLCOM(NTP1, 1))
              ALLOCATE (BTRCOM(1, NTP1))

              DO 1 IK=1,NK
              !KOHN DEFINITIONS----------------------
              UMAT = [1.0D0,0.0D0,0.0D0,1.0D0]
              CLCTIL = CLC
              CLSTIL = CLS
              FILCTIL = FILC
              FILSTIL = FILS
              SLSTIL = SLS

              DO 2 I=1,1
              DO 20 J=1,1
              A(I,J) = CLCTIL(I,J,IK)
              B(I,J) = CLSTIL(I,J,IK)
   20         CONTINUE

              DO 3 J=2,NTP1
              A(I,J) = FILCTIL(J-1,I,IK)
              A(J,I) = A(I,J)
              B(J,I) = FILSTIL(J-1,I,IK)
    3         CONTINUE
    2         CONTINUE

              AKS = AK(IK)**2

              DO 40 I=1,NT
              IP1 = I+1
              DO 4 J=I,NT
              JP1 = J+1
              A(IP1,JP1) = FILFI(I,J)
              A(JP1,IP1) = A(IP1,JP1)
    4         CONTINUE
   40         CONTINUE

              !INTRODUCING LAPACK ROUTINE HERE
              DO 999 I=1,NTP1
              DO 998 J=1,NTP1
              AL(I,J) = A(I,J)
  998         CONTINUE
  999         CONTINUE

              DO 997 I=1,1
              DO 996 J=1,NTP1
              BL(J,I) = -B(J,I)
  996         CONTINUE
  997         CONTINUE

              INFO = INT(0)
              NTP1L = INT(NTP1)
              NRHS = INT(1)
              LDA = INT(NTP1)
              LDB = INT(NTP1)

              CALL DGESV(NTP1L, NRHS, AL, LDA, IPIV, BL, LDB, INFO)

              WRITE(6,21)
   21         FORMAT(" LINEAR PARAMETERS")
              WRITE(6,22)(BL(I,1),I=1,NTP1)
   22         FORMAT(1D20.8)

              !COMPUTING TRIAL PHASE SHIFT
              ELLT = BL(1,1)
              PST = (UMAT(2)+UMAT(4)*ELLT)/(UMAT(1)+UMAT(3)*ELLT)

              AKS = AK(IK)**2
              ENERGY = 13.605693122994*AKS*1.5D0
              WRITE(6,213)
  213         FORMAT(" ENERGY | K-VALUE")
              WRITE(6,222) ENERGY, AK(IK)
  222         FORMAT(2D20.8)
              WRITE(6,211)
  211         FORMAT(" TAN OF KOHN TRIAL PHASE SHIFT: ")
              WRITE(6,22) PST

              !COMPUTING STATIONARY PHASE SHIFT
              DETU = UMAT(1)*UMAT(4) - UMAT(2)*UMAT(3)

              DO 995 I=1,1
              DO 994 J=1,NTP1
              BTR(I,J) = B(J,I)
  994         CONTINUE
  995         CONTINUE

              ELLV = 0.0D0

              DO 993 I=1,NTP1
              ELLV = ELLV + BTR(1,I)*BL(I,1)
  993         CONTINUE

              ELLV = ELLV + SLSTIL(1,1,IK)
              ELLV = ELLV*(-1.0D0/DETU)
              PSS = (UMAT(2)+UMAT(4)*ELLV)/(UMAT(1)+UMAT(3)*ELLV)

              WRITE(6,212)
  212         FORMAT(" TAN OF KOHN STATIONARY PHASE SHIFT")
              WRITE(6,22) PSS
              WRITE(6,214)
  214         FORMAT("-------------------------------------------")

              !INVERSE KOHN DEFINITIONS-------------------------------
              UMAT = [0.0D0, 1.0D0, -1.0D0, 0.0D0]
              SLSTIL = CLC
              CLCTIL = SLS
              CLSTIL = -SLC
              SLCTIL = -CLS
              FILSTIL = FILC
              FILCTIL = -FILS

              DO 90 I=1,1
              DO 91 J=1,1
              A(I,J) = CLCTIL(I,J,IK)
              B(I,J) = CLSTIL(I,J,IK)
   91         CONTINUE

              DO 93 J=2,NTP1
              A(I,J) = FILCTIL(J-1,I,IK)
              A(J,I) = A(I,J)
              B(J,I) = FILSTIL(J-1,I,IK)
   93         CONTINUE
   90         CONTINUE

              AKS = AK(IK)**2

              DO 94 I=1,NT
              IP1 = I+1
              DO 95 J=I,NT
              JP1 = J+1
              A(IP1,JP1) = FILFI(I,J)
              A(JP1,IP1) = A(IP1,JP1)
   95         CONTINUE
   94         CONTINUE

              !INTRODUCING LAPACK ROUTINE HERE
              DO 899 I=1,NTP1
              DO 898 J=1,NTP1
              AL(I,J) = A(I,J)
  898         CONTINUE
  899         CONTINUE

              DO 897 I=1,1
              DO 896 J=1,NTP1
              BL(J,I) = -B(J,I)
  896         CONTINUE
  897         CONTINUE

              INFO = INT(0)
              NTP1L = INT(NTP1)
              NRHS = INT(1)
              LDA = INT(NTP1)
              LDB = INT(NTP1)

              CALL DGESV(NTP1L, NRHS, AL, LDA, IPIV, BL, LDB, INFO)

              WRITE(6,21)
              WRITE(6,22)(BL(I,1),I=1,NTP1)

              !COMPUTING TRIAL PHASE SHIFT
              ELLT = BL(1,1)
              PST = (UMAT(2)+UMAT(4)*ELLT)/(UMAT(1)+UMAT(3)*ELLT)
              AKS = AK(IK)**2
              ENERGY = 13.605693122994*AKS*1.5D0

              WRITE(6,213)
              WRITE(6,222) ENERGY, AK(IK)
              WRITE(6,314)
  314         FORMAT(" TAN OF INVERSE KOHN TRIAL PHASE SHIFT")
              WRITE(6,22) PST

              !COMPUTING STATIONARY PHASE SHIFT
              DETU = UMAT(1)*UMAT(4) - UMAT(2)*UMAT(3)

              DO 895 I=1,1
              DO 894 J=1,NTP1
              BTR(I,J) = B(J,I)
  894         CONTINUE
  895         CONTINUE

              ELLV = 0.0D0

              DO 893 I=1,NTP1
              ELLV = ELLV + BTR(1,I)*BL(I,1)
  893         CONTINUE

              ELLV = ELLV + SLSTIL(1,1,IK)
              ELLV = ELLV*(-1.0D0/DETU)
              PSS = (UMAT(2)+UMAT(4)*ELLV)/(UMAT(1)+UMAT(3)*ELLV)

              WRITE(6,313)
  313         FORMAT(" TAN OF INVERSE KOHN STATIONARY PHASE SHIFT")
              WRITE(6,22) PSS
              WRITE(6,214)

              !COMPLEX KOHN DEFINITIONS-----------------------------
              UMATCOM = [CMPLX(0.0D0, -1.0D0), CMPLX(1.0D0, 0.0D0),
     1                   CMPLX(0.0D0, 1.0D0), CMPLX(1.0D0, 0.0D0)]
              CLCCOM = CMPLX(CLC-SLS, CLS+SLC)
              SLSCOM = CMPLX(CLC-SLS, -CLS-SLC)
              CLSCOM = CMPLX(CLC+SLS, SLC-CLS)
              SLCCOM = CMPLX(CLC+SLS, CLS-SLC)
              FILCCOM = CMPLX(FILC, FILS)
              FILSCOM = CMPLX(FILC, -FILS)

              DO 190 I=1,1
              DO 191 J=1,1
              ACOM(I,J) = CLCCOM(I,J,IK)
              BCOM(I,J) = CLSCOM(I,J,IK)
  191         CONTINUE

              DO 193 J=2,NTP1
              ACOM(I,J) = FILCCOM(J-1,I,IK)
              ACOM(J,I) = ACOM(I,J)
              BCOM(J,I) = FILSCOM(J-1,I,IK)
  193         CONTINUE
  190         CONTINUE

              AKS = AK(IK)**2

              DO 194 I=1,NT
              IP1 = I+1
              DO 195 J=I,NT
              JP1 = J+1
              ACOM(IP1,JP1) = FILFI(I,J)
              ACOM(JP1,IP1) = ACOM(IP1,JP1)
  195         CONTINUE
  194         CONTINUE

              !INTRODUCING LAPACK ROUTINE HERE
              DO 799 I=1,NTP1
              DO 798 J=1,NTP1
              ALCOM(I,J) = ACOM(I,J)
  798         CONTINUE
  799         CONTINUE

              DO 797 I=1,1
              DO 796 J=1,NTP1
              BLCOM(J,I) = -BCOM(J,I)
  796         CONTINUE
  797         CONTINUE

              INFO = INT(0)
              NTP1L = INT(NTP1)
              NRHS = INT(1)
              LDA = INT(NTP1)
              LDB = INT(NTP1)

              CALL ZGESV(NTP1L,NRHS,ALCOM,LDA,IPIV,BLCOM,LDB,INFO)

              WRITE(6,21)
              WRITE(6,222)(BLCOM(I,1),I=1,NTP1)

              !COMPUTING TRIAL PHASE SHIFT
              ELLTCOM = BLCOM(1,1)
              PSTCOM = (UMATCOM(2)+UMATCOM(4)*ELLTCOM)/
     1                 (UMATCOM(1)+UMATCOM(3)*ELLTCOM)
              AKS = AK(IK)**2
              ENERGY = 13.605693122994*AKS*1.5D0

              WRITE(6,213)
              WRITE(6,222) ENERGY, AK(IK)
              WRITE(6,414)
  414         FORMAT(" TAN OF COMPLEX KOHN TRIAL PHASE SHIFT")
              WRITE(6,22) REAL(PSTCOM)

              !COMPUTING STATIONARY PHASE SHIFT
              DETUCOM = UMATCOM(1)*UMATCOM(4) - UMATCOM(2)*UMATCOM(3)

              DO 795 I=1,1
              DO 794 J=1,NTP1
              BTRCOM(I,J) = BCOM(J,I)
  794         CONTINUE
  795         CONTINUE

              ELLVCOM = CMPLX(0.0D0, 0.0D0)

              DO 793 I=1,NTP1
              ELLVCOM = ELLVCOM + BTRCOM(1,I)*BLCOM(I,1)
  793         CONTINUE

              ELLVCOM = ELLVCOM + SLSCOM(1,1,IK)
              ELLVCOM = ELLVCOM*(-1.0D0/DETUCOM)

              PSSCOM = (UMATCOM(2)+UMATCOM(4)*ELLVCOM)/
     1                 (UMATCOM(1)+UMATCOM(3)*ELLVCOM)

              WRITE(6,413)
  413         FORMAT(" TAN OF COMPLEX KOHN STATIONARY PHASE SHIFT")
              WRITE(6,22) REAL(PSSCOM)
              WRITE(6,214)

    1         CONTINUE
              RETURN
              END

      SUBROUTINE LWTS(N,Y,W,A,B)
              IMPLICIT REAL * 8 (A-H,O-Z)

              DIMENSION Y(N),W(N)

     0DATA D1, D3, D15, D19, D24, D25, D255/1.0D0, 3.0D0, 
     1         15.0D0, 1.9D0, 2.4D0, 2.5D0, 2.55D0/

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

              IF (ABS(DX).LE.1.0D-12) GO TO 7

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
              IMPLICIT REAL * 8 (A-H,O-Z)

              REAL * 8 L

              DIMENSION L(3)

              L(1)=0.0D0
              L(2)=1.0D0

              DO 1 I=1,N
              L(3)=((2*I-1-X)*L(2)-(I-1)*L(1))/I
              L(1)=L(2)
              L(2)=L(3)
    1         CONTINUE

              F=L(2)
              D=N/X*(F-L(1))

              RETURN
              END
