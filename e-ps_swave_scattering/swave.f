      ! THIS CODE WAS ORIGINALLY WRITTEN BY DR. SANDRA J. WARD, WHO
      ! ADAPTED IT FROM DR. JOHN HUMBERTSON'S CODE. THIS VERSION HAS
      ! BEEN MODERNIZED AND DIGITIZED BY
      ! WILLIAM MITCHELL (UNIVERSITY OF NORTH TEXAS)
      
      IMPLICIT REAL * 8 (A-H,O-Z)

      COMMON /CJB4/FILFI(969,969),FIFI(969,969),
     1FILS(969,2,10),FILC(969,2,10)
      COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
      COMMON AK(10),AKAP(10)
      COMMON RK(969),RL(969),RM(969)
      COMMON SP(18),TP(18),R3P(18)
      COMMON XSIN(100),CCSIN(100),XCOL(100),CCOL(100),XSQ(20),CSQ(20)
      COMMON ALPHA2, GAMMA2, ZETA, PI
      COMMON K(969),L(969),M(969)
      COMMON NK,NS,NC,N1F,N2F,IHDPTF,N1FP1,N1PN2F
      COMMON/EPS1/RNEPS(10)
      COMMON/EPS2/ISPIN

      ! THE INPUT FILE IS PS.DATA, OUTPUT IS PS.TXT
      OPEN(UNIT=5,FILE='Ps.data',STATUS='old')
      OPEN(UNIT=6,FILE='Ps.txt',STATUS='unknown')
      
      CALL CONST
      WRITE(6,12)
   12 FORMAT(" GAUSSIAN POINTS")
      print *, XSIN(1)
      WRITE(6,13) (XSIN(I),CCSIN(I),I,I=1,NS)
      WRITE(6,13) (XCOL(I),CCOL(I),I,I=1,NC)
   13 FORMAT(2D20.10,2X,I2)
      WRITE(6,39)NS,NC
   39 FORMAT(3X,"NS=",I2,2X,"NC=",I2)

      CALL SINGLE
      WRITE(6,179)
  179 FORMAT(3X," SUBROUTINE SINGLE COMPLETE")
    3 WRITE(6,14)
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
      WRITE(6,181)
  181 FORMAT(3X,"SUBROUTINE COLUMN COMPLETE")
      WRITE(6,20)
   20 FORMAT(18X,"FILS",36X,"FILC")
      DO 21 IK=1,NK
      DO 22 I=1,N1PN2F
      WRITE(6,23) FILS(I,1,IK),FILS(I,2,IK),FILC(I,1,IK),FILC(I,2,IK)
   23 FORMAT(4D20.10)
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
      DO 32 I=1,N1PN2F
      WRITE(6,33) I,(FILFI(I,J),J=1,N1PN2F)
   33 FORMAT(I5/(10D12.4))
   32 CONTINUE
      WRITE(6,31)
   31 FORMAT(" FIFI")
      DO 34 I=1,N1PN2F
      WRITE(6,33) I,(FIFI(I,J),J=1,N1PN2F)
   34 CONTINUE
    5 CONTINUE

      CALL RMAT

      STOP
      END

      SUBROUTINE CONST
              ! THIS SUBROUTINE HANDLES A LOT OF THE PRELIMINARY
              ! NUMERICAL WORK AND ORDERING FOR LATER SUBROUTINES
              IMPLICIT REAL * 8 (A-H,O-Z)

              COMMON /CJB4/FILFI(969,969),FIFI(969,969),
     1        FILS(969,2,10),FILC(969,2,10)
              COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
              COMMON AK(10),AKAP(10)
              COMMON RK(969),RL(969),RM(969)
              COMMON SP(18),TP(18),R3P(18)
              COMMON XSIN(100),CCSIN(100),XCOL(100),CCOL(100),XSQ(20),
     1               CSQ(20)
              COMMON ALPHA2, GAMMA2, ZETA, PI
              COMMON K(969),L(969),M(969)
              COMMON NK,NS,NC,N1F,N2F,IHDPTF,N1FP1,N1PN2F
              COMMON/EPS1/RNEPS(10)
              COMMON/EPS2/ISPIN

              N1FP1 = 0.0D0
              DO 28 I=1,969
              K(I) = 0.0D0
              L(I) = 0.0D0
              M(I) = 0.0D0
              RK(I) = 0.0D0
              RL(I) = 0.0D0
              RM(I) = 0.0D0
   28         CONTINUE

              PI = 4.0D0*ATAN(1.0D0)

              ! THIS READS IN OUR NONLINEAR PARAMETERS
    1         FORMAT(3F6.3)
              READ(5,1) ALPHA2,GAMMA2,ZETA

              ! THESE ARE OUR INTEGRATION POINT NUMBERS
              READ(5,3) NS,NC
    3         FORMAT(2I3)
              WRITE(6,44)NS,NC
   44         FORMAT(3X,"NS=",I2,3X,"NC=",I2)

              N = NS
              CALL LWTS(N,XSIN,CCSIN,1.0D0,0.0D0)
              CALL LWTS(NC,XCOL,CCOL,1.0D0,0.0D0)

              ! THIS IS OUR OMEGA VALUE
              READ(5,2) IHDPTF
    2         FORMAT(5I3)
              WRITE(6,10) ALPHA2,GAMMA2,ZETA
   10         FORMAT(" A2=",F6.2," G2=",F6.2," ZT=",F6.2)
              WRITE(6,11) IHDPTF
   11         FORMAT(" DEGREE OF POLYNOMIAL=",I3)
              READ(5,2)ISPIN
              WRITE(6,88)ISPIN
   88         FORMAT(2X,"ISPIN=",I1)
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
              K(IND) = I1
              L(IND) = I2
              M(IND) = I12
              RK(IND) = K(IND)
              RL(IND) = L(IND)
              RM(IND) = M(IND)
   67         CONTINUE
   66         CONTINUE
   65         CONTINUE

              NT = IND
              SP(1) = 1.0D0
              TP(1) = 1.0D0
              R3P(1) = 1.0D0
              WRITE(6,72)(K(I),L(I),M(I),I,I=1,NT)
   72         FORMAT(4I5)

              KC = 0

              IF (ISPIN.NE.0) GO TO 90 

              DO 100 I=1,NT
              IRES = MOD(L(I),2)
              IF (IRES.EQ.0) THEN
                      KC = KC+1
                      K(KC) = K(I)
                      L(KC) = L(I)
                      M(KC) = M(I)
                      RK(KC) = RK(I)
                      RL(KC) = RL(I)
                      RM(KC) = RM(I)
              END IF
  100         CONTINUE

              GO TO 140 
   90         CONTINUE

              DO 1000 I=1,NT
              IRES = MOD(L(I),2)
              IF (IRES.NE.0) THEN
                      KC = KC+1
                      K(KC) = K(I)
                      L(KC) = L(I)
                      M(KC) = M(I)
                      RK(KC) = RK(I)
                      RL(KC) = RL(I)
                      RM(KC) = RM(I)
              END IF
 1000         CONTINUE
  140         CONTINUE

              KK = KC
              N1F = KK
              N2F = 0
              N1PN2F = N1F+N2F

              WRITE(6,12)N1F,N2F,N1PN2F
   12         FORMAT(3X,"N1F=",I3,"N2F=",I3,3X,"N1PN2F=",I3)
              WRITE(6,72)(K(KC),L(KC),M(KC),KC,KC=1,KK)

              ! THIS IS OUR NUMBER OF ENERGIES
              READ(5,2)NK
              WRITE(6,83)NK
   83         FORMAT(1X,"NUMBER OF ENERGIES=",I3)

              !THIS IS OUR ENERGIES
              DO 73 I=1,NK
              READ(5,1)AK(I)
              WRITE(6,84)I,AK(I)
   84         FORMAT(1X,"I=",I2,3X,"AK(I)=",D14.6)
   73         CONTINUE

              DO 4 IK=1,NK
              RNEPS(IK) = 2.0D0/(3.0D0*SQRT(AK(IK)))
              WRITE(6,29)IK,NK,RNEPS(IK)
   29         FORMAT(3X,"IK=",I2,3X,"NK=",I2,3X,"RNEPS(IK)=",D12.4)
    4         CONTINUE

              RETURN
              END

      SUBROUTINE SINGLE
              !THIS SUBROUTINE COMPUTES THE LONG-RANGE--LONG-RANGE
              !MATRIX ELEMENTS

              IMPLICIT REAL * 8 (A-H,O-Z)

              COMMON /CJB4/FILFI(969,969),FIFI(969,969),
     1        FILS(969,2,10),FILC(969,2,10)
              COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
              COMMON AK(10),AKAP(10)
              COMMON RK(969),RL(969),RM(969)
              COMMON SP(18),TP(18),R3P(18)
              COMMON XSIN(100),CCSIN(100),XCOL(100),CCOL(100),XSQ(20),
     1               CSQ(20)
              COMMON ALPHA2,GAMMA2,ZETA,PI
              COMMON K(969),L(969),M(969)
              COMMON NK,NS,NC,N1F,N2F,IHDPTF,N1FP1,N1PN2F
              COMMON/EPS1/RNEPS(10)
              COMMON/EPS2/ISPIN
              
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
              SLC11(IK) = 0.5D0
              CLC11(IK) = 0.0D0
              SLC21(IK) = 0.0D0
              SLS21(IK) = 0.0D0
              CLS21(IK) = 0.0D0
              CLC21(IK) = 0.0D0
   10         CONTINUE

              ! THIS CALCULATES CLC DIRECT-DIRECT
              ZETA2 = 2.0D0*ZETA
              RECZ = 1.0D0/ZETA
              NN = NS

              DO 20 IRHO=1,NN
              RHOO = XSIN(IRHO)*RECZ
              EXZ = EXP(-ZETA*RHOO)
              SHF = (1.0D0-EXP(-ZETA*RHOO))
              SHF4 = SHF**4

              DO 25 IK=1,NK
              AKK = AK(IK)
              AKREC = 1.0D0/AKK
              COS1 = COS(AKK*RHOO)
              SIN1 = SIN(AKK*RHOO)
              FIND1 = -3.0D0*ZETA*COS1*EXZ
              FIND2 = ZETA*COS1
              FIND3 = 2.0D0*AKK*SIN1*SHF
              FIND = SHF4*COS1*(FIND1+FIND2+FIND3)
              CLC11(IK) = CLC11(IK)+FIND*CCSIN(IRHO)*
     1                    RECZ*3.0D0*ZETA/2.0D0*AKREC
   25         CONTINUE
   20         CONTINUE
              
              WRITE(6,100)
  100         FORMAT(3X,"CLC11(IK)")

              DO 130 IK=1,NK
              WRITE(6,120)IK,CLC11(IK)
  120         FORMAT(3X,I2,3X,D12.4)
  130         CONTINUE

              ! THIS COMPUTES ALL THE DIRECT-EXCHANGE INTEGRALS
              CX = 0.5D0
              CY = 0.25D0
              CZ = 0.25D0
              CXINV = 1.0D0/CX
              CYINV = 1.0D0/CY
              CZINV = 1.0D0/CZ
              EXFAC = 0.03125D0*CXINV*CYINV*CZINV

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
              
              RECR1 = 1.0D0/R1
              RECR2 = 1.0D0/R2
              RECR3 = 1.0D0/R3
              
              RECRHO = 1.0D0/RHO
              RRHOEX = 1.0D0/RHOEX
              RHOS = RHO*RHO
              RHOEXS = RHOEX*RHOEX
              POT13 = 2.0D0*(RECR1-RECR3)
              
              COEF = R1*R2*R3*CCSIN(IX)*CCSIN(IY)*CCSIN(IZ)

              IF(COEF.LT.ALIMIT)GO TO 2

              EX1D = EXP(-ZETA*RHO)
              EX1EX = EXP(-ZETA*RHOEX)
              EXC2D = EXP(-2.0D0*ZETA*RHO)
              EXC2EX = EXP(-2.0D0*ZETA*RHOEX)

              DO 4 IK=1,NK
              ARG1 = AK(IK)*RHO
              ARG2 = AK(IK)*RHOEX              
              SIN1 = SIN(ARG1)
              SIN2 = SIN(ARG2)
              COS1 = COS(ARG1)
              COS2 = COS(ARG2)
              
              SS21 = SIN2*SIN1*RECRHO*RRHOEX*(-POT13)
              SLS21(IK) = SLS21(IK)+SS21*2.0D0/3.0D0*COEF
              CS21 = SIN1*COS2*RECRHO*RRHOEX*(-POT13)*(FAC1EX**3)
              CLS21(IK) = CLS21(IK)+CS21*2.0D0/3.0D0*COEF
              
              FINT1 = 2.0D0/3.0D0*(-POT13)*COS1*(FAC1D**2)
              FINT2 = 6.0D0*ZETA*ZETA*EXC2D*COS1
              FINT3 = 3.0D0*ZETA*EX1D*(2.0*AK(IK)*SIN1+ZETA*COS1)*FAC1D
              
              FINT4 = FINT1-FINT2+FINT3
              CC21 = COS2*RECRHO*RRHOEX*(FAC1EX**3)*FAC1D*FINT4
              
              CLC21(IK) = CLC21(IK)+CC21*COEF
              
    4         CONTINUE
    3         CONTINUE
    2         CONTINUE
    1         CONTINUE

              DO 5 IK=1,NK
              FAC = EXFAC/AK(IK)
              SLS21(IK) = SLS21(IK)*FAC
              CLS21(IK) = CLS21(IK)*FAC
              CLC21(IK) = CLC21(IK)*FAC
              SLC21(IK) = CLS21(IK)
    5         CONTINUE

              SPINN = ((-1.0D0)**ISPIN)

              DO 7 IK=1,NK
              SLC(1,1,IK) = 2.0D0*SLC11(IK)+SPINN*2.0D0*SLC21(IK)
              SLS(1,1,IK) = 2.0D0*SLS11(IK)+SPINN*2.0D0*SLS21(IK)
              CLS(1,1,IK) = 2.0D0*CLS11(IK)+SPINN*2.0D0*CLS21(IK)
              CLC(1,1,IK) = 2.0D0*CLC11(IK)+SPINN*2.0D0*CLC21(IK)
    7         CONTINUE

              ! THIS IS AN OLD TEST AND IS NO LONGER DIRECTLY USEFUL
!              PRINT *, "KOHN LONG-RANGE TERMS: -(CLS)/(CLC)=TAN(DELTA)"
!              PRINT *, -CLS(1,1,1)/CLC(1,1,1)
!              PRINT *, "INVERSE KOHN LONG-RANGE TERMS: -(SLC)/(SLS)=
!     1 COT(DELTA)"
!              PRINT *, -SLC(1,1,1)/SLS(1,1,1)
!
!              PRINT *, "COMPLEX KOHN LONG-RANGE: "
!              CPS = -CMPLX(CLC(1,1,1)+SLS(1,1,1),
!     1 SLC(1,1,1)-CLS(1,1,1))/CMPLX(CLC(1,1,1)-SLS(1,1,1),
!     1 CLS(1,1,1)+SLC(1,1,1))
!
!              CPS = CMPLX(0.0D0, 1.0D0)*(1+CPS)/(1-CPS)
!
!              PRINT *, REAL(CPS)

              RETURN
              END

      SUBROUTINE COLUMN
              ! THIS SUBROUTINE CALCULATES MIXED-RANGE MATRIX ELEMENTS
              IMPLICIT REAL * 8 (A-H,O-Z)

              COMMON /CJB4/FILFI(969,969),FIFI(969,969),
     1        FILS(969,2,10),FILC(969,2,10)
              COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
              COMMON AK(10),AKAP(10)
              COMMON RK(969),RL(969),RM(969)
              COMMON SP(18),TP(18),R3P(18)
              COMMON XSIN(100),CCSIN(100),XCOL(100),CCOL(100),XSQ(20),
     1               CSQ(20)
              COMMON ALPHA2,GAMMA2,ZETA,PI
              COMMON K(969),L(969),M(969)
              COMMON NK,NS,NC,N1F,N2F,IHDPTF,N1FP1,N1PN2F
              COMMON/EPS1/RNEPS(10)
              COMMON/EPS2/ISPIN

              DIMENSION FILS1(15504),FILC1(15504)
              DIMENSION FILS2(15504),FILC2(15504)
              DIMENSION FICOEF(969)

              EQUIVALENCE(FILS1(1),FILS2(1)),(FILC1(1),FILC2(1))

              ALIMIT = 1.0D-40

              DO 20 I=1,15504
              FILS1(I) = 0.0D0
              FILC1(I) = 0.0D0
   20         CONTINUE

    1         XA = ALPHA2
              XB = ALPHA2+0.5D0
              XC = GAMMA2
              EXFAC = 0.25D0*SQRT(PI*0.5D0)
    4         CX = 0.5D0*(XA+XB)
              CY = 0.5D0*(XB+XC)
              CZ = 0.5D0*(XC+XA)
              CXINV = 1.0D0/CX
              CYINV = 1.0D0/CY
              CZINV = 1.0D0/CZ
              EXFAC = EXFAC*CXINV*CYINV*CZINV

              ! THIS IS THE START OF THE GAUSS-LAGUERRE INTEGRATION
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
              RHO = 0.5D0*SQRT(2.0D0*(R3*R3+R1*R1)-R2*R2)
              RECRHO = 1.0D0/RHO
              RHOS = RHO*RHO
              POT13 = 2.0D0*(RECR1-RECR3)

              COEF = R1*R2*R3*CCOL(IX)*CCOL(IY)*CCOL(IZ)*EXFAC

              IF(COEF.LT.ALIMIT) GO TO 6

              DO 8 I=2,18
              SP(I) = SP(I-1)*S
              TP(I) = TP(I-1)*T
              R3P(I) = R3P(I-1)*R3
    8         CONTINUE

              DO 9 I=1,N1F
              FICOEF(I) = COEF*SP(K(I)+1)*TP(L(I)+1)*R3P(M(I)+1)
    9         CONTINUE

              SHFN = (1-EXP(-ZETA*RHO))
              EX2Z = EXP(-2.0D0*ZETA*RHO)
              EXZ = EXP(-ZETA*RHO) 
              II = 0

              DO 15 IK=1,NK
              AKK = AK(IK)
              RNEPSS = RNEPS(IK)
              ARG2 = AKK*RHO
              AKS = AKK*AKK
 
              S2 = -(SIN(ARG2)*RECRHO*POT13*RNEPSS)
              
              SINRHO = SIN(ARG2)*RECRHO
              COSRHO = COS(ARG2)*RECRHO
        
              C2 = RNEPSS*(-POT13*COSRHO*(SHFN**3)-4.5D0*ZETA*SHFN*
     1        (2.0D0*ZETA*EX2Z*COSRHO-EXZ*SHFN*(2.0D0*AKK*SINRHO+ZETA
     2        *COSRHO))) 
              
              DO 16 I=1,N1F
              II = II+1
              FICF = FICOEF(I)
              FILS2(II) = FILS2(II)+FICF*S2
              FILC2(II) = FILC2(II)+FICF*C2
   16         CONTINUE
   15         CONTINUE
    7         CONTINUE
    6         CONTINUE
    5         CONTINUE

              II = 0

              DO 17 IK=1,NK
              DO 18 I=1,N1F
              II = II+1
              FILS(I,1,IK) = FILS1(II)
              FILC(I,1,IK) = FILC1(II)
!-----FILS1 AND FILS2 ARE IN EQUIVALENCE
   18         CONTINUE
   17         CONTINUE

              RETURN
              END
              
      SUBROUTINE SQUARE
              ! THIS SUBROUTINE COMPUTES THE SHORT-RANGE--SHORT-RANGE
              ! MATRIX ELEMENTS

              IMPLICIT REAL * 8 (A-H,O-Z)

              COMMON /CJB4/FILFI(969,969),FIFI(969,969),
     1        FILS(969,2,10),FILC(969,2,10)
              COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
              COMMON AK(10),AKAP(10)
              COMMON RK(969),RL(969),RM(969)
              COMMON SP(18),TP(18),R3P(18)
              COMMON XSIN(100),CCSIN(100),XCOL(100),CCOL(100),XSQ(20),
     1               CSQ(20)
              COMMON ALPHA2,GAMMA2,ZETA,PI
              COMMON K(969),L(969),M(969)
              COMMON NK,NS,NC,N1F,N2F,IHDPTF,N1FP1,N1PN2F
              COMMON/EPS1/RNEPS(10)
              COMMON/EPS2/ISPIN

              DIMENSION FI(969)

              DO 77 JA=1,N1F
              DO 67 IA=1,JA              
              FILFI(IA,JA) = 0.0D0
              FIFI(IA,JA) = 0.0D0
   67         CONTINUE
   77         CONTINUE

              WRITE(6,161)
  161         FORMAT(3X,"BEGINNING OF SUBROUTINE SQU")
              NS = IHDPTF+2

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

              DO 8 I=2,18
              SP(I) = SP(I-1)*S
              TP(I) = TP(I-1)*T
              R3P(I) = R3P(I-1)*R3
    8         CONTINUE

              COEF = EXFAC*CSQ(IX)*CSQ(IY)*CSQ(IZ)*R1*R2*R3

              DO 11 IA=1,N1F
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
              SR38 = 8.0D0*S*R3

              DO 9 JA=1,N1F
              FIJCF = SP(K(JA)+1)*TP(L(JA)+1)*R3P(M(JA)+1)*COEF
              TERM27 = RK(JA)*RECS-ALPHA2
              TERM67 = RM(JA)*RECR3-GAMMA2
              TERM11 = ALPHA2*ALPHA2-2.0D0*ALPHA2*RK(JA)*RECS+K(JA)*
     1                 (K(JA)-1.0)*RECS2
              TERM1 = STR31*R3*TERM11
              TERM2 = SR38*TERM27

              IF(L(JA).EQ.0.OR.L(JA).EQ.1) THEN
                      TERM3 = 0.0D0
              ELSE
                      TERM3 = -4.0D0/3.0D0*RECST*RECR3*STR32*R3*RL(JA)*
     1                        (RL(JA)-1.0)*SP(K(JA)+1)*TP(L(JA)-1)*
     1                        R3P(M(JA)+1)*COEF
              END IF

              TERM4 = -8.0D0*R3*RL(JA)
              TERM55 = GAMMA2*GAMMA2-2.0D0*GAMMA2*RM(JA)*RECR3+RM(JA)*
     1                 (RM(JA)-1.0)*RECR32
              TERM5 = ST*R3*TERM55
              TERM6 = 2.0D0*ST*TERM67
              TERM7 = R3T*2.0D0*S*TERM67*TERM27
              TERM8 = SR3*2.0D0*L(JA)*TERM67
              POTTOL = -(RECR1+RECR2-RECR3)
              TERM = 4.0D0/3.0D0*(-RECST*RECR3*(TERM1+TERM2+TERM4+
     1               TERM5+TERM6+TERM7+TERM8)+POTTOL+0.25)

              DO 10 IA=1,JA
              FILFI(IA,JA) = FILFI(IA,JA)+FI(IA)*(TERM*FIJCF+TERM3)
              FIFI(IA,JA) = FIFI(IA,JA)+FI(IA)*FIJCF
   10         CONTINUE
    9         CONTINUE
    7         CONTINUE
    6         CONTINUE
    5         CONTINUE

              RETURN
              END

      SUBROUTINE RMAT
              !THIS SUBROUTINE SOLVES THE MATRIX EQUATION AND PRODUCES
              !OUR PHASE SHIFTS
              IMPLICIT REAL * 8 (A-H,O-Z)

              COMMON /CJB4/FILFI(969,969),FIFI(969,969),
     1        FILS(969,2,10),FILC(969,2,10)
              COMMON SLS(2,2,10),SLC(2,2,10),CLS(2,2,10),CLC(2,2,10)
              COMMON AK(10),AKAP(10)
              COMMON RK(969),RL(969),RM(969)
              COMMON SP(18),TP(18),R3P(18)
              COMMON XSIN(100),CCSIN(100),XCOL(100),CCOL(100),XSQ(20),
     1               CSQ(20)
              COMMON ALPHA2,GAMMA2,ZETA,PI
              COMMON K(969),L(969),M(969)
              COMMON NK,NS,NC,N1F,N2F,IHDPTF,N1FP1,N1PN2F
              COMMON/EPS1/RNEPS(10)
              COMMON/EPS2/ISPIN

              DIMENSION A(1000,1000),X(1000,2),B(1000,2),RBAR(2,2)
              DIMENSION DELR(2,2)
              DIMENSION DELU(2,2),UBAR(2,2)

              DIMENSION CLCTIL(2,2,10), CLSTIL(2,2,10)
              DIMENSION FILCTIL(969,2,10), FILSTIL(969,2,10)
              DIMENSION SLSTIL(2,2,10), SLCTIL(2,2,10)
              DIMENSION UMAT(4)

              ALLOCATABLE BL(:,:), AL(:,:), BTR(:,:), IPIV(:)
              ALLOCATABLE BLCOM(:,:), ALCOM(:,:), BTRCOM(:,:)

              COMPLEX*16 BLCOM, ALCOM, BTRCOM
              COMPLEX*16 CLCCOM(2,2,10), CLSCOM(2,2,10)
              COMPLEX*16 FILCCOM(969,2,10), FILSCOM(969,2,10)
              COMPLEX*16 SLSCOM(2,2,10), SLCCOM(2,2,10)
              COMPLEX*16 ACOM(1000,1000), XCOM(1000,2), BCOM(1000,2)
              COMPLEX*16 UMATCOM(4), ELLTCOM, ELLVCOM, PSTCOM, PSSCOM
              COMPLEX*16 DETUCOM

              NDIM = 24
              NT = N1F+N2F
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
              !KOHN DEFINITIONS------------------------------
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
              A(IP1,JP1) = FILFI(I,J)-AKS*FIFI(I,J)
              A(JP1,IP1) = A(IP1,JP1)
    4         CONTINUE
   40         CONTINUE

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
   21         FORMAT("  LINEAR PARAMETERS")
              WRITE(6,22)(BL(I,1),I=1,NTP1)
   22         FORMAT(1D20.8)

              !COMPUTING TRIAL PHASE SHIFT
              ELLT = BL(1,1)
              PST = (UMAT(2)+UMAT(4)*ELLT)/(UMAT(1)+UMAT(3)*ELLT)

              AKS = AK(IK)**2
              ENERGY = 13.605693122994*AKS*1.5D0
              WRITE(6,11) ENERGY,AK(IK)
   11         FORMAT(" ENERGY | K-VALUE",2D12.4)
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
              PRINT *, "KOHN PS: ", PSS

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

              DO 32 I=1,1
              DO 50 J=1,1
              A(I,J) = CLCTIL(I,J,IK)
              B(I,J) = CLSTIL(I,J,IK)
   50         CONTINUE

              DO 33 J=2,NTP1
              A(I,J) = FILCTIL(J-1,I,IK)
              A(J,I) = A(I,J)
              B(J,I) = FILSTIL(J-1,I,IK)
   33         CONTINUE
   32         CONTINUE

              AKS = AK(IK)**2

              DO 94 I=1,NT
              IP1 = I+1
              DO 95 J=I,NT
              JP1 = J+1
              A(IP1,JP1) = FILFI(I,J)-AKS*FIFI(I,J)
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
              WRITE(6,11) ENERGY, AK(IK)
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
              PRINT *, "INVERSE KOHN PS: ", PSS

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

              !---TESTING---!
              ! COMPUTING TRIAL/STATIONARY PHASE SHIFTS USING ONLY
              ! LONG-RANGE TERMS
!
!              PRINT *, "LR SWAVE KOHN TRIAL: ", -CLS(1,1,1)/
!     1                                           CLC(1,1,1)
!              AHOLD = -CLS(1,1,1)/CLC(1,1,1)
!              PRINT *, "LR SWAVE KOHN STAT: ",-CLS(1,1,1)*AHOLD-
!     1                                         SLS(1,1,1)
!
!              PRINT *, "LR SWAVE CK TRIAL: ", REAL(-CLSCOM(1,1,1)/
!     1                                              CLCCOM(1,1,1))
!              ACOM(1,1) = -CLSCOM(1,1,1)/CLCCOM(1,1,1)
!              PRINT *,"LR SWAVE CK STAT: ", REAL(-CLSCOM(1,1,1)*
!     1                                      ACOM(1,1)-SLSCOM(1,1,1))
              !---TESTING---!

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
              ACOM(IP1,JP1) = FILFI(I,J)-AKS*FIFI(I,J)
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
  222         FORMAT(2D20.8)

              !COMPUTING TRIAL PHASE SHIFT
              ELLTCOM = BLCOM(1,1)
              PSTCOM = (UMATCOM(2)+UMATCOM(4)*ELLTCOM)/
     1                 (UMATCOM(1)+UMATCOM(3)*ELLTCOM)

              AKS = AK(IK)**2
              ENERGY = 13.605693122994*AKS*1.5D0
              WRITE(6,11) ENERGY, AK(IK)
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
              PRINT *, "COMPLEX KOHN PS: ", REAL(PSSCOM)

              WRITE(6,214)

    1         CONTINUE

              RETURN
              END

      SUBROUTINE LWTS(N,Y,W,A,B)
              IMPLICIT REAL * 8 (a-h,o-z)
 
              DIMENSION Y(N), W(N)
              
     0DATA D1,D3,D15,D19,D24,D25,D255/1.0D0,3.0D0,15.0D0,1.9D0,
     1              2.4D0,2.5D0,2.55D0/
              
              FN = N
              C1 = D1/D19
              C2 = D255/D19

              DO 1 I=1,N
              IF(I.GT.2) GO TO 5
              GO TO (3,4),I
    3         X = D3/(D1+D24*FN)

              GO TO 6

    4         X = D15/(D25*FN+D1)+Y(I-1)

              GO TO 6

    5         X = (C1/I+D1+C2)*Y(I-1)-(C1/I+C2)*Y(I-2)

    6         CALL LN(F,D,X,N)

              DX = F/D
              X = X-DX

              IF (ABS(DX).LE.1.0D-12) GO TO 7
              GO TO 6

    7         W(I) = D1/(X*D*D)
              Y(I) = X
    1         CONTINUE
              
              CAB = EXP(-A*B)/A

              DO 8 I=1,N
              Y(I) = Y(I)/A+B
              W(I) = CAB*W(I)
    8         CONTINUE
              
              RETURN
              END

      SUBROUTINE LN(F,D,X,N)
              IMPLICIT REAL * 8 (a-h,o-z)
              
              REAL * 8 L
              
              DIMENSION L(3)
              
              L(1) = 0.0D0
              L(2) = 1.0D0

              DO 1 I=1,N
              L(3) = ((2*I-1-X)*L(2)-(I-1)*L(1))/I
              L(1) = L(2)
              L(2) = L(3)
    1         CONTINUE

              F = L(2)
              D = N/X*(F-L(1))
              
              RETURN
              END
