c*******************************************************************
c PROGRAM FOR FILLING UP THE MISSING DATA USING DOUBLE MASS METHOD *
c (by linear equation Y=A+B.X) 			  						   *
c TO FIL UP BY ROW AND CULUMN CAN BE CHOOSEN                       *
c J.FARAJZADEH & A.F.FARD                                          *
c*******************************************************************
      PROGRAM MASSlin
      CHARACTER*25 LIST(237)
      DIMENSION A(800,800),C(800),X(800),Y(800),E(800),Q(800,800),
     *AFIL(800,800)
      OPEN(1,FILE='rainfall.dat')
      OPEN(2,FILE='rainfall.fil')
      OPEN(3,FILE='REGRESSION.OUT')
      OPEN(4,FILE='DATALENGTH.OUT')
      OPEN(5,FILE='SITE.NAM')
      WRITE(*,112)
  112 FORMAT(///' IF THERE ARE MISSING DATA WITHIN DATA MATRIX YOU HAVE'
     */' TO FILL IT UP BY A NEGATIVE NUMBER'///)
      WRITE(*,100)
  100 FORMAT('  INTER THE NUMBER OF CULUMNS (NCULUMN)')
      READ(*,*)NCULUMN
      WRITE(*,110)
  110 FORMAT('  INTER THE NUMBER OF ROWS(NROW)')
      READ(*,*)NROW
      READ(1,*)((A(I,J),J=1,NCULUMN),I=1,NROW)
C------------------------------------------------------------------------
      READ(5,129)(LIST(I),I=1,NCULUMN)
129   FORMAT(A25)
      DO 300 J=1,NCULUMN
      M=0
      L=0
      DO 315 I=1,NROW
      IF(A(I,J).LT.0)THEN
      L=L+1
      ELSE
      M=M+1
      ENDIF
  315 CONTINUE
C	WRITE(4,225)J,M
C  225 FORMAT(2X,'SITE NO.:',I4,2X,'NUMBER OF DATA:',I4)
      WRITE(4,225)LIST(J),J,M
  225 FORMAT(3X,A25,2X,'SITE NO.:',I4,2X,'NUMBER OF DATA:',I4)
  300 CONTINUE
C--------------------------------------------------------------------
      WRITE(*,111)
111   FORMAT(//2X,'TO FIL WITH CULUMNS RELATIONS ENTER 1'/
     */2X,'TO FIL WITH ROWS RELATIONS ENTER 2')
      READ(*,*)CHOICE
      IF(CHOICE.EQ.1)THEN
      DO 43 I=1,NROW
      DO 43 J=1,NCULUMN
      Q(J,I)=A(I,J)
   43 CONTINUE
      N=NCULUMN
      NP=NROW
      ELSE 
      DO 47 I=1,NROW
      DO 47 J=1,NCULUMN
      Q(I,J)=A(I,J)
   47 CONTINUE
      N=NROW
      NP=NCULUMN
      ENDIF
c*********************************************************
      DO 11 K=1,N
      RMAX=0
      DO 10 I=1,N
      IF(I.EQ.K)GOTO 10
C---------------------------------------------------------
C ELIMINATION OF DATA PARALLEL WITH MISSING DATA POSITIONS
C---------------------------------------------------------
      L=0
      DO 20 J=1,NP
      IF(Q(K,J).LT.0.AND.Q(I,J).LT.0)GOTO 10
C      IF(Q(K,J).LT.0.OR.Q(I,J).GE.0)GOTO 20
      IF(Q(K,J).LT.0)GOTO 20
      L=L+1
      E(L)=Q(K,J)
      C(L)=Q(I,J)
   20 CONTINUE
      IF(L.EQ.0)GOTO 10
c----------------------------
c cumulative data making
c----------------------------
      sume=0
      sumc=0
      DO 25 M=1,L
      sume=sume+e(m)
      sumc=sumc+c(m)
      X(M)=sumc
      Y(M)=sume
   25 CONTINUE
c--------------------------------------------
c to restore the last summation in matrices
c--------------------------------------------
      SUMX=X(L)
      SUMY=Y(L)
c--------------------------------------------
      CALL CORR(L,X,Y,R,AC,BC)
C----------------------------------------------------------
C TO RESTORE THE NO. OF EXPLETIVE AND EXPLETIVED SITES AND 
C RELATED REGRESSION COEFFICIENTS
C----------------------------------------------------------   
      IF(R.GE.RMAX)THEN
      KMAX=K
      IMAX=I
      RMAX=R
      BCF=BC
      ACF=AC
      ENDIF
C----------------------------------------------------------
   10 CONTINUE
C----------------------------------------------------------
C TO WRITE THE NO. OF EXPLETIVE AND EXPLETIVED SITES AND 
C RELATED REGRESSION COEFFICIENTS AND EQUATION
C----------------------------------------------------------
      WRITE(3,200)LIST(KMAX),LIST(IMAX),RMAX,ACF,BCF
  200 FORMAT(2X,'THE SITE NO.',1X,A25,1X,'HAS BEEN FILLED BY SITE NO.',
     *1X,A25,1X,'WITH CORRELATION COEFFICIENT R=',F8.5/3X,'THE EQUATION 
     *FOR CUMMULATIVE CASE IS:'/5X,'Y=',F8.3,'+(',F8.4,'X)'//)
C----------------------------------------------------------
      DO 11 J=1,NP
      IF(Q(KMAX,J).LT.0)THEN
      SUX=Q(IMAX,J)+SUMX
      SUY=ACF+BCF*SUMX
      YCUM=ACF+BCF*SUX
      Q(KMAX,J)=YCUM-SUY
      ENDIF
   11 CONTINUE
      IF(CHOICE.EQ.1)THEN
      DO 500 I=1,N
      DO 500 J=1,NP
      AFIL(J,I)=Q(I,J)
  500 CONTINUE
      ELSE 
      DO 501 I=1,N
      DO 501 J=1,NP
      AFIL(I,J)=Q(I,J)
  501 CONTINUE
      ENDIF
      N=NROW
      NP=NCULUMN
      WRITE(2,310)((AFIL(I,J),J=1,NP),I=1,N)
C 	WRITE(*,310)((AFIL(I,J),J=1,NP),I=1,N)
  310 format(237(f8.2,2x))
      STOP
      END
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C SUBROUTINE FOR CALCULATION OF CORRELATION AND REGRESSION COEFFICIENTS
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      SUBROUTINE CORR(N,X,Y,R2,A,B)
      DIMENSION X(800),Y(800)
      SX=0
      SY=0
      SXY=0
      S2X=0
      S2Y=0
      DO 26 I=1,N
      SX=SX+X(I)
      SY=SY+Y(I)
      SXY=SXY+X(I)*Y(I)
      S2X=S2X+X(I)**2
      S2Y=S2Y+Y(I)**2
   26 CONTINUE
      Z1=S2X-SX**2/N
      Z2=S2Y-SY**2/N
      IF(Z1.EQ.0.OR.Z2.EQ.0)GOTO 10
      B=(SXY-SX*SY/N)/(Z1)
      A=(SY-B*SX)/N
      R2=(A*SY+B*SXY-SY**2/N)/(Z2)
10    RETURN
      END
