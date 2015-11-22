C
C*    Group  Linear Solver subroutines (Code DECCON/SOLCON)
C
      SUBROUTINE BLDECC (A,NROW,NCOL,MCON,M,N,IRANK,COND,D,
     1                                            PIVOT,KRED,AH,V)
C     ------------------------------------------------------------
C
C*  Title
C
C*    Deccon - Constrained Least Squares QR-Decomposition
C
C*  Written by        P. Deuflhard
C*  Purpose           Solution of least squares problems, optionally
C                     with equality constraints.
C*  Method            Constrained Least Squares QR-Decomposition
C                     (see references below)
C*  Category          D9b1. -  Singular, overdetermined or
C                              underdetermined systems of linear 
C                              equations, generalized inverses. 
C                              Constrained Least Squares solution
C*  Keywords          Linear Least Square Problems, constrained, 
C                     QR-decomposition, pseudo inverse.
C*  Version           0.9
C*  Revision          April 1984
C*  Latest Change     January 1991
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C*  Copyright     (c) Konrad-Zuse-Zentrum fuer
C                     Informationstechnik Berlin (ZIB)
C                     Takustrasse 7, D-14195 Berlin-Dahlem
C                     phone : + 49/30/84185-0
C                     fax   : + 49/30/84185-125
C*  Contact           Lutz Weimann
C                     ZIB, Division Scientific Computing, 
C                          Department Numerical Analysis and Modelling
C                     phone : + 49/30/84185-185
C                     fax   : + 49/30/84185-107
C                     e-mail: weimann@zib.de
C
C*    References:
C     ===========
C
C       /1/ P.Deuflhard, V.Apostolescu:
C           An underrelaxed Gauss-Newton method for equality
C           constrained nonlinear least squares problems.
C           Lecture Notes Control Inform. Sci. vol. 7, p.
C           22-32 (1978)
C       /2/ P.Deuflhard, W.Sautter:
C           On rank-deficient pseudoinverses.
C           J. Lin. Alg. Appl. vol. 29, p. 91-111 (1980)
C    
C*    Related Programs:     SOLCON (BLSOLC)
C
C  ---------------------------------------------------------------
C
C* Licence
C    You may use or modify this code for your own non commercial
C    purposes for an unlimited time. 
C    In any case you should not deliver this code without a special 
C    permission of ZIB.
C    In case you intend to use the code commercially, we oblige you
C    to sign an according licence agreement with ZIB.
C
C* Warranty 
C    This code has been tested up to a certain level. Defects and
C    weaknesses, which may be included in the code, do not establish
C    any warranties by ZIB. ZIB does not take over any liabilities
C    which may follow from aquisition or application of this code.
C
C* Software status 
C    This code is under partial care of ZIB and belongs to ZIB software
C    class 2.
C
C     ------------------------------------------------------------
C
C*    Summary:
C     ========
C     Constrained QR-decomposition of (M,N)-system  with
C     computation of pseudoinverse in case of rank-defeciency .
C     First MCON rows belong to equality constraints.
C
C     ------------------------------------------------------------
C
C     INPUT PARAMETERS (* MARKS INOUT PARAMETERS)
C     -----------------------------------------------
C
C
C      * A(NROW,NCOL)  INPUT MATRIX
C                      A(M,N) CONTAINS ACTUAL INPUT
C        NROW          DECLARED NUMBER OF ROWS OF A AND AH
C        NCOL          DECLARED NUMBER OF COLUMNS OF A AND AH
C     (*)MCON          NUMBER OF EQUALITY CONSTRAINTS (MCON<=N)
C                      INTERNALLY REDUCED IF EQUALITY CONSTRAINTS
C                      ARE LINEARLY DEPENDENT
C        M             TREATED NUMBER OF ROWS OF MATRIX A
C        N             TREATED NUMBER OF COLUMNS OF MATRIX A
C     (*)IRANK         PSEUDO-RANK OF MATRIX A
C      * COND          PERMITTED UPPER BOUND OF DABS(D(1)/D(IRANKC))
C                      AND OF DABS(D(IRANKC+1))/D(IRANK))
C                      (SUB-CONDITION NUMBERS OF A)
C        KRED          >=0    HOUSEHOLDER TRIANGULARIZATION
C                             (BUILD UP OF PSEUDO-INVERSE,IF IRANK<N )
C                      < 0    REDUCTION OF PSEUDO-RANK OF MATRIX A
C                             SKIPPING HOUSEHOLDER TRIANGULARIZATION
C                             BUILD-UP OF NEW PSEUDO-INVERSE
C        V(N)          REAL WORK ARRAY
C
C     OUTPUT PARAMETERS
C     -----------------
C
C        A(M,N)        OUTPUT MATRIX UPDATING PRODUCT OF HOUSEHOLDER
C                      TRANSFORMATIONS AND UPPER TRIANGULAR MATRIX
C        MCON          PSEUDO-RANK OF CONSTRAINED PART OF MATRIX A
C        IRANK         PSEUDO-RANK OF TOTAL MATRIX A
C        D(IRANK)      DIAGONAL ELEMENTS OF UPPER TRIANGULAR MATRIX
C        PIVOT(N)      INDEX VECTOR STORING PERMUTATION OF COLUMNS
C                      DUE TO PIVOTING
C        COND          SUB-CONDITION NUMBER OF A
C                      (IN CASE OF RANK REDUCTION: SUB-CONDITION NUMBER
C                      WHICH LED TO RANK REDUCTION)
C        AH(N,N)       UPDATING MATRIX FOR PART OF PSEUDO INVERSE
C
C----------------------------------------------------------------------
C
      INTEGER  IRANK, KRED, MCON, M, N, NROW, NCOL, PIVOT(N)
      INTEGER  I, II, IRK1, I1, J, JD, JJ, K, K1, MH, ISUB
      DOUBLE PRECISION    A(NROW,NCOL), AH(NCOL,NCOL), D(N), V(N)
      DOUBLE PRECISION    COND, ONE , DD, DABS, DSQRT
      DOUBLE PRECISION    H, HMAX, S, T, SMALL, ZERO, EPMACH
C     COMMON /MACHIN/ EPMACH, SMALL
C
      PARAMETER( ZERO=0.D0, ONE=1.D0 )
C
      CALL ZIBCONST(EPMACH,SMALL)
      SMALL = DSQRT(EPMACH*1.D1)
C
      IF(IRANK.GT.N) IRANK=N
      IF(IRANK.GT.M) IRANK=M
C
C---1.0 SPECIAL CASE M=1 AND N=1
C
      IF(M.EQ.1 .AND. N.EQ.1) THEN
         PIVOT(1)=1
         D(1)=A(1,1)
         COND=1.D0
         RETURN
      ENDIF
C
C---1.1 INITIALIZE PIVOT-ARRAY
      IF  (KRED.GE.0)  THEN
         DO 1100 J=1,N
1100        PIVOT(J) = J
C        ENDDO
C
C
C---2. CONSTRAINED HOUSEHOLDER TRIANGULARIZATION
C
         JD = 1
         ISUB = 1
         MH = MCON
         IF (MH.EQ.0) MH=M
         K1 = 1
2000     K = K1
         IF (K.NE.N)  THEN
            K1 = K+1
2100        IF (JD.NE.0)  THEN
               DO  2110 J=K,N
                  S = ZERO
                  DO 2111 I=K,MH
2111                 S = S+A(I,J)*A(I,J)
C                 ENDDO
2110              D(J) = S
C              ENDDO
            ENDIF
C
C---2.1     COLUMN PIVOTING
            H = D(K)
            JJ = K
            DO   2120 J=K1,N
               IF (D(J).GT.H)  THEN
                  H = D(J)
                  JJ = J
               ENDIF
2120        CONTINUE
C           ENDDO
            IF (JD.EQ.1)  HMAX = H * SMALL
            JD = 0
            IF (H.LT.HMAX)  THEN
               JD = 1
               GOTO 2100
            ENDIF
            IF (JJ.NE.K)  THEN
C
C---2.2        COLUMN INTERCHANGE
               I = PIVOT(K)
               PIVOT(K) = PIVOT(JJ)
               PIVOT(JJ) = I
               D(JJ) = D(K)
               DO  2210 I=1,M
                  T = A(I,K)
                  A(I,K) = A(I,JJ)
2210              A(I,JJ) = T
C              ENDDO
            ENDIF
         ENDIF
C
         H = ZERO
         DO  2220 I=K,MH
2220        H = H+A(I,K)*A(I,K)
C        ENDDO
         T = DSQRT(H)
C
C---2.3.0  A PRIORI TEST ON PSEUDO-RANK
C
         IF (ISUB.GT.0) DD = T/COND
         ISUB = 0
         IF (T.LE.DD) THEN
C
C---2.3.1 RANK REDUCTION
C
            IF (K.LE.MCON) THEN
C              CONSTRAINTS ARE LINEARLY DEPENDENT
               MCON = K-1
               K1 = K
               MH = M
               JD = 1
               ISUB = 1
               GOTO 2000
            ENDIF
C
            IRANK = K - 1
            IF (IRANK.EQ.0)  THEN
               GOTO 4000
            ELSE
               GOTO 3000
            ENDIF
         ENDIF
C
         S = A(K,K)
         IF (S.GT.ZERO) T = -T
         D(K) = T
         A(K,K) = S-T
         IF (K.EQ.N)  GOTO 4000
C
         T = ONE/(H-S*T)
         DO  2300 J=K1,N
            S = ZERO
            DO  2310 I=K,MH
2310           S = S+A(I,K)*A(I,J)
C           ENDDO
            S = S*T
            DO  2320 I=K,M
2320           A(I,J) = A(I,J)-A(I,K)*S
C           ENDDO
2300        D(J) = D(J)-A(K,J)*A(K,J)
C        ENDDO
C
         IF (K.EQ.IRANK) GOTO 3000
         IF (K.EQ.MCON) THEN
            MH = M
            JD = 1
            ISUB = 1
         ENDIF
         GOTO 2000
      ENDIF
C
C---3. RANK-DEFICIENT PSEUDO-INVERSE
C
3000  IRK1 = IRANK+1
      DO  3300 J=IRK1,N
         DO  3100 II=1,IRANK
            I = IRK1-II
            S = A(I,J)
            IF (II.NE.1)  THEN
               DO  3110 JJ=I1,IRANK
3110              S = S-A(I,JJ)*V(JJ)
C              ENDDO
            ENDIF
            I1 = I
            V(I) = S/D(I)
3100        AH(I,J) = V(I)
C        ENDDO
         DO  3200 I=IRK1,J
            S = ZERO
            I1 = I-1
            DO  3210 JJ=1,I1
3210           S = S+AH(JJ,I)*V(JJ)
C           ENDDO
            IF (I.NE.J)  THEN
               V(I) = -S/D(I)
               AH(I,J) = -V(I)
            ENDIF
3200     CONTINUE
C        ENDDO
3300     D(J) = DSQRT(S+ONE)
C     ENDDO
C
C---4.  EXIT
C
4000  IF (K.EQ.IRANK) T=D(IRANK)
      IF (T.NE.0.D0) COND=DABS(D(1)/T)
      RETURN
C
C     **********  LAST CARD OF DECCON  **********
C
      END
C
      SUBROUTINE BLSOLC (A,NROW,NCOL,MCON,M,N,X,B,IRANK,D,
     @                   PIVOT,KRED,AH,V)
C
C
C     BEST CONSTRAINED LINEAR LEAST SQUARES SOLUTION OF (M,N)-SYSTEM
C     FIRST MCON ROWS COMPRISE MCON EQUALITY CONSTRAINTS
C
C *********************************************************************
C
C     TO BE USED IN CONNECTION WITH SUBROUTINE DECCON
C
C     RESEARCH CODE FOR GENERAL (M,N)-MATRICES     V 19.01.1984
C
C     INPUT PARAMETERS (* MARKS INOUT PARAMETERS)
C     -----------------------------------------------
C
C        A(M,N)      SEE OUTPUT OF DECCON
C        NROW        SEE OUTPUT OF DECCON
C        NCOL        SEE OUTPUT OF DECCON
C        M           SEE OUTPUT OF DECCON
C        N           SEE OUTPUT OF DECCON
C        MCON        SEE OUTPUT OF DECCON
C        IRANK       SEE OUTPUT OF DECCON
C        D(N)        SEE OUTPUT OF DECCON
C        PIVOT(N)    SEE OUTPUT OF DECCON
C        AH(N,N)     SEE OUTPUT OF DECCON
C        KRED        SEE OUTPUT OF DECCON
C      * B(M)        RIGHT-HAND SIDE OF LINEAR SYSTEM, IF (KRED.GE.0)
C                    RIGHT-HAND SIDE OF UPPER LINEAR SYSTEM,
C                                                      IF (KRED.LT.0)
C        V(N)        REAL WORK ARRAY
C
C     OUTPUT PARAMETERS
C     -----------------
C
C        X(N)        BEST LSQ-SOLUTION OF LINEAR SYSTEM
C        B(M)        RIGHT-HAND OF UPPER TRIGULAR SYSTEM
C                    (TRANSFORMED RIGHT-HAND SIDE OF LINEAR SYSTEM)
C
C
      INTEGER  I, II, I1, IH, IRK1, J, JJ, J1, MH
      INTEGER  IRANK, KRED, M, MCON, N, NROW, NCOL, PIVOT(N)
      DOUBLE PRECISION A(NROW,NCOL), AH(NCOL,NCOL)
      DOUBLE PRECISION B(M), D(N), V(N), X(N), S, ZERO
C
C     COMMON /MACHIN/ EPMACH, SMALL
C
C
      PARAMETER( ZERO=0.D0 )
C
C---1. SOLUTION FOR PSEUDO-RANK ZERO
C
      IF (IRANK.EQ.0)  THEN
         DO 1000 I=1,N
1000        X(I) = ZERO
C        ENDDO
         RETURN
      ENDIF
C
      IF (KRED.GE.0 .AND. (M.NE.1 .OR. N.NE.1) ) THEN
C
C---2. CONSTRAINED HOUSEHOLDER TRANSFORMATIONS OF RIGHT-HAND SIDE
C
         MH = MCON
         IF (MH.EQ.0)  MH = M
         DO  2100 J=1,IRANK
            S = ZERO
            DO  2110 I=J,MH
2110           S = S+A(I,J)*B(I)
C           ENDDO
            S = S/(D(J)*A(J,J))
            DO  2120 I=J,M
2120           B(I) = B(I)+A(I,J)*S
C           ENDDO
            IF (J.EQ.MCON)  MH = M
2100     CONTINUE
C        ENDDO
      ENDIF
C
C---3.1  SOLUTION OF UPPER TRIANGULAR SYSTEM
C
      IRK1 = IRANK+1
      DO  3100 II=1,IRANK
         I = IRK1-II
         I1 = I + 1
         S = B(I)
         IF (I1.LE.IRANK)  THEN
            DO  3111  JJ=I1,IRANK
3111           S = S-A(I,JJ)*V(JJ)
C           ENDDO
         ENDIF
3100     V(I) = S/D(I)
C     ENDDO
      IF (IRK1.LE.N) THEN
C
C---3.2  COMPUTATION OF THE BEST CONSTRAINED LSQ-SOLUTION
C
         DO  3210 J=IRK1,N
            S = ZERO
            J1 = J-1
            DO  3211  I=1,J1
3211           S = S+AH(I,J)*V(I)
C           ENDDO
3210        V(J) = -S/D(J)
C        ENDDO
         DO  3220 JJ=1,N
            J = N-JJ+1
            S = ZERO
            IF (JJ.NE.1) THEN
               DO  3221  I=J1,N
3221              S = S+AH(J,I)*V(I)
C              ENDDO
               IF (J.LE.IRANK) THEN
                  V(J) = V(J)-S
                  GOTO 3220
               ENDIF
            ENDIF
            J1=J
            V(J)=-(V(J)+S)/D(J)
3220     CONTINUE
C        ENDDO
      ENDIF
C
C---4. BACK-PERMUTATION OF SOLUTION COMPONENTS
C
      DO  4000 J=1,N
         IH=PIVOT(J)
4000     X(IH) = V(J)
C     ENDDO
      RETURN
C
C     **********  LAST CARD OF SOLCON  **********
C
      END
