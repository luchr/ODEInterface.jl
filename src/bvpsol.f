      SUBROUTINE BVPSOL(FCN,BC,IVPSOL,N,M,T,X,EPS,IOPT,INFO,IRW,RW,
     $                  IIW,IW)
C*    Begin Prologue BVPSOL
      IMPLICIT DOUBLEPRECISION(S)
      EXTERNAL FCN,BC,IVPSOL
      INTEGER N,M
      DOUBLE PRECISION T(M)
      DOUBLE PRECISION X(N,M)
      DOUBLE PRECISION EPS
      INTEGER INFO,IOPT(5)
      INTEGER IRW
      DOUBLE PRECISION RW(IRW)
      INTEGER IIW
      INTEGER IW(IIW)
C
C     ------------------------------------------------------------
C
C*  Title
C
C     (B)oundary (V)alue (P)roblem (Sol)lver for highly nonlinear
C     two point boundary value problems using a local linear
C     solver (condensing algorithm) or a global sparse linear solver
C     for the solution of the arising 
C     linear subproblems.
C
C*  Written by        P. Deuflhard, G.Bader, L. Weimann
C*  Purpose           Solution of nonlinear two-point boundary value
C                     problems.
C*  Method            Local and Global Nonlinear two-point Boundary Value
C                     Problems solver (Multiple shooting approach)
C*  Category          I1b2a - Differential and integral equations
C                             Two point boundary value problems
C*  Keywords          Nonlinear boundary value problems, Multiple
C                     shooting, Newton methods
C*  Version           1.2
C*  Revision          February 2002
C*  Latest Change     January 2004
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
C
C     /1/ R.Bulirsch:
C         Die Mehrzielmethode zur numerischen Loesung von
C         nichtlinearen Randwertproblemen und Aufgaben der
C         optimalen Steuerung.
C         Carl-Cranz-Gesellschaft: Tech.Rep. (Oct.1971)
C
C     /2/ J.Stoer, R.Bulirsch:
C         Einfuehrung in die Numerische Mathematik II.
C         Berlin, Heidelberg, New York: Springer (1st Ed. 1973)
C
C     /3/ P.Deuflhard:
C         A Modified Newton Method for the Solution of
C         Ill-Conditioned Systems of Nonlinear Equations with
C         Application to Multiple Shooting.
C         Numer. Math. 22, 289-315 (1974)
C
C     /4/ P.Deuflhard:
C         Recent Advances in Multiple Shooting Techniques.
C         (Survey Article including further References)
C         In: I.Gladwell, D.K.Sayers (Ed.): Computational
C         Techniques for Ordinary Differential Equations.
C         Section 10, P.217-272.
C         London, New York: Academic Press (1980)
C
C     /5/ P.Deuflhard, G.Bader:
C         Multiple Shooting Techniques Revisited.
C         Univ. Heidelberg, SFB 123, Tech. Rep. 163 (1982)
C
C     /6/ P. Deuflhard:
C         Newton Methods for Nonlinear Problems. -
C         Affine Invariance and Adaptive Algorithms.
C         Series Computational Mathematics 35, Springer (2004)
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
C     External subroutines (to be supplied by the user)
C     =================================================
C
C       FCN(N,T,Y,DY)         Right-hand side of system of
C                             first-order differential equations
C         N                   Input: Number of first order ODE's
C         T                   Input: Actual position in the
C                             interval ³ A ,  B ü
C         Y(N)                Input: Values at T
C         DY(N)               Output: Derivatives at T
C
C       BC(YA,YB,R)           Two-point boundary conditions at ( A
C                             = T(1),  B = T(M))
C         YA(N)               Input: Values at A = T(1)
C         YB(N)               Input: Values at B = T(M)
C         R(N)                Output: Values of
C                             boundary conditions function
C
C       IVPSOL(N,FCN,T,Y,TEND,TOL,HMAX,H,KFLAG)
C                             Initial value problem (IVP)
C                             integrator
C         N                   Number of first-order ODE's
C         FCN                 Right-hand side of the ODE's system
C                             ( see above )
C         T                   Input: Starting point of integration
C                             T.LT.TEND
C                             Output: Achieved final point of
C                             integration
C         Y(N)                Input and Output: Values at T
C         TEND                Input: Prescribed final point of
C                             integration
C         TOL                 Input: Prescribed relative precision
C                             (>0)
C         HMAX                Input: Maximum permitted stepsize
C         H                   Input: Initial stepsize guess
C                             Output: Stepsize proposal for next
C                             integration step ( H.EQ.0 ,  if
C                             IVPSOL fails to proceed )
C         KFLAG               Input: Print parameter
C                             Output: Error flag ( KFLAG.LT.0
C                             indicates an error ) .
C                             For further details, see IVPSOL .
C
C     Input parameters (* marks inout parameters)
C     ===========================================
C
C       N                     Number of first-order ordinary
C                             differential equations.
C       M                     Number of Shooting nodes.
C                             =2    Single Shooting
C                             >2    Multiple Shooting
C       T(M)                  Single/Multiple Shooting Nodes
C                             ( T(1)= A ,  T(M)= B )
C     * X(N,M)                Start data for Newton iteration.
C       EPS                   Required relative precision of
C                             solution.
C       IOPT(5)               Options array. The options have the
C                               meanings as described below:
C         IOPT(1) - MAXIT     Maximum permitted number of
C                             iteration steps.
C         IOPT(2) - NONLIN    Boundary value problem
C                             classification by user:
C                             0     Linear boundary value problem.
C                             1     Nonlinear boundary value
C                                   problem. Good initial data
C                                   available.
C                             2     Highly nonlinear boundary
C                                   value problem. Only bad
C                                   initial data available. Small
C                                   initial damping factor in
C                                   Gauss Newton method.
C                             3     Highly nonlinear boundary
C                                   value problem. Only bad
C                                   initial data available. Small
C                                   initial damping factor in
C                                   Gauss Newton method.
C                                   Additionally initial rank
C                                   reduction to separable linear
C                                   boundary conditions.
C         IOPT(3) - IBVPSL    Solution method switch:
C                             0     use local linear solver with
C                                   condensing algorithm
C                             1     use global sparse linear solver
C                                   (formerly realized as "BVPSOG")
C         IOPT(4) - IPRINT    Print parameter:
C                             -1    No print
C                              0    Print initial data, iterative
C                                   values of level functions,
C                                   solution data (or final data,
C                                   respectively)
C                             +1    Additionally print iterates
C                                   T(J),X(I,J),  I = 1,...,N ,  J
C                                   = 1,...,M
C         IOPT(5) - LUPRI     Print output unit:
C                             The unit where the output, controlled
C                             by IOPT(4), is written to.
C                             If set to 0, output unit will be 6.
C       NRW                   Dimension of real workspace RW
C                             If IOPT(3) = 0: 
C                               NRW.GE.N*N*(M+5)+10*M*N+10*N+M
C                             If IOPT(3) = 1:
C                               NRW.GE.N*N*(M+1)+12*N*M+4*N+M-1+LICN
C                             with:
C                               NZ = N*N*(M+1)+N*(M-1)
C                               LICN = 2*NZ
C                               LIRN = DMIN1(DMAX1(1.5*NZ,NZ+4*M*N),
C                                            LICN)
C       RW(NRW)               Real workspace
C
C       NIW                   Dimension of integer workspace IW
C                             If IOPT(3) = 0: 
C                               NIW.GE.2*N*N+4*N
C                             If IOPT(3) = 1:
C                               NIW.GE.2*N*N+3*N+8*M*N+8*M*N
C                                      +2*NZ+LICN+LIRN
C       IW(NIW)               Integer workspace
C
C
C
C     Output parameters:
C     ==================
C
C       X(N,M)                Solution data ( or final data,
C                             respectively )
C       INFO                  Information output parameter
C                              >0   Number of iterations performed
C                                   to obtain the solution
C                              <0   BVPSOL termination
C                              -1   If IOPT(3) = 0:
C                                     Iteration stops at stationary
C                                     point
C                                   If IOPT(3) = 1:
C                                     Gaussian elimination failed due 
C                                     to singular Jacobian
C                              -2   Iteration stops after ITMAX
C                                   iteration steps ( as indicated
C                                   by option ITMAX=IOPT(1) )
C                              -3   Integrator failed to complete
C                                   the trajectory computation
C                              -4   Gauss Newton method failed to
C                                   converge
C                              -5   Given initial values
C                                   inconsistent with separable
C                                   linear boundary conditions
C                              -6   If IOPT(3) = 0:
C                                     Iterative refinement failed to
C                                     converge
C                                   If IOPT(3) = 1:
C                                     Termination since Multiple
C                                     Shooting condition or
C                                     condition of Jacobian is too bad
C                              -7   Reliable relative accuracy
C                                   greater than 1.0D-2
C                              -8   Condensing algorithm for
C                                   linear block system fails, use
C                                   global linear solver in
C                                   boundary value problem routine,
C                                   i.e. set IOPT(3)=1
C                              -9   If IOPT(3) = 1:
C                                     Sparse linear solver failed,
C                                     possibly workspace is
C                                     too small
C                             -10   Real or integer work-space
C                                   exhausted
C                             -11   If IOPT(3) = 0:
C                                     Rank reduction failed -
C                                     resulting rank is zero
C
C     ------------------------------------------------------------
C
C*    End Prologue
C:    SMALL = squareroot of "smallest positive machine number
C     divided by relative machine precision"
      DOUBLE PRECISION SMALL
      PARAMETER (SMALL=4.94D-32)
      INTEGER M1,NM,NM1,NN,NRW,NIW,ITMAX,NONLIN,IBVPSL
      DOUBLE PRECISION RELDIF,TOL,XTHR,V
C:    Begin
C     ------------------------------------------------------------
C     1 Internal parameters
C     Standard values fixed below
C     Monitor output unit
      LUMON = IOPT(5)
      IF (LUMON.LE.0.OR.LUMON.GE.100) LUMON=6
C     Scaling threshold
      XTHR = SMALL
C     Prescribed relative precision for numerical integration
      TOL = EPS*1.0D-2
C     Prescribed relative deviation for numerical differentiation
      RELDIF = DSQRT(TOL)
      ITMAX  = IOPT(1)
      NONLIN = IOPT(2)
      IBVPSL = IOPT(3)
      INFO   = IOPT(4)
      IF(INFO.GE.0)THEN
C       Print BVPSOL heading lines
1       FORMAT('1',2X,'B V P S O L',2X,5('*'),2X,'V e r s i o n',2
     *  X,'1 . 2',1X,3('*'),//,1X,'Newton',1('-'),'Method ','for ',
     *  'the ','solution ','of ','boundary ','value ','problems',/
     *  /)
        WRITE(LUMON,1)
      ENDIF
      IF (IBVPSL.EQ.0) THEN
C     Starting value for pseudo - rank of sensitivity matrix E
      IRANK = N
C     Initial preparations
      M1 = M-1
      NN = N*N
      NM = N*M
      NM1 = N*M1
C:    WorkSpace: IW
        L4=1
        L5=L4+N
        L6=L5+N
        L7=L6+N
        L8=L7+N
        L9=L8+N*N
        L10=L9+N*N
        NIW=L10-1
C.    End WorkSpace at NIW
C:    WorkSpace: RW
        L11=1
        L12=L11+N*N*M1
        L13=L12+N*N
        L14=L13+N*N
        L15=L14+N*N
        L16=L15+N*N
        L17=L16+N*N
        L18=L17+NM
        L19=L18+NM
        L20=L19+NM
        L21=L20+NM
        L22=L21+NM
        L23=L22+NM
        L24=L23+NM1
        L25=L24+NM1
        L26=L25+NM1
        L27=L26+NM1
        L28=L27+N
        L29=L28+N
        L30=L29+N
        L31=L30+N
        L32=L31+N
        L33=L32+N
        L34=L33+N
        L35=L34+N
        L36=L35+N
        L37=L36+N
        L38=L37+N
        L39=L38+N
        L40=L39+N
        L41=L40+M
        L42=L41+N
        L43=L42+N*N
        NRW=L43-1
C.    End WorkSpace at NRW
C     ------------------------------------------------------------
C     2 Check for sufficient real/integer workspace
      IF (INFO.GE.0) THEN
2       FORMAT('0','Minimal ','required ','work-space ',':',/,'0',
     *  'Real    ','array ','RW(',I4,')',/,'0','Integer ','array ',
     *  'IW(',I4,')')
        WRITE(LUMON,2)NRW,NIW
      ENDIF
      IF(NRW.LE.IRW.AND.NIW.LE.IIW)THEN
        CALL BVPL(FCN,BC,IVPSOL,N,M,M1,NM,NM1,T,X,EPS,TOL,RELDIF,
     *  NONLIN,IRANK,ITMAX,INFO,XTHR,IW(L4),IW(L5),IW(L6),IW(L7),
     *  IW(L8),IW(L9),RW(L11),RW(L12),RW(L13),RW(L14),RW(L15),RW(
     *  L16),RW(L17),RW(L18),RW(L19),RW(L20),RW(L21),RW(L22),RW(
     *  L23),RW(L24),RW(L25),RW(L26),RW(L27),RW(L28),RW(L29),RW(
     *  L30),RW(L31),RW(L32),RW(L33),RW(L34),RW(L35),RW(L36),RW(
     *  L37),RW(L38),RW(L39),RW(L40),RW(L41),RW(L42),LUMON)
      ELSE
C       Fail exit work-space exhausted
        IF(INFO.GE.0.AND.NRW.GT.IRW)THEN
3         FORMAT('0','Error: ','real    ','work ','- ','space ',
     *    'exhausted',/)
          WRITE(LUMON,3)
        ENDIF
        IF(INFO.GE.0.AND.NIW.GT.IIW)THEN
4         FORMAT('0','Error: ','integer ','work ','- ','space ',
     *    'exhausted',/)
          WRITE(LUMON,4)
        ENDIF
        INFO = -10
      ENDIF
      ELSE IF (IBVPSL.EQ.1) THEN
C     Initial preparations
      M1 = M-1
      NN = N*N
      NM = N*M
      NM1 = N*M1
      NMX8 = 8*NM
      NZ = N*N*(M+1)+N*(M-1)
      LICNQ = 2*NZ
      LIRNQ = MAX0(3*NZ/2,NZ+4*M*N)
      LIRNQ = MIN0(LIRNQ,LICNQ)
      V = DBLE(LIRNQ)/DBLE(LICNQ)
      NRW = N*N*(M+1)+12*M*N+4*N+M-1
      NI2W = 2*NZ+8*M*N
      NIW = NMX8+2*N*N+3*N
      L1 = IRW-NRW
      L2 = IDINT(DBLE(IIW-NIW-NI2W)/(V+1.0D0))
      LICN = MIN0(L1,L2)
      LIRN = IDINT(V*DBLE(LICN))
      LISNQ = LICNQ+LIRNQ
      NKEEP = NMX8
      MINI2W = NI2W + LISNQ
C:    WorkSpace: I2W
        L4=1
        L5=L4+LIRN
        L6=L5+LICN
        L7=L6+NZ
        L8=L7+NZ
        L9=L8+NKEEP
        NI2W=L9-1
C.    End WorkSpace at NI2W
C:    WorkSpace: IW
        L10=L9
        L11=L10+N
        L12=L11+N
        L13=L12+N
        L14=L13+N*N
        L15=L14+N*N
        L16=L15+NMX8
        NIW=L16-1
C.    End WorkSpace at NIW
C:    WorkSpace: RW
        L17=1
        L18=L17+N*N*M1
        L19=L18+N*N
        L20=L19+N*N
        L21=L20+LICN
        L22=L21+NM
        L23=L22+NM
        L24=L23+NM
        L25=L24+NM
        L26=L25+NM
        L27=L26+NM
        L28=L27+NM1
        L29=L28+NM1
        L30=L29+NM1
        L31=L30+NM1
        L32=L31+N
        L33=L32+N
        L34=L33+N
        L35=L34+N
        L36=L35+NM
        L37=L36+NM
        L38=L37+N
        L39=L38+N
        L40=L39+N
        L41=L40+N
        L42=L41+M1
        NRW=L42-1
C.    End WorkSpace at NRW
C     ------------------------------------------------------------
C     2 Check for sufficient real/integer workspace
      IF (INFO.GE.0) THEN
5       FORMAT('0','Minimal ','required ','work-space ',':',/,'0',
     *  'Real          ','array ','RW( ',I5,')',/,'0',
     *  'Integer       ','array ','IW( ',I5,')')
        WRITE(LUMON,5)NRW-LICN+LICNQ,NIW-LIRN-LICN+LISNQ
      ENDIF
      IF(NRW.LE.IRW.AND.NIW.LE.IIW
     *   .AND.LICNQ.LE.LICN+1.AND.LIRNQ.LE.LIRN+1)THEN
        CALL BVPG(FCN,BC,IVPSOL,N,M,M1,NM,NM1,NMX8,NZ,LICN,LIRN,
     *  LISNQ,NKEEP,T,X,EPS,TOL,RELDIF,NONLIN,ITMAX,INFO,XTHR,IW(
     *  L10),IW(L11),IW(L12),IW(L13),IW(L14),IW(L15),RW(L17),RW(
     *  L18),RW(L19),RW(L20),RW(L21),RW(L22),RW(L23),RW(L24),RW(
     *  L25),RW(L26),RW(L27),RW(L28),RW(L29),RW(L30),RW(L31),RW(
     *  L32),RW(L33),RW(L34),RW(L35),RW(L36),RW(L37),RW(L38),RW(
     *  L39),RW(L40),RW(L41),IW(L6),IW(L7),IW(L5),IW(L4),IW(
     *  L8),LUMON)
      ELSE
C       Fail exit work-space exhausted
        IF(INFO.GE.0.AND.NRW.GT.IRW)THEN
6         FORMAT('0','Error: ',A,'work ','- ',
     *    'space ','exhausted',/)
          WRITE(LUMON,6) 'real          '
        ENDIF
        IF(INFO.GE.0.AND.NIW.GT.IIW)THEN
          WRITE(LUMON,6) 'integer       '
        ENDIF
        INFO = -10
      ENDIF
      ELSE
        WRITE(LUMON,*) ' Invalid option IOPT(3) set - must be 0 or 1'
      ENDIF
      RETURN
C     End of driver routine BVPSOL
      END
      SUBROUTINE BVPL(FCN,BC,IVPSOL,N,M,M1,NM,NM1,T,X,EPS,TOL,
     *RELDIF,NONLIN,IRANK,ITMAX,INFO,XTHR,IROW,ICOL,ICOLB,PIVOT,IA,
     *IB,G,A,B,BG,E,QE,DX,DDX,DXQ,DXQA,XA,XW,XU,HH,DHH,HHA,D,DE,R,
     *DR,RA,U,DU,QU,X1,XM,T1,T2,DX1,RF,US,EH,LUMON)
      IMPLICIT DOUBLEPRECISION(S)
      EXTERNAL FCN,BC,IVPSOL
      INTEGER N,M,M1,NM,NM1,LUMON
      DOUBLE PRECISION T(M),X(NM)
      DOUBLE PRECISION EPS
      DOUBLE PRECISION TOL,RELDIF
      INTEGER NONLIN,ITMAX,IRANK
      INTEGER INFO
      DOUBLE PRECISION XTHR
      INTEGER IROW(N),ICOL(N),ICOLB(N),PIVOT(N)
      INTEGER IA(N,N),IB(N,N)
      DOUBLE PRECISION G(N,N,M1)
      DOUBLE PRECISION A(N,N),B(N,N),BG(N,N),E(N,N),QE(N,N)
      DOUBLE PRECISION DX(NM),DDX(NM),DXQ(NM),DXQA(NM),XA(NM),XW(
     *NM),XU(NM1),HH(NM1),DHH(NM1),HHA(NM1),D(N),DE(N),R(N),DR(N),
     *RA(N),U(N),DU(N),QU(N),X1(N),XM(N),T1(N),T2(N),DX1(N),RF(M),
     *US(N)
      DOUBLE PRECISION EH(N,N)
C
C     Additional dimensional integer variables
C     ========================================
C
C       M1                M-1
C       NM                N*M
C       NM1               N*(M-1)
C
C     Internal real arrays (workspace) :
C     ==================================
C
C       G(N,N,M1)        (N,N) -Wronskian Matrices G(1),...,G(M-1)
C                         .
C       A(N,N)            Wronskian Matrix on left boundary
C                         dBC/dX(X(1,...,N),T(1)).
C       B(N,N)            Wronskian Matrix on right boundary
C                         dBC/dX(X((N-1)*M+1,...,N*M),T(M)).
C       BG(N,N)           Workspace for subroutine RHS1. Holds
C                         subsequently the matrices B,B*G(M-1),...,
C                         B*G(M-1)*...*G(2)during computation of
C                         the right hand side of the condensed
C                         linear system.
C       E(N,N)            Sensitivity Matrix of the boundary value
C                         problem: E = A+B*G(M-1)*...*G(1).
C       EH(N,N)           Holds a copy of the row- and column
C                         scaled, but not decomposed sensitivity
C                         matrix E needed for iterative refinement
C                         computations.
C       QE(N,N)           Workspace for DECCON and SOLCON to hold
C                         the updating part of the pseudoinverse
C                         of E in case of it's rank defeciency.
C       D(N)              Diagonal elements of upper triangular
C                         matrix of QR-decomposition computed by
C                         DECCON .
C       DDX(NM)           Workspace for subroutine BLSOLI - gets
C                         the iterative refinement vector of DXQ .
C       DE(N)             Holds row scaling factors for the
C                         sensitivity matrix.
C       DHH(NM1)          Holds the recursive refinement of HH
C                         computed in BLSOLI .
C       DR(N)             Workspace for subroutine BLSOLI to hold
C                         the boundary residual
C                         BC(DXQ(1,...,N),DXQ((M-1)*N+1,...,M*N))+
C                         (A*DXQ(1,...,N))+B*DXQ((M-1)*N+1,...,M*N)
C                         .
C       DU(N)             Workspace for subroutine BLSOLI to hold
C                         the right hand side of the linear system
C                         E*T2 = U-E*DX1 solved to compute the
C                         iterative refinement for DX1 .
C       DX(NM)            Actual newton correction.
C       DXQ(NM)           Simplified Newton correction J(k-1)*X(k)
C                         with the Jacobian J(k) and the iterate
C                         vector X(k) at the k-th iterate.
C       DXQA(NM)          Previous simplified Newton correction
C                         J(k-2)*X(k-1).
C       DX1(N)            Workspace to receive the solution output
C                         of SOLCON within computation of DXQ and
C                         it's iterative refinement.
C       HH(NM1)           Elements (J-1)*N+1 to J*N are holding
C                         the values
C                         Y(T(J+1),X((J-1)*N+1,...,J*N))-X(J*N+1,
C                         ...,(J+1)*N)
C                         ( with the trajectory Y in
C                         ³ T(J),T(J+1)ü , J = 1,...,M-1 ).
C       HHA(NM1)          Holds the previous values of HH .
C       QU(N)             Savespace to hold U(N)for restoring it
C                         in the case of rank reduction.
C       R(N)              Value of the boundary condition function
C                         BC for the current iterate.
C       RA(N)             Previous values of R .
C       RF(M)             Workspace for subroutine BLSOLI - R(J)
C                          gets the maximum norm of the
C                         components (J-1)*N,...,J*N of the
C                         iterative refinement vector to DXQ for
C                         determination of the sweep index for
C                         subsequent iterations.
C       T1(N)             Workspace used for miscellaneous
C                         purposes temporarely.
C       T2(N)             Workspace used for miscellaneous
C                         purposes temporarely.
C       U(N)              Holds the right hand side of the
C                         condensed linear system computed by
C                         subroutine BLRHS1 .
C       US(N)             Workspace for subroutine BLSOLI to save
C                         U(N).
C       XA(NM)            Previous Newton iterate.
C       XU(NM1)           Elements (J-1)*N+1 to J*N are holding
C                         the values Y(T(J+1),X((J-1)*N+1,...,J*N))
C                         of the trajectory in the interval ³ T(J),
C                         T(J+1)ü , (for J = 1,...,M-1 ).
C       XW(NM)            Scaling factors for iteration vector.
C       X1(N)             Components of the iteration vector
C                         corresponding to the left boundary
C                         A = T(1).
C       XM(N)             Components of the iteration vector
C                         corresponding to the right boundary
C                         B = T(M).
C
C     Internal integer arrays (workspace)
C     ===================================
C
C       IROW(N)           Row permutations of sensitivity matrix.
C       ICOL(N)           Column permutations of matrix A
C                         (left boundary).
C       ICOLB(N)          Column permutations of matrix B
C                         (right boundary).
C       PIVOT(N)          Workspace for DECCON and SOLCON to hold
C                         permutations performed during
C                         QR-decomposition.
C       IA(N,N)           Reflects the sparse structure of matrix
C                         A by values 0, 1.
C       IB(N,N)           Reflects the sparse structure of matrix
C                         B by values 0, 1.
C
C     Internal real variables
C     =======================
C
C       COND              Gets the condition of the sensitivity
C                         matrix determined by DECCON . Not
C                         further used.
C       CONDE             Maximum permitted subcondition of
C                         sensitivity matrix E .
C       COND1             Condition of boundary conditions part of
C                         the decomposed sensitivity matrix E .
C       COND2             Condition of the remaining part of the
C                         sensitivity matrix E .
C       CONV              Scaled maximum norm of DXQ computed by
C                         subroutine BLLVLS . Used for convergence
C                         test.
C       CONVA             Holds the previous value of CONV .
C       DEL               Becomes the Euklid's normsquare of the
C                         defect belonging to the condensed linear
C                         system in case of rank defeciency. Used
C                         for computation of damping factor
C                         predictor.
C       EPSMIN            Smallest reasonable permitted accuracy
C                         EPS that can be prescribed by the user.
C       FC                Actual Gauss Newton iteration damping
C                         factor.
C       FCA               Previous Gauss Newton iteration damping
C                         factor.
C       FCDNM             Used to compute the denominator of the
C                         damping factor FC during computation of
C                         it's predictor, corrector and
C                         aposteriori estimate (in the case of
C                         performing a Rank1 update) .
C       FCH               Temporarely used for storing the new FC
C                         when computing aposteriori estimate.
C       FCMIN             Minimum permitted relaxation factor. If
C                         FC becomes smaller than this value, one
C                         of the following may occur:
C                         a.    Recomputation of the sensitivity
C                               matrix by means of difference
C                               approximation (instead of Rank1
C                               update), if Rank1 - update
C                               previously was used
C                         b.    Rank reduction of sensitivity
C                               matrix E ,  if difference
C                               approximation was used previously
C                               and Rank(E).NE.0
C                         c.    Fail exit otherwise
C       FCMINH            DSQRT(FCMIN). Used for rank
C                         determination computations.
C       FCMIN2            FCMIN**2 . Used for FC-predictor
C                         computation.
C       FCNUM             Gets the numerator of the aposteriori
C                         estimate of FC .
C       FCNUMK            Gets the numerator of the corrector
C                         computation of FC .
C       FCNUMP            Gets the numerator of the predictor
C                         computation of FC .
C       H                 Actual integrator stepsize.
C       HMAX              Maximum permitted integrator stepsize.
C                         Set to the length of the integration
C                         interval, e.g. the distance of the
C                         effected Shooting points.
C       HSAVE             Stepsize saved across the call of the
C                         integrator.
C       HSTART            Start stepsize for integration used by
C                         subroutines BLFCNI and BLDERG .
C       MUE               Temporary value used during computation
C                         of damping factors predictor.
C       REDH              Multi purpose reduction factor. (???)
C       RELDIF            Relative deviation for numerical
C                         differentation.
C       SENS1             Sensitivity of boundary conditions part
C                         of the decomposed sensitivity matrix E .
C       SENS2             Sensitivity of the remaining part of the
C                         matrix E .
C       SIGMA             Decision parameter for Jacobian Rank1
C                         updates (SIGMA.GT.1) . Rank1 updates are
C                         inhibited, if SIGMA.GT.1/FCMIN is set.
C       SKAP              Used to compute and print out the
C                         incompatibility factor of the nonlinear
C                         boundary value (e.g. least squares)
C                         problem.
C       SUMF              Standard level of the current iterate,
C                         e.g. Norm2(F(X))**2
C                         with the nonlinear model function F on
C                         which Newton iteration is performed,
C                         arising from the Multiple Shooting
C                         approach.
C       SUMX              Natural level of the current iterate,
C                         e.g. Norm2(DX)
C                         with the Newton correcture DX
C                         (see above).
C       SUMXA             Natural level of the previous iterate.
C       TFAIL             Used to get and print out in case of an
C                         integrator failure the last reached T
C                         value as a proposal for insertion of a
C                         new Shooting point.
C       TOL               Prescribed relative precision for
C                         numerical integration.
C       TOLH              Temporary used for computation of TOL
C                         (may be obmitted|).
C       TOLMIN            Lower bound value for TOL .
C       XTHR              Threshold for scaling.
C       TJ                Used by BLFCNI to hold T(J).
C       TJ1               Used by BLFCNI to hold T(J+1).
C       CORR              Used by BLSOLI to compute RF(J),  J = 1,
C                         ...,M .
C       EPX1              Used by BLSOLI to get the maximum norm
C                         of DX1 corrected by iterative
C                         refinement.
C       EPDX1             Used by BLSOLI to get the maximum norm
C                         of the iterative refinement correcture
C                         for DX1 .
C       EPX1H             Relative accuracy of DX1 after iterative
C                         refinement (= EPDX1/EPX1 ). If EPX1H.GT.
C                         1/2 ,  refinement is considered to be
C                         failed.
C       EPH               Used by BLSOLI as tolerance for the
C                         iterative refinement vectors maximum
C                         norm. If it exceeds at index J the
C                         tolerance, refinement will be performed
C                         next time for the components associated
C                         to Shooting points T(J),...,T(M).
C       SIGDEL            Used by BLSOLI to compute the required
C                         integrator accuracy from the multiple
C                         shooting condition.
C       SIGDLH            Used by BLSOLI to compute the multiple
C                         shooting condition
C                         Max( RF(J+1)/RF(J) ) ,  ( J = 1,...,M-1
C                         ). See above for RF .
C
C     Internal integer variables
C     ==========================
C
C       IC                Permutated index. Used by BLSOLI .
C       ICA               Temporarely used during search for
C                         separable linear boundary conditions at
C                         T(1).
C       ICB               Temporarely used during search for
C                         separable linear boundary conditions at
C                         T(M).
C       IH                Temporarely used during search for
C                         separable linear boundary conditions at
C                         T(1) and T(M).
C       IR                Temporary storage for a row permutation
C                         index.
C       IRANK             Rank of sensitivity matrix E of current
C                         iteration step.
C       IRANKA            Rank of sensitivity matrix E of previous
C                         iteration step.
C       IRANKB            Rank of the decomposed sensitivity
C                         matrix part belonging to the boundary
C                         conditions.
C       IREPET            Parameter of subroutines DECCON and
C                         SOLCON indicating the mode of operation
C                         of these routines:
C                         >=0   Full QR-decomposition required
C                          <0   Rank reduction required
C       IRKMAX            Holds the maximum applied rank of all
C                         previous iterates. Must be necessary
C                         .EQ. NE for convergence.
C       IS                Additional DO loop index.
C       ISUM              Used for determination of sparse
C                         structure of matrices A and B as
C                         nonzeros location counter.
C       ITER              Iteration count.
C       JA                Previous sweep index. Used in subroutine
C                         BLSOLI .
C       JIN               Current sweep index. Used in subroutine
C                         BLSOLI .
C       JJ                Used as "reverse DO loop" index:
C                         JJ = IUPB-J in a loop like DO J = 1,IUPB
C                         ...
C       JN                New sweep index determined in subroutine
C                         BLSOLI .
C       JRED              Damping factor reduction count during an
C                         iterate.
C       KC                Temporary storage for a column
C                         permutation index.
C       KFLAG             Gets the subintervall number of the
C                         failure from subroutine BLDERG
C                         if the integrator failed.
C       KOUNT             Trajectory evaluations count.
C       KPRINT            Print parameter - copy of input
C                         parameter INFO .
C       LEVEL             Flow control parameter needed by
C                         subroutine BLSOLI :
C                         0     indicates computation of Gauss
C                               Newton correcture,
C                         1     indicates computation of
C                               simplified Gauss Newton correcture
C                               (after computation of the
C                               preliminary new iterate)
C       NB                Number of separable boundary conditions
C                         at T(M)
C       NE                Number of not separable boundary
C                         conditions at T(1) (and number of rows
C                         and columns of sensitivity matrix)
C       NEW               Count of subsequent performed Rank1
C                         (Broyden) updates.
C       NY                Iterative refinement sweep count for an
C                         iterate ( used in subroutine BLSOLI ) .
C       NYMAX             Highest allowed iterative refinement
C                         sweep index.
C:    End Parameter
C:    EPMACH = relative machine precision
      DOUBLE PRECISION EPMACH,SMALL
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      DOUBLE PRECISION REDH
      PARAMETER (REDH=1.0D-2)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
      DOUBLE PRECISION EIGHT
      PARAMETER (EIGHT=8.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION FCMIN
      PARAMETER (FCMIN=1.0D-2)
      DOUBLE PRECISION FCMIN2
      PARAMETER (FCMIN2=1.0D-4)
      DOUBLE PRECISION FCNLIN
      PARAMETER (FCNLIN=1.0D-2)
      DOUBLE PRECISION SIGMA
      PARAMETER (SIGMA=2.0D0)
      INTEGER I,ICA,ICB,IH,IR,IRANKA,IRANKB,IREPET,IRKMAX,IS,
     *ISUM,ITER,I0,J,JJ,JN,JRED,J1,K,KC,KFLAG,
     *KOUNT,KPRINT,L,LEVEL,NB,NE,NEW,NY,NYMAX
      DOUBLE PRECISION COND,COND1,COND2,CONV,CONVA,DEL,EPH,
     *EPSMIN,FC,FCA,FCDNM,FCH,FCMINH,FCNMP2,FCNUM,FCNUMK,FCNUMP,
     *HSTART,MUE,S,SENS1,SENS2,SIGDEL,SIGDLH,SKAP,
     *SMALIN,ST,SUMF,SUMX,SUMXA,TFAIL,TH,TOLH,TOLMIN,
     *EPX1H,CONDE
      LOGICAL DIFAPP,FCOMPT,JACRFR,NEXT,REDUCT
      INTEGER L1,L2
      DOUBLE PRECISION S1
C:    Begin
C:    Begin of Segment Bvpsol.Body
C       ----------------------------------------------------------
C       1 Initialization
C       ----------------------------------------------------------
C       1.1 Internal parameters
C       Standard values fixed below
C       Minimum relative precision of integrator ( to be adapted )
        CALL ZIBCONST(EPMACH,SMALL)
        TOLMIN = EPMACH*TEN*TEN
C       Maximum permitted number of iterative refinements sweeps
C       Maximum permitted sub - condition number of senitivity
C         matrix E
        CONDE = ONE/EPMACH
        NYMAX = M-1
C
        FCMINH = DSQRT(FCMIN)
C       ----------------------------------------------------------
C       1.1.1 Common parameters
C       Starting value of relaxation factor (FCMIN.LE.FC.LE.1.0)
        IF(NONLIN.LE.1)THEN
C         for linear or mildly nonlinear problems
          FC = ONE
        ELSE
C         for highly nonlinear problems
          FC = FCNLIN
        ENDIF
C       Starting value for pseudo - rank of matrix A
        IRANK = N
C       Minimum reasonable value for EPS
        EPSMIN = DSQRT(TEN*EPMACH)
        IF(EPS.LT.EPSMIN) EPS = EPSMIN
C       ----------------------------------------------------------
C       1.2 Initial preparations
        IF(FC.LT.FCMIN) FC = FCMIN
        IF(FC.GT.ONE) FC = ONE
        KPRINT = INFO
        ITER = 0
        KOUNT = 0
        IREPET = 0
        INFO = -1000
        FCA = FC
        CONV = ONE
        JACRFR = .FALSE.
C:      Begin SetVec.Vec
        DO 5 L1=1,NM
          XA(L1)=X(L1)
5       CONTINUE
C.      End SetVec.Vec
        IF(TOL.LE.ZERO) TOL = EPS/TEN
        IF(TOL.LT.TOLMIN) TOL = TOLMIN
        DIFAPP = .TRUE.
        HSTART =(T(2)-T(1))*REDH
        SENS1 = ZERO
        COND1 = ONE
        SENS2 = ZERO
        COND2 = ONE
        IRKMAX = 0
        IRANKB = 0
        SIGDLH = ZERO
        SUMF = ONE
C:      Mat IA = Scalar (Rows 1,N ; Cols 1,N)
        L1 = 0
        DO 6 L2=1,N
        DO 6 L44=1,N
          IA(L2,L44)=L1
6       CONTINUE
C.      End SetIntMat.S
C:      Mat IB = Scalar (Rows 1,N ; Cols 1,N)
        L1 = 0
        DO 7 L2=1,N
        DO 7 L44=1,N
          IB(L2,L44)=L1
7       CONTINUE
C.      End SetIntMat.S
C:      CubeMat G (layer 1)= Scalar (Rows 1,N ; Cols 1,N)
        S1 = ZERO
        DO 8 L1=1,N
        DO 8 L2=1,N
          G(L1,L2,1)=S1
8       CONTINUE
C.      End SetCubeMat.S
        IF(KPRINT.GE.0)THEN
C         Print Start vector data, predescribed precision and max
C         iteration steps
9         FORMAT('0','Initial ','data',//)
          WRITE(LUMON,9)
          DO 10 J=1,M
11          FORMAT(D13.5,2X)
            WRITE(LUMON,11)T(J)
12          FORMAT((14X,3(D20.10,1X)))
            WRITE(LUMON,12)(X(L1),L1=(J-1)*N+1,J*N)
10        CONTINUE
13        FORMAT('0','N ','=',I2,2X,'M ','=',I2,/,'0',
     *    'Prescribed ','relative ','precision',D10.2,2X,/,'0',
     *    'Maximum ','permitted ','number ','of ','iteration ',
     *    'steps',1X,I3,//,'1')
          WRITE(LUMON,13)N,M,EPS,ITMAX
          IF(KPRINT.EQ.0)THEN
14          FORMAT('0',1X,66('*'))
            WRITE(LUMON,14)
15          FORMAT('0',4X,'It',4X,'Ny',7X,'Levelf',10X,'Levelx',8X,
     *      'Rel.Fc.')
            WRITE(LUMON,15)
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       1.3 Startup step
C       ----------------------------------------------------------
C       1.3.1 Computation of the residual vector
        CALL BLFCNI(IVPSOL,FCN,BC,N,M,NM,NM1,ITER,KPRINT,HSTART,
     *  FCMIN,T,X,X1,XM,T1,XU,HH,R,TOL,FC,FCOMPT,REDUCT,KFLAG,
     *  KOUNT,INFO,LUMON)
C
C       Main iteration loop
C       ===================
C
C:      While (expression)
16      IF(INFO.EQ.-1000)THEN
C:          Begin of Segment Bvpsol.Core
C             ----------------------------------------------------
C             2 Startup of iteration step
              IF(.NOT.JACRFR)THEN
                LEVEL = 0
C               --------------------------------------------------
C               2.1 Scaling of variables X(NM)
                CALL BLSCLE(N,M,NM,NM1,X,XU,XW,XTHR)
                IF(ITER.NE.0)THEN
                  IRANKA = IRANK
C:                Begin SetVec.Vec
                  DO 17 L1=1,NM
                    DXQA(L1)=DXQ(L1)
17                CONTINUE
C.                End SetVec.Vec
C:                FCNUM = Sum of Formula Elements (for 1,NM)
                  FCNUM = 0.0D0
                  DO 18 L1=1,NM
                    FCNUM=FCNUM+((DX(L1)/XW(L1))**2)
18                CONTINUE
C.                End MakeSum.Comp
C:                FCNMP2 = Sum of Formula Elements (for 1,NM)
                  FCNMP2 = 0.0D0
                  DO 19 L1=1,NM
                    FCNMP2=FCNMP2+((DXQ(L1)/XW(L1))**2)
19                CONTINUE
C.                End MakeSum.Comp
                  FCNUMP = FCNUM*FCNMP2
                ENDIF
                IF(ITER.NE.0)THEN
                  IF(IRANK.GE.NB.AND.FC.GT.FCMINH) IRANK = NE
                  IF(IRANK.LT.NB) IRANK = NB
                  TH = FC-ONE
C:                FCDNM = Sum of Formula Elements (for 1,NM)
                  FCDNM = 0.0D0
                  DO 20 L1=1,NM
                    FCDNM=FCDNM+(((DXQ(L1)+TH*DX(L1))/XW(L1))**2)
20                CONTINUE
C.                End MakeSum.Comp
                  FCH = DSQRT(FCNUM/FCDNM)*FC*FC*HALF
C                 ------------------------------------------------
C                 2.1.1 Decision criterion for Jacobian updating
C                       technique:
C                       DIFAPP.EQ..TRUE. numerical
C                       differentiation,
C                       DIFAPP.EQ..FALSE. rank1 updating
                  DIFAPP = FC.LT.FCA.AND.NEW.GT.0.OR.FCH.LT.FC*
     *            SIGMA.OR.IRANKA.LT.IRANK.OR.EPH*REDH.GT.EPS
                  FCA = FC
                  IF(NONLIN.GT.0) FC = DMIN1(FCH,ONE)
                ENDIF
C               --------------------------------------------------
C               2.2 Difference approximation of jacobian matrix A
C                   ( If Difapp.EQ..TRUE. ) or
C                   Rank-1 update of jacobian matrix A ( If Difapp
C                   .EQ..FALSE. )
                CALL BLDERA(BC,N,M,NM,XW,X1,XM,R,T2,A,B,RELDIF)
C               --------------------------------------------------
C               2.3 Determination of sparse structure of matrices
C                   A and B and determination of internal row
C                   scaling of sensitivity matrix E
                ISUM = 0
                DO 21 I=1,N
                  S = ZERO
                  DO 22 K=1,N
                    TH = DABS(A(I,K))*XW(K)
                    IF(S.LT.TH) S = TH
                    TH = DABS(B(I,K))*XW(K+NM1)
                    IF(S.LT.TH) S = TH
22                CONTINUE
                  IF(S.LT.XTHR) S = XTHR
                  DE(I)=SMALL/S
                  DO 23 K=1,N
                    IF(IA(I,K).LE.0)THEN
                      IF(A(I,K).NE.ZERO)THEN
                        IA(I,K)=1
                        ISUM = 1
                      ENDIF
                    ENDIF
                    IF(IB(I,K).LE.0)THEN
                      IF(B(I,K).NE.ZERO)THEN
                        IB(I,K)=1
                        ISUM = 1
                      ENDIF
                    ENDIF
23                CONTINUE
21              CONTINUE
                IF(ISUM.NE.0)THEN
C                 ------------------------------------------------
C                 2.3.1 Determination of row and column
C                       permutation vectors
                  DO 24 I=1,N
                    ICOL(I)=I
                    ICOLB(I)=I
                    IROW(I)=I
24                CONTINUE
C                 ------------------------------------------------
C                 2.3.2 Search for separable linear boundary
C                       conditions at T(1)
                  NE = N
                  DO 25 I=1,N
                      DO 26 K=1,N
                        IF(IB(I,K).NE.0) GOTO 9996
26                    CONTINUE
                      ISUM = 0
                      DO 27 K=1,N
                        IF(IA(I,K).NE.0)THEN
                          ISUM = ISUM+1
                          ICA = K
                        ENDIF
27                    CONTINUE
                      IF(ISUM.LE.1)THEN
                        DO 28 IS=1,N
                          IH = ICOL(IS)
                          IF(IH.EQ.ICA) ICOL(IS)=ICOL(NE)
                          IH = IROW(IS)
                          IF(IH.EQ.I) IROW(IS)=IROW(NE)
28                      CONTINUE
                        ICOL(NE)=ICA
                        IROW(NE)=I
                        NE = NE-1
                        IF(DABS(R(I)).GT.TEN*EPMACH*DABS(X(ICA)))
     *                  THEN
                          INFO = -5
                          GOTO 9998
                        ENDIF
                      ENDIF
9996                CONTINUE
25                CONTINUE
                  IF(KPRINT.GE.0.AND.NE.EQ.0)THEN
29                  FORMAT('0','Warning: ','attempt ','to ',
     *              'solve ','initial ','value ','problem')
                    WRITE(LUMON,29)
                  ENDIF
                  IF(IRANK.GT.NE) IRANK = NE
                  IRANKA = IRANK
C                 ------------------------------------------------
C                 2.3.3 Search for separable linear boundary
C                       conditions at T(M)
                  NB = 0
                ENDIF
                IF(ISUM.NE.0.AND.NE.NE.0)THEN
                  DO 30 I=1,NE
                      IR = IROW(I)
                      DO 31 K=1,N
                        IF(IA(IR,K).NE.0) GOTO 9995
31                    CONTINUE
                      ISUM = 0
                      DO 32 K=1,N
                        IF(IB(IR,K).NE.0)THEN
                          ISUM = ISUM+1
                          ICB = K
                        ENDIF
32                    CONTINUE
                      IF(ISUM.LE.1)THEN
                        NB = NB+1
                        DO 33 IS=1,N
                          IH = ICOLB(IS)
                          IF(IH.EQ.ICB) ICOLB(IS)=ICOLB(NB)
33                      CONTINUE
                        ICOLB(NB)=ICB
                        IROW(I)=IROW(NB)
                        IROW(NB)=IR
                        IF(DABS(R(IR)).GT.TEN*EPMACH*DABS(X(ICB+
     *                  NM1)))THEN
                          INFO = -5
                          GOTO 9998
                        ENDIF
                      ENDIF
9995                CONTINUE
30                CONTINUE
                  IF(KPRINT.GE.0.AND.NB.EQ.N)THEN
34                  FORMAT('0','Warning: ','attempt ','to ',
     *              'solve ','initial ','value ','problem')
                    WRITE(LUMON,34)
                  ENDIF
C                 Initial rank strategy for highly nonlinear
C                   problems
                  IF(NB.LT.NE.AND.ITER.EQ.0.AND.NONLIN.GT.2) IRANK
     *            = NB
                ENDIF
              ENDIF
              JACRFR = .FALSE.
              IF(DIFAPP)THEN
                NEW = 0
                KFLAG = 0
                CALL BLDERG(FCN,N,NE,M,M1,NM,NM1,T,X,XU,XW,T2,
     *          TFAIL,G,ICOL,IVPSOL,HSTART,TOL,RELDIF,KFLAG)
                IF(KFLAG.LT.0)THEN
                  INFO = -3
                  GOTO 9998
                ENDIF
                IF(M.GT.2) KOUNT = KOUNT+N
                IF(M.EQ.2) KOUNT = KOUNT+NE
              ELSE
                NEW = NEW+1
                CALL BLRK1G(N,M,M1,NM,NM1,XW,DX,HH,HHA,T1,G,FCA)
              ENDIF
C             ----------------------------------------------------
C             2.3.4 Computation of sensitivity matrix E =-A+B*G(M1
C                   *..*G(1))
C                   (projections included)
              IF(IRANK.NE.0)THEN
                DO 35 I=1,NE
                  IR = IROW(I)
C:                Mat E(Row I)= Mat B(Row IR)* Scalar  (for 1,N)
                  S1 = DE(IR)
                  DO 36 L1=1,N
                    E(I,L1)=B(IR,L1)*S1
36                CONTINUE
C.                End SetRow.RowxS
35              CONTINUE
                DO 37 JJ=1,M1
                  J = M-JJ
                  DO 38 I=1,NE
                    DO 39 K=1,N
                      S = ZERO
                      DO 40 L=1,N
                        S = S+E(I,L)*G(L,K,J)
40                    CONTINUE
                      T1(K)=S
39                  CONTINUE
                    DO 41 K=1,N
                      E(I,K)=T1(K)
41                  CONTINUE
38                CONTINUE
37              CONTINUE
C               --------------------------------------------------
C               2.4 Prepare solution of the linear system
C               --------------------------------------------------
C               2.4.1 internal row and column scaling of matrix A
C               INTERNAL ROW AND COLUMN SCALING AND PERMUTATION OF
C                 MATRIX E
                DO 42 K=1,NE
                  KC = ICOL(K)
                  S = XW(KC)
                  DO 43 I=1,NE
                    IR = IROW(I)
                    E(I,K)=-(A(IR,KC)*DE(IR)+E(I,KC))*S
43                CONTINUE
42              CONTINUE
C               --------------------------------------------------
C               2.4.2 Save matrix E on EH
C:              Mat EH = Mat E (Rows 1,NE ; Cols 1,NE)
                DO 44 L1=1,NE
                DO 44 L2=1,NE
                  EH(L1,L2)=E(L1,L2)
44              CONTINUE
C.              End SetMat.Mat
              ENDIF
              IRANKB = NB
C             ----------------------------------------------------
C             2.4.3 Monitor for actually applied maximum rank
              IF(IRKMAX.LT.IRANK) IRKMAX = IRANK
C             ----------------------------------------------------
C             2.5 Save values of R(N)and HH((M-1)*N)
              IF(IREPET.EQ.0)THEN
C:              Begin SetVec.Vec
                DO 45 L1=1,N
                  RA(L1)=R(L1)
45              CONTINUE
C.              End SetVec.Vec
C:              Begin SetVec.Vec
                DO 46 L1=1,NM1
                  HHA(L1)=HH(L1)
46              CONTINUE
C.              End SetVec.Vec
              ENDIF
              NEXT = .FALSE.
C
C             Pseudo-rank reduction loop
C             ==========================
C
C:            DO (Until)
47            CONTINUE
C               --------------------------------------------------
C               3 Main-part of iteration step
C               --------------------------------------------------
C               3.1 Solution of the linear system
C               --------------------------------------------------
C               3.1.1 Constrained QR-decomposition of ( ( COMMA NE
C                     NE ) ) - matrix E
                COND = CONDE
                IF(IRANK.GT.0) CALL BLDECC(E,N,N,IRANKB,NE,NE,
     *          IRANK,COND,D,PIVOT,IREPET,QE,T1)
                IF(NONLIN.EQ.0.AND.IRANK.LT.NE)THEN
                  INFO = -8
                  GOTO 9998
                ENDIF
C               --------------------------------------------------
C               3.1.2 evaluation of subcondition and sensitivity
C                     numbers
                COND1 = ONE
                COND2 = ONE
                SENS1 = ZERO
                SENS2 = ZERO
                IF(IRANKB.NE.0)THEN
                  SENS1 = DABS(D(1))
                  COND1 = SENS1/DABS(D(IRANKB))
                ENDIF
                IF(IRANKB.NE.IRANK)THEN
                  SENS2 = DABS(D(IRANKB+1))
                  COND2 = SENS2/DABS(D(IRANK))
                ENDIF
                IF(FCA.GE.1.0D0.AND.FC.GE.1.0D0.AND.ITER.NE.0)THEN
                  IF(IRANKB.NE.IRANK.AND.SENS2.LT.(EPS/REDH)*SMALL)
     *            IRANK = IRANKB
                  IF(IRANKB.NE.0.AND.SENS1.LT.(EPS/REDH)*SMALL)
     *            IRANK = 0
                ENDIF
C               --------------------------------------------------
C               3.1.3 (best) (least squares) solution of linear (N,
C                     N)-system
                CALL BLSOLI(N,M,M1,NM,NM1,LEVEL,NE,NB,IRANK,IRANKB,
     *          IREPET,NYMAX,KPRINT,EPS,REDH,TOLMIN,TOL,RELDIF,EPH,
     *          EPX1H,SIGDEL,SIGDLH,E,EH,HH,DHH,R,A,B,BG,G,QE,U,QU,
     *          DE,DU,T1,T2,US,DX1,D,DDX,DXQ,XW,DR,RF,IROW,ICOL,
     *          ICOLB,PIVOT,NY,INFO,LUMON)
                IF(INFO.NE.-1000) GOTO 9998
C               --------------------------------------------------
C               3.2 Evaluation of scaled natural level function
C                   SUMX
C                   scaled maximum error norm CONV
C                   evaluation of (scaled) standard level function
C                   SUMF ( SUMF only, if KPRINT.GE.0 )
C                   and computation of ordinary newton corrections
C                   DX(N)
                CALL BLLVLS(N,M,NM,NM1,XW,DXQ,HH,R,DE,CONV,SUMX,
     *          SUMF,KPRINT)
C:              Begin SetVec.Vec
                DO 48 L1=1,NM
                  DX(L1)=DXQ(L1)
48              CONTINUE
C.              End SetVec.Vec
                IF(IREPET.EQ.0.AND.IRANK.NE.0)THEN
C:                Begin SetVec.Vec
                  DO 49 L1=1,IRANK
                    QU(L1)=U(L1)
49                CONTINUE
C.                End SetVec.Vec
                ENDIF
C:              Begin SetVec.Vec
                DO 50 L1=1,NM
                  XA(L1)=X(L1)
50              CONTINUE
C.              End SetVec.Vec
                SUMXA = SUMX
                CONVA = CONV
C               --------------------------------------------------
C               3.3 a - priori estimate of relaxation factor FC
                JRED = 0
                REDUCT = .FALSE.
                IF(ITER.NE.0.AND.NONLIN.NE.0)THEN
                  IF(NEW.LE.0.AND.(IRANK.GE.NE.OR.IRANKA.GE.NE)
     *            .OR.IREPET.NE.0)THEN
C                   ----------------------------------------------
C                   3.3.1 Full rank case (independent of preceding
C                         rank) computation of the denominator of
C                         a-priori estimate
C:                  FCDNM = Sum of Formula Elements (for 1,NM)
                    FCDNM = 0.0D0
                    DO 51 L1=1,NM
                      FCDNM=FCDNM+(((DX(L1)-DXQA(L1))/XW(L1))**2)
51                  CONTINUE
C.                  End MakeSum.Comp
                    IF(IRANK.NE.N)THEN
C                     --------------------------------------------
C                     3.4 Rank - deficient case ( if previous rank
C                         was full )
C                     --------------------------------------------
C                     3.4.1 Computation of the projected
C                           denominator of a-priori estimate
C:                    Vec T1 = Scalar (for 1,N)
                      S1 = ZERO
                      DO 52 L1=1,N
                        T1(L1)=S1
52                    CONTINUE
C.                    End SetVec.S
                      IF(IRANK.NE.0)THEN
                        DO 53 L=1,NE
                          K = ICOL(L)
                          DX1(L)=DXQA(K)/XW(K)
53                      CONTINUE
C                       ------------------------------------------
C                       3.4.2 Projection for reduced component DX1
C                             (NE)
                        CALL BLPRJC(N,NE,IRANK,DEL,DX1,D,T2,QE,
     *                  PIVOT)
                        DO 54 L=1,NE
                          K = ICOL(L)
                          T1(K)=DX1(L)*XW(K)
54                      CONTINUE
                      ENDIF
                      DO 55 J=1,M1
                        DO 56 I=1,N
                          S = ZERO
                          DO 57 K=1,N
                            S = S+T1(K)*G(I,K,J)
57                        CONTINUE
                          T2(I)=S
56                      CONTINUE
C:                      Begin SetVec.Vec
                        DO 58 L1=1,N
                          T1(L1)=T2(L1)
58                      CONTINUE
C.                      End SetVec.Vec
                        I0 = J*N
                        DO 59 I=1,N
                          ST = ONE/XW(I+I0)
                          S = T1(I)
                          DEL = DEL+S*ST*ST*(S+(DX(I+I0)-DXQA(I+I0))
     *                    *TWO)
59                      CONTINUE
55                    CONTINUE
                      FCDNM = FCDNM+DEL
                    ENDIF
                    FCDNM = FCDNM*SUMX
C                   ----------------------------------------------
C                   3.4.3 New relaxation factor
                    IF(FCDNM.GE.FCNUMP*FCMIN2)THEN
                      MUE = FCA*DSQRT(FCNUMP/FCDNM)
                      FC = DMIN1(MUE,ONE)
                    ELSE
                      FC = ONE
                    ENDIF
                  ENDIF
                  IREPET = 0
                  REDUCT = FC.LT.FCMIN
                ENDIF
                LEVEL = 1
                  IF(.NOT.REDUCT)THEN
C                   ----------------------------------------------
C                   3.5 Save natural level for later computations
C                       of corrector and print iterate
                    FCNUMK = SUMX
                    IF(KPRINT.GE.0)THEN
C                     Print Standard - and natural level
                      IF(KPRINT.GT.0)THEN
60                      FORMAT('0',1X,66('*'))
                        WRITE(LUMON,60)
61                      FORMAT('0',4X,'It',4X,'Ny',7X,'Levelf',10X,
     *                  'Levelx',18X,'New',4X,'Rank')
                        WRITE(LUMON,61)
                      ENDIF
62                    FORMAT('0',4X,I2,4X,I2,5X,D10.3,2X,4X,D10.3
     *                ,2X,13X,I2,6X,I2)
                      WRITE(LUMON,62)ITER,NY,SUMF,SUMXA,NEW,IRANK
                      IF(KPRINT.GT.0)THEN
63                      FORMAT('0',1X,66('*'))
                        WRITE(LUMON,63)
                      ENDIF
                    ENDIF
C
C                   Relaxation-factor reduction loop
C                   ================================
C
C:                  DO (Until)
64                  CONTINUE
C                     --------------------------------------------
C                     3.6 Preliminary new iterate
C:                    DO (Until)
65                    CONTINUE
                        FCOMPT = .FALSE.
C:                      Vec X = Vec XA + Vec DX * Scalar (for 1,NM)
                        S1 = FC
                        DO 66 L1=1,NM
                          X(L1)=XA(L1)+DX(L1)*S1
66                      CONTINUE
C.                      End SetVec.Vec&VecxS
                        IF(ITER.GT.ITMAX)THEN
                          INFO = -2
                          GOTO 9997
                        ENDIF
C                       ------------------------------------------
C                       3.6.1 Computation of the residual vector
                        CALL BLFCNI(IVPSOL,FCN,BC,N,M,NM,NM1,ITER,
     *                  KPRINT,HSTART,FCMIN,T,X,X1,XM,T1,XU,HH,R,
     *                  TOL,FC,FCOMPT,REDUCT,KFLAG,KOUNT,INFO,LUMON)
                        IF(INFO.NE.-1000) GOTO 9997
                        IF(REDUCT) GOTO 9994
                      IF(.NOT.(FCOMPT)) GOTO  65
C.                    UNTIL ( expression - negated above)
C                     --------------------------------------------
C                     3.6.2 (best) (least squares) solution of
C                           linear (N,N) -system
                      CALL BLSOLI(N,M,M1,NM,NM1,LEVEL,NE,NB,IRANK,
     *                IRANKB,IREPET,NYMAX,KPRINT,EPS,REDH,TOLMIN,
     *                TOL,RELDIF,EPH,EPX1H,SIGDEL,SIGDLH,E,EH,HH,
     *                DHH,R,A,B,BG,G,QE,U,QU,DE,DU,T1,T2,US,DX1,D,
     *                DDX,DXQ,XW,DR,RF,IROW,ICOL,ICOLB,PIVOT,NY,
     *                INFO,LUMON)
                      IF(INFO.NE.-1000) GOTO 9998
C                     --------------------------------------------
C                     3.6.3 Evaluation of scaled natural level
C                           function SUMX
C                           scaled maximum error norm CONV and
C                           evaluation of (scaled) standard level
C                           function SUMF
                      CALL BLLVLS(N,M,NM,NM1,XW,DXQ,HH,R,DE,CONV,
     *                SUMX,SUMF,KPRINT)
C                     --------------------------------------------
C                     3.7 Rank independent convergence test
                      IF(CONV.LE.EPS.AND.IRKMAX.EQ.NE)THEN
                        INFO = 0
                        GOTO 9997
                      ENDIF
C                     --------------------------------------------
C                     3.8 Natural monotonicity test
                      IF(SUMX.GT.SUMXA)THEN
C                       ------------------------------------------
C                       3.9 Output of iterate
                        IF(KPRINT.GE.0)THEN
C                         Print Standard - and natural level, and
C                         damping factor
                          IF(KPRINT.GT.0)THEN
67                          FORMAT('0',1X,66('*'))
                            WRITE(LUMON,67)
68                          FORMAT('0',4X,'It',4X,'Ny',7X,'Levelf',
     *                      10X,'Levelx',8X,'Rel.Fc.')
                            WRITE(LUMON,68)
                          ENDIF
69                        FORMAT('0',4X,I2,4X,I2,5X,D10.3,2X,4X,D10.3
     *                    ,2X,4X,F5.3)
                          WRITE(LUMON,69)ITER,NY,SUMF,SUMX,FC
                          IF(KPRINT.GT.0)THEN
70                          FORMAT('0',1X,66('*'))
                            WRITE(LUMON,70)
                          ENDIF
                        ENDIF
                        JRED = JRED+1
                        IF(NONLIN.EQ.0)THEN
                          INFO = -4
                          GOTO 9997
                        ENDIF
C                       ------------------------------------------
C                       3.10 Compute reduced relaxation factor FC
                        TH = FC-ONE
C:                      FCDNM = Sum of Formula Elements (for 1,NM)
                        FCDNM = 0.0D0
                        DO 71 L1=1,NM
                          FCDNM=FCDNM+(((DXQ(L1)+TH*DX(L1))/XW(L1))
     *                    **2)
71                      CONTINUE
C.                      End MakeSum.Comp
                        FC = DSQRT(FCNUMK/FCDNM)*FC*FC*HALF
C                       Rank reduction, if relaxation factor to
C                         small
                        REDUCT = FC.LT.FCMIN.OR.NEW.GT.0.AND.JRED
     *                  .GT.1
                      ELSE
                        NEXT = .TRUE.
                      ENDIF
                    IF(.NOT.(NEXT.OR.REDUCT)) GOTO  64
C.                  UNTIL ( expression - negated above)
C
C                   End of relaxation-factor reduction loop
C                   =======================================
C
                  ENDIF
9994            CONTINUE
                IF(.NOT.NEXT)THEN
C                 ------------------------------------------------
C                 3.11 Restore former values for repeting
C                      iteration step
                  IREPET = 1
C                 Restore former values
                  LEVEL = 0
C:                Begin SetVec.Vec
                  DO 72 L1=1,N
                    R(L1)=RA(L1)
72                CONTINUE
C.                End SetVec.Vec
C:                Begin SetVec.Vec
                  DO 73 L1=1,N
                    X1(L1)=XA(L1)
73                CONTINUE
C.                End SetVec.Vec
C:                Begin SetVec.Vec
                  DO 74 L1=1,N
                    XM(L1)=XA(L1+NM1)
74                CONTINUE
C.                End SetVec.Vec
C:                Begin SetVec.Vec
                  DO 75 L1=1,NM
                    X(L1)=XA(L1)
75                CONTINUE
C.                End SetVec.Vec
C:                Begin SetVec.Vec&Vec
                  DO 76 L1=1,NM1
                    XU(L1)=X(L1+N)+HHA(L1)
76                CONTINUE
C.                End SetVec.Vec&Vec
C:                Begin SetVec.Vec
                  DO 77 L1=1,NM1
                    HH(L1)=HHA(L1)
77                CONTINUE
C.                End SetVec.Vec
                  IF(KPRINT.GE.0)THEN
78                  FORMAT('0',5X,I2,1X,'Not ','accepted ',
     *              'relaxation ','factor',5X,F5.3,12X,I2)
                    WRITE(LUMON,78)ITER,FC,IRANK
                  ENDIF
                  IF(ITER.EQ.0)THEN
                    FC = FCMIN
                  ENDIF
                  IF(NEW.GT.0)THEN
                    DIFAPP = .TRUE.
                    JACRFR = .TRUE.
                    IRANK = NE
                    IF(IRANK.LT.NB) IRANK = NB
                    GOTO 9998
                  ENDIF
C                 ------------------------------------------------
C                 3.12 Pseudo-rank reduction
                  IREPET = -1
                  IF(IRANK.EQ.0)THEN
                    INFO = -11
                    GOTO 9997
                  ENDIF
C:                Begin SetVec.Vec
                  DO 79 L1=1,IRANK
                    U(L1)=QU(L1)
79                CONTINUE
C.                End SetVec.Vec
                  IRANK = IRANK-1
                  IF(IRANKB.GT.IRANK) IRANKB = IRANK
                ENDIF
              IF(.NOT.(NEXT)) GOTO  47
C.            UNTIL ( expression - negated above)
C
C             End of pseudo-rank reduction loop
C             =================================
C
C             ----------------------------------------------------
C             4 Preparations to start the following iteration step
              ITER = ITER+1
              IRANKA = IRANK
C             Preliminary pseudo-rank
              IF(IRANK.GE.NB.AND.FC.GT.FCMINH) IRANK = NE
              IF(IRANK.LT.NB) IRANK = NB
C             ----------------------------------------------------
C             4.1 Print values
              IF(KPRINT.GE.0)THEN
C               Print Standard - and natural level, and damping
C               factor
                IF(KPRINT.GT.0)THEN
80                FORMAT('0',1X,66('*'))
                  WRITE(LUMON,80)
81                FORMAT('0',4X,'It',4X,'Ny',7X,'Levelf',10X,
     *            'Levelx',8X,'Rel.Fc.')
                  WRITE(LUMON,81)
                ENDIF
82              FORMAT('0',4X,I2,4X,I2,5X,D10.3,2X,4X,D10.3,2X,4X,F5.3)
                WRITE(LUMON,82)ITER,NY,SUMF,SUMX,FC
                IF(KPRINT.GT.0)THEN
83                FORMAT('0',1X,66('*'))
                  WRITE(LUMON,83)
                  DO 84 J=1,M
85                  FORMAT(D13.5,2X)
                    WRITE(LUMON,85)T(J)
86                  FORMAT((14X,3(D20.10,1X)))
                    WRITE(LUMON,86)(X(L1),L1=(J-1)*N+1,J*N)
84                CONTINUE
                ENDIF
              ENDIF
9997        CONTINUE
C.          End of Segment Bvpsol.Core
9998      CONTINUE
        GOTO 16
        ENDIF
C.      EndWhile
C
C       End of main iteration loop
C       ==========================
C
C       ----------------------------------------------------------
C       5 Exits
C       ----------------------------------------------------------
C       5.1 Solution exit
        IF(INFO.EQ.0)THEN
          ITER = ITER+1
C:        Vec X = Vec X + Vec DXQ (for 1,NM)
          DO 87 L1=1,NM
            X(L1)=X(L1)+DXQ(L1)
87        CONTINUE
C.        End SetVec.&Vec
          INFO = ITER
          IF(KPRINT.LT.0)THEN
            GOTO 9999
          ENDIF
          IF(KPRINT.GT.0)THEN
C           Print levels, damping factor of last iteration step
88          FORMAT('0',1X,66('*'))
            WRITE(LUMON,88)
89          FORMAT('0',4X,'It',4X,'Ny',7X,'Levelf',10X,'Levelx',8X,
     *      'Rel.Fc.')
            WRITE(LUMON,89)
          ENDIF
90        FORMAT('0',4X,I2,4X,I2,5X,D10.3,2X,4X,D10.3,2X,4X,F5.3)
          WRITE(LUMON,90)ITER,NY,SUMF,SUMX,FC
91        FORMAT('0',1X,66('*'))
          WRITE(LUMON,91)
92        FORMAT('1')
          WRITE(LUMON,92)
          IF(IRANK.LT.NE)THEN
            INFO = -1
          ELSE
C           Print solution info
93          FORMAT('0','Solution ','of',1X,'boundary ','value ',
     *      'problem',' obtained',/,'0','BVPSOL',' required',I3,1X,
     *      'Iteration ','steps ','with',I4,1X,'trajectory',
     *      ' evaluations',//)
            WRITE(LUMON,93)ITER,KOUNT
            CALL BLPRCV(LUMON,CONV,EPH)
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       5.2 Fail exit messages
C       ----------------------------------------------------------
C       5.2.1 Rank-deficiency : best least squares solution of bvp
C             obtained
        IF(INFO.EQ.-1.AND.KPRINT.GE.0)THEN
94        FORMAT('0','Iteration ','terminates ','at ',
     *    'stationary ','point',/)
          WRITE(LUMON,94)
          CALL BLPRCV(LUMON,CONVA,EPH)
          IF(ITER.NE.0)THEN
            SKAP = ZERO
            IF(FCA.EQ.ONE.AND.FC.EQ.ONE.AND.IRANKA.EQ.IRANK) SKAP
     *      = DSQRT(SUMXA/FCNUMK)
            IF(SKAP.GT.ZERO)THEN
95            FORMAT('0','Incompatibility ','factor ','kappa',D10.3
     *        ,2X,/)
              WRITE(LUMON,95)SKAP
            ENDIF
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       5.2.2 Termination after more than itmax iterations
        IF(INFO.EQ.-2.AND.KPRINT.GE.0)THEN
96        FORMAT('0','Iteration ','terminates ','after ','itmax ',
     *    '=',I3,2X,'iteration ','steps')
          WRITE(LUMON,96)ITMAX
        ENDIF
C       ----------------------------------------------------------
C       5.2.3 Integrator failed in trajectory computation
        IF(INFO.EQ.-3.AND.KPRINT.GE.0)THEN
97        FORMAT('0','Integrator ','failed ','in ','trajectory ',
     *    'computation ',/)
          WRITE(LUMON,97)
          J1 =-KFLAG
98        FORMAT('0','BVPSOL ','terminates',/,'Subinterval',I3,1X,
     *    'possibly ','insert ','new ','node',D20.11,2X,/)
          WRITE(LUMON,98)J1,TFAIL
        ENDIF
C       ----------------------------------------------------------
C       5.2.4 Convergence fail of Gauss - Newton method
        IF(INFO.EQ.-4.AND.KPRINT.GE.0)THEN
99        FORMAT('0','Gauss ','Newton ','method ','fails ','to ',
     *    'converge',/)
          WRITE(LUMON,99)
        ENDIF
C       ----------------------------------------------------------
C       5.2.5 Inconsistent initial data
        IF(INFO.EQ.-5.AND.KPRINT.GE.0)THEN
100       FORMAT('0','Error: ','initial ','data ','and ',
     *    'boundary ','conditions ','are ','inconsistent',/)
          WRITE(LUMON,100)
        ENDIF
C       ----------------------------------------------------------
C       5.2.6 Convergence fail of iterative refinement sweeps
        IF(INFO.EQ.-6)THEN
          IF(KPRINT.GE.0)THEN
101         FORMAT('0','Termination ','since ','iterative ',
     *      'refinement ','fails ','to ','converge',/,2X,'Insert ',
     *      'new ','nodes',/)
            WRITE(LUMON,101)
          ENDIF
          JN = JN-1
          IF(JN.GT.0)THEN
102         FORMAT('0',8X,'in ','subinterval',2X,I3,/)
            WRITE(LUMON,102)JN
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       5.2.7 Insufficient error tolerance for integrator
        IF(INFO.EQ.-7.AND.KPRINT.GE.0)THEN
          TOLH = EPS/SIGDEL
          RELDIF = DSQRT(TOLH/SIGDEL)
103       FORMAT('0','Suggested ','integrator ','accuracy',D10.1
     *    ,2X,/,'0','Suggested ','relative ','deviation ',
     *    'parameter',D10.1,2X,/)
          WRITE(LUMON,103)TOLH,RELDIF
104       FORMAT('0','Reduce ','relative ','error ','tolerance ',
     *    'for ','integrator ','to',D10.1,2X,/,2X,'or ','insert ',
     *    'new ','nodes',/)
          WRITE(LUMON,104)TOLH
          S = REDH/TOL
          DO 105 J=1,M1
            IF(RF(J).GT.S)THEN
106           FORMAT(2X,'in ','subinterval',I3,/)
              WRITE(LUMON,106)J
            ENDIF
105       CONTINUE
107       FORMAT('0','Reliable ','relative ','accuracy ',
     *    'greater ','than',1X,D6.1,2X,/)
          WRITE(LUMON,107)1.0D-2
        ENDIF
C       ----------------------------------------------------------
C       5.2.8 ill - conditioned condensed linear system
        IF(INFO.EQ.-8.AND.KPRINT.GE.0)THEN
108       FORMAT('0','Gaussian ','block ','elimination ','fails',/,
     *    2X,'by ','ill ','- ','conditioned ','condensed ',
     *    'linear ','system',/)
          WRITE(LUMON,108)
          IF(IRANK.EQ.NE)THEN
109         FORMAT('0','Relative ','accuracy ','of ','DX1',D10.3
     *      ,2X)
            WRITE(LUMON,109)EPX1H
          ENDIF
110       FORMAT('0','Possibly set IOPT(3)=1',/)
          WRITE(LUMON,110)
        ENDIF
C       ----------------------------------------------------------
C       5.2.9 rank reduced to zero
        IF(INFO.EQ.-11.AND.KPRINT.GE.0)THEN
11001     FORMAT('0','rank ','reduction ','failed - ',
     *    2X,'resulting ','rank ','is ','zero ',/)
          WRITE(LUMON,11001)
        ENDIF
C       ----------------------------------------------------------
C       5.3 Common exit
        IF(KPRINT.GE.0)THEN
C
          J1 = 1
          SMALIN = ONE/SMALL
          IF(IRANKB.NE.0) 
     $      CALL BLPRCD(LUMON,COND1,SENS1,SMALIN,J1,IRANKB)
          IF(IRANKB.NE.IRANK)THEN
            J1 = IRANKB+1
            CALL BLPRCD(LUMON,COND2,SENS2,SMALIN,J1,IRANK)
          ENDIF
111       FORMAT('0','Multiple ','shooting ','condition',D10.3,2X,
     *    /,'1')
          WRITE(LUMON,111)SIGDLH
          IF(INFO.GT.0)THEN
112         FORMAT('0','Solution ','data',/)
            WRITE(LUMON,112)
          ENDIF
          IF(INFO.LT.0)THEN
113         FORMAT('0','Final ','data',/)
            WRITE(LUMON,113)
          ENDIF
          DO 114 J=1,M
115         FORMAT(D13.5,2X)
            WRITE(LUMON,115)T(J)
116         FORMAT((14X,3(D20.10,1X)))
            WRITE(LUMON,116)(X(L1),L1=(J-1)*N+1,J*N)
114       CONTINUE
        ENDIF
C       End of exits
C       End of subroutine BVPSOL
9999  CONTINUE
C.    End of Segment Bvpsol.Body
      RETURN
      END
      SUBROUTINE BVPG(FCN,BC,IVPSOL,N,M,M1,NM,NM1,NMX8,NZ,LICN,
     *LIRN,LISNQ,NKEEP,T,X,EPS,TOL,RELDIF,NONLIN,ITMAX,INFO,XTHR,
     *IROW,ICOLA,ICOLB,IA,IB,IW,G,A,B,E,WO,DX,DXQ,DXQA,XA,XW,XU,HH,
     *DHH,HHA,DE,R,DR,RA,U,DU,X1,XM,T1,T2,RF,IVECT,JVECT,ICN,IRN,
     *IKEEP,LUMON)
      IMPLICIT DOUBLEPRECISION(S)
      EXTERNAL FCN,BC,IVPSOL
      INTEGER N,M,M1,NM,NM1,NMX8,LISNQ,NKEEP,LUMON
      DOUBLE PRECISION T(M),X(NM)
      DOUBLE PRECISION EPS
      DOUBLE PRECISION TOL,RELDIF
      INTEGER NONLIN,ITMAX
      INTEGER INFO
      DOUBLE PRECISION XTHR
      INTEGER IRN(LIRN),ICN(LICN),IVECT(NZ),JVECT(NZ),IKEEP(
     *NKEEP)
      INTEGER IROW(N),ICOLA(N),ICOLB(N)
      INTEGER IA(N,N),IB(N,N)
      INTEGER IW(NMX8)
      DOUBLE PRECISION G(N,N,M1)
      DOUBLE PRECISION A(N,N),B(N,N)
      DOUBLE PRECISION E(LICN),WO(NM)
      DOUBLE PRECISION DX(NM),DXQ(NM),DXQA(NM),XA(NM),XW(NM),XU(
     *NM1),HH(NM1),DHH(NM1),HHA(NM1),DE(N),R(N),DR(N),RA(N),U(NM),
     *DU(NM),X1(N),XM(N),T1(N),T2(N),RF(M1)
C
C     Addtional dimensional integer variables:
C     ========================================
C
C       M1                M-1
C       NM                N*M
C       NM1               N*(M-1)
C       NMX8              8*N*M
C       LIRN,LICN,NKEEP,NZ
C                         See driver routine BVPSOL
C
C     Internal real arrays (workspace) :
C     ==================================
C
C       G(N,N,M1)        (N,N) -Wronskian Matrices G(1),...,G(M-1)
C                         .
C       A(N,N)            Wronskian Matrix on left boundary
C                         dBC/dX(X(1,...,N),T(1)).
C       B(N,N)            Wronskian Matrix on right boundary
C                         dBC/dX(X((N-1)*M+1,...,N*M),T(M)).
C       E(LICN)           Holds the values of the Jacobian stored
C                         in sparse mode.
C       DE(N)             Holds row scaling factors for the
C                         boundary conditions part of the Jacobian
C                         matrix.
C       DHH(NM1)          Holds the continuity residuals computed
C                         in BGSOLI .
C       DR(N)             Workspace for subroutine BGSOLI to hold
C                         the boundary residual
C                         BC(DXQ(1,...,N),DXQ((M-1)*N+1,...,M*N))+
C                         (A*DXQ(1,...,N))+B*DXQ((M-1)*N+1,...,M*N)
C                         .
C       DU(NM)            Used by BGSOLI . Gets the total residual
C                         for the current iterate.
C       DX(NM)            Actual newton correction.
C       DXQ(NM)           Simplified Newton correction J(k-1)*X(k)
C                         with the Jacobian J(k) and the iterate
C                         vector X(k) at the k-th iterate.
C       DXQA(NM)          Previous simplified Newton correction
C                         J(k-2)*X(k-1).
C       HH(NM1)           Elements (J-1)*N+1 to J*N are holding
C                         the values
C                         Y(T(J+1),X((J-1)*N+1,...,J*N))-X(J*N+1,
C                         ...,(J+1)*N)
C                         ( with the trajectory Y in
C                         [T(J),T(J+1)] , J = 1,...,M-1 ).
C       HHA(NM1)          Holds the previous value of HH .
C       R(N)              Value of the boundary condition function
C                         BC for the current iterate.
C       RA(N)             Previous values of R .
C       RF(M1)            Used by BGSOLI . Gets the norms of the
C                         Wronskian matrices.
C       T1(N)             Workspace used for miscellaneous
C                         purposes temporarely.
C       T2(N)             Workspace used for miscellaneous
C                         purposes temporarely.
C       U(NM)             Gets the right hand side of the linear
C                         system to be solved in each iteration
C                         step. Used in BGSOLI .
C       WO(NM)            Workspace needed for sparse solver. Must
C                         not be altered outside the sparse packet
C                         routines.
C       XA(NM)            Previous Newton iterate.
C       XU(NM1)           Elements (J-1)*N+1 to J*N are holding
C                         the values Y(T(J+1),X((J-1)*N+1,...,J*N))
C                         of the trajectory in the interval
C                         [T(J),T(J+1)] , (for J = 1,...,M-1 ).
C       XW(NM)            Scaling factors for iteration vector.
C       X1(N)             Components of the iteration vector
C                         corresponding to the left boundary
C                         A = T(1).
C       XM(N)             Components of the iteration vector
C                         corresponding to the right boundary
C                         B = T(M).
C
C     Internal integer arrays (workspace)
C     ===================================
C
C       IROW(N)           Row permutations of boundary derivative
C                         matrices A and B .
C       ICOLA(N)          Column permutations of matrix A
C                         (left boundary).
C       ICOLB(N)          Column permutations of matrix B
C                         (right boundary).
C       IA(N,N)           Reflects the sparse structure of matrix
C                         A by values 0, 1.
C       IB(N,N)           Reflects the sparse structure of matrix
C                         B by values 0, 1.
C       IW(NMX8)          Workspace needed for sparse solver
C                         package.
C
C     Internal short integer arrays (workspace)
C     =========================================
C
C       IRN(LIRN)         Workspace for MA28/MA30 sparse package.
C                         On Input to routine MA28A, it must hold
C                         the row indices of the sparse matrix.
C       ICN(LICN)         Workspace for MA28/MA30 sparse package.
C                         On Input to routine MA28A, it must hold
C                         the column indices of the sparse matrix.
C       IVECT(NZ)         Input to routine MA28B: must hold the
C                         row indices of the sparse matrix.
C       JVECT(NZ)         Input to routine MA28B: must hold the
C                         column indices of the sparse matrix.
C       IKEEP(NKEEP)      Workspace array for MA28 sparse package.
C                         To be preserved across the calls of the
C                         routines MA28A,MA28B,MA28C .
C
C     Internal real variables:
C     ========================
C
C       COND              Gets the condition of the Jacobian
C                         matrix computed by BGSOLI .
C       CORR              Gets the 1-norm of the residual DU .
C                         Computed by BGSOLI .
C       CONV              Scaled maximum norm of DXQ computed by
C                         subroutine BGLVLS . Used for convergence
C                         test.
C       CONVA             Holds the previous value of CONV .
C       EPSMIN            Smallest reasonable permitted accuracy
C                         EPS that can be prescribed by the user.
C       FC                Actual Gauss Newton iteration damping
C                         factor.
C       FCA               Previous Gauss Newton iteration damping
C                         factor.
C       FCDNM             Used to compute the denominator of the
C                         damping factor FC during computation of
C                         it's predictor, corrector and
C                         aposteriori estimate (in the case of
C                         performing a Rank1 update) .
C       FCH               Temporarely used for storing the new FC
C                         when computing aposteriori estimate.
C       FCMIN             Minimum permitted relaxation factor. If
C                         FC becomes smaller than this value, one
C                         of the following may occur:
C                         a.    Recomputation of the sensitivity
C                               matrix by means of difference
C                               approximation (instead of Rank1
C                               update), if Rank1 - update
C                               previously was used
C                         b.    Rank reduction of sensitivity
C                               matrix E ,  if difference
C                               approximation was used previously
C                               and Rank(E).NE.0
C                         c.    Fail exit otherwise
C       FCMIN2            FCMIN**2 . Used for FC-predictor
C                         computation.
C       FCNUM             Gets the numerator of the aposteriori
C                         estimate of FC .
C       FCNUMP            Gets the numerator of the predictor
C                         computation of FC .
C       FCNUMK            Gets the numerator of the corrector
C                         computation of FC .
C       H                 Actual integrator stepsize.
C       HMAX              Maximum permitted integrator stepsize.
C                         Set to the length of the integration
C                         interval, e.g. the distance of the
C                         effected Shooting points.
C       HSAVE             Stepsize saved across the call of the
C                         integrator.
C       HSTART            Start stepsize for integration used by
C                         subroutines BLFCNI and BLDERG .
C       MUE               Temporary value used during computation
C                         of damping factors predictor.
C       REDH              Multi purpose reduction factor. (???)
C       RELDIF            Relative deviation for numerical
C                         differentation.
C       SIGMA             Decision parameter for Jacobian Rank1
C                         updates (SIGMA.GT.1) . Rank1 updates are
C                         inhibited, if SIGMA.GT.1/FCMIN is set.
C       SKAP              Used to compute and print out the
C                         incompatibility factor of the nonlinear
C                         boundary value (e.g. least squares)
C                         problem.
C       SUMF              Standard level of the current iterate,
C                         e.g. Norm2(F(X))**2
C                         with the nonlinear model function F on
C                         which Newton iteration is performed,
C                         arising from the Multiple Shooting
C                         approach.
C       SUMX              Natural level of the current iterate,
C                         e.g. Norm2(DX)
C                         with the Newton correcture DX
C                         (see above).
C       SUMXA             Natural level of the previous iterate.
C       TFAIL             Used to get and print out in case of an
C                         integrator failure the last reached T
C                         value as a proposal for insertion of a
C                         new Shooting point.
C       TOL               Prescribed relative precision for
C                         numerical integration.
C       TOLH              Temporary used for computation of TOL
C                         (may be obmitted|).
C       TOLMIN            Lower bound value for TOL .
C       XTHR              Threshold for scaling.
C       TJ                Used by BLFCNI to hold T(J).
C       TJ1               Used by BLFCNI to hold T(J+1).
C       EPH               Gets TOL*SIGDEL by BGSOLI . If EPH.GT.
C                         REDH ,  termination occurs, since
C                         Multiple Shooting condition is too bad.
C       SIGDEL            Used by BGSOLI to compute the required
C                         integrator accuracy from the multiple
C                         shooting condition.
C       SIGDLH            Used by BGSOLI temporary during
C                         determination of SIGDEL .
C
C     Internal integer variables
C     ==========================
C
C       IAF               Indicates, if sparse structure of
C                         Jacobian matrix must be reordered by
C                         sparse solver or not:
C                         0     Not necessary, can call MA28B to
C                               decompose Jacobian
C                         1     Must be done, so use MA28A .
C       IC                Permutated index. Used by BGSOLI .
C       ICA               Temporarely used during search for
C                         separable linear boundary conditions at
C                         T(1).
C       ICB               Temporarely used during search for
C                         separable linear boundary conditions at
C                         T(M).
C       IH                Temporarely used during search for
C                         separable linear boundary conditions at
C                         T(1)and T(M).
C       INZ               Count of nonzero Jacobian matrix
C                         elements.
C       IR                Temporary storage for a row permutation
C                         index.
C       IRANK             Rank of Jacobian matrix E of current
C                         iteration step estimated by sparse
C                         solver package (Common block variable).
C       IS                Additional DO loop index.
C       ISUM              Used for determination of sparse
C                         structure of matrices A and B as
C                         nonzeros location counter.
C       ITER              Iteration count.
C       JJ                Used as "reverse DO loop" index:
C                         JJ = IUPB-J in a loop like DO J = 1,IUPB
C                         ...
C       JRED              Damping factor reduction count during an
C                         iterate.
C       KC                Temporary storage for a column
C                         permutation index.
C       KFLAG             Gets the subintervall number of the
C                         failure from subroutine BLDERG ,  if the
C                         integrator failed.
C       KOUNT             Trajectory evaluations count.
C       KPRINT            Print parameter - copy of input
C                         parameter INFO .
C       LEVEL             Flow control parameter needed by
C                         subroutine BGSOLI :
C                         0     indicates computation of Newton
C                               correcture,
C                         1     indicates computation of
C                               simplified Newton correcture
C                               (after computation of the
C                               preliminary new iterate)
C       NA                Number of separable boundary conditions
C                         at T(1): N-NAQ
C       NAQ               Number of not separable boundary
C                         conditions at T(1)
C       NB                Number of separable boundary conditions
C                         at T(M)
C       NBQ               Number of not separable boundary
C                         conditions at T(M): N-NB
C       NEW               Count of subsequent performed Rank1
C                         (Broyden) updates.
C       NM                Number of rows and columns of the
C                         Jacobian matrix part to be decomposed by
C                         the sparse solver.
C       NRS               N-(NA+NB)
C:    End Parameter
C:    EPMACH = relative machine precision
      DOUBLE PRECISION EPMACH,SMALL
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      DOUBLE PRECISION REDH
      PARAMETER (REDH=1.0D-2)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
      DOUBLE PRECISION EIGHT
      PARAMETER (EIGHT=8.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION FCMIN
      PARAMETER (FCMIN=1.0D-2)
      DOUBLE PRECISION FCMIN2
      PARAMETER (FCMIN2=1.0D-4)
      DOUBLE PRECISION FCNLIN
      PARAMETER (FCNLIN=1.0D-2)
      DOUBLE PRECISION SIGMA
      PARAMETER (SIGMA=2.0D0)
      INTEGER I,ICA,ICB,ICNCP,IH,IR,IRNCP,IRANK,IS,ISUM,ITER,
     *J,JRED,J0,J1,K,KFLAG,KOUNT,KPRINT,
     *L,LEVEL,MINIRN,MINICN,NB,NDIM,NAQ,NEW
      DOUBLE PRECISION COND,CORR,CONV,CONVA,EPH,EPSMIN,
     *EPX1H,EPSQ,FC,FCA,FCDNM,FCH,FCNMP2,FCNUM,FCNUMK,FCNUMP,
     *HSTART,MUE,RESID,RMIN,S,SIGDEL,SIGDLH,SUMF,
     *SUMX,SUMXA,TFAIL,TH,TOLH,TOLMIN,UQ
      LOGICAL ABORT1,ABORT2,DIFAPP,FCOMPT,GROW,IBLOCK,IVFAIL,
     *JACRFR,JACRST,NEXT
      COMMON /MA28ED/ LP,MP,IBLOCK,GROW
      COMMON /MA28FD/ EPSQ,RMIN,RESID,IRNCP,ICNCP,MINIRN,MINICN,
     *IRANK,ABORT1,ABORT2
      INTEGER L1,L2
      DOUBLE PRECISION S1
C:    Begin
C:    Begin of Segment BVPSOG.Body
C       ----------------------------------------------------------
C       1 Initialization
C       ----------------------------------------------------------
C       1.1 Internal parameters
C       Standard values fixed below
C       Minimum relative precision of integrator ( to be adapted )
        CALL ZIBCONST(EPMACH,SMALL)
        TOLMIN = EPMACH*TEN*TEN
C       Maximum permitted number of iterative refinements sweeps
C       ----------------------------------------------------------
C       1.1.1 Common parameters
C       Starting value of relaxation factor (FCMIN.LE.FC.LE.1.0)
        IF(NONLIN.LE.1)THEN
C         for linear or mildly nonlinear problems
          FC = ONE
        ELSE
C         for highly nonlinear problems
          FC = FCNLIN
        ENDIF
C       Minimum reasonable value for EPS
        EPSMIN = DSQRT(TEN*EPMACH)
        IF(EPS.LT.EPSMIN) EPS = EPSMIN
C       ----------------------------------------------------------
C       1.2 Initial preparations
        IF(FC.LT.FCMIN) FC = FCMIN
        IF(FC.GT.ONE) FC = ONE
        KPRINT = INFO
        ITER = 0
        KOUNT = 0
        INFO = -1000
        FCA = FC
        CONV = ZERO
        JACRFR = .FALSE.
        JACRST = .FALSE.
C:      Begin SetVec.Vec
        DO 6 L1=1,NM
          XA(L1)=X(L1)
6       CONTINUE
C.      End SetVec.Vec
        IF(TOL.LE.ZERO) TOL = EPS/TEN
        IF(TOL.LT.TOLMIN) TOL = TOLMIN
        DIFAPP = .TRUE.
        HSTART =(T(2)-T(1))*REDH
        LP = 0
        MP = 0
        IBLOCK = .FALSE.
        EPSQ = 0.1D0
        UQ = 0.1D0
        SIGDLH = ZERO
        SUMF = ONE
C:      Mat IA = Scalar (Rows 1,N ; Cols 1,N)
        L1 = 0
        DO 7 L2=1,N
        DO 7 L43=1,N
          IA(L2,L43)=L1
7       CONTINUE
C.      End SetIntMat.S
C:      Mat IB = Scalar (Rows 1,N ; Cols 1,N)
        L1 = 0
        DO 8 L2=1,N
        DO 8 L43=1,N
          IB(L2,L43)=L1
8       CONTINUE
C.      End SetIntMat.S
C:      CubeMat G (layer 1)= Scalar (Rows 1,N ; Cols 1,N)
        S1 = ZERO
        DO 9 L1=1,N
        DO 9 L2=1,N
          G(L1,L2,1)=S1
9       CONTINUE
C.      End SetCubeMat.S
        IF(KPRINT.GE.0)THEN
C         Print Start vector data, predescribed precision and max
C         iteration steps
10        FORMAT('0','Initial ','data',//)
          WRITE(LUMON,10)
          DO 11 J=1,M
12          FORMAT(D13.5,2X)
            WRITE(LUMON,12)T(J)
13          FORMAT((14X,3(D20.10,1X)))
            WRITE(LUMON,13)(X(L1),L1=(J-1)*N+1,J*N)
11        CONTINUE
14        FORMAT('0','N ','=',I2,2X,'M ','=',I2,/,'0',
     *    'Prescribed ','relative ','precision',D10.2,2X,/,'0',
     *    'Maximum ','permitted ','number ','of ','iteration ',
     *    'steps',1X,I3,//,'1')
          WRITE(LUMON,14)N,M,EPS,ITMAX
          IF(KPRINT.EQ.0)THEN
15          FORMAT('0',1X,66('*'))
            WRITE(LUMON,15)
16          FORMAT('0',4X,'It',7X,'Levelf',10X,'Levelx',8X,
     *      'Rel.Fc.')
            WRITE(LUMON,16)
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       1.3 Startup step
C       ----------------------------------------------------------
C       1.3.1 Computation of the residual vector
        CALL BLFCNI(IVPSOL,FCN,BC,N,M,NM,NM1,ITER,KPRINT,HSTART,
     *  FCMIN,T,X,X1,XM,T1,XU,HH,R,TOL,FC,FCOMPT,IVFAIL,KFLAG,
     *  KOUNT,INFO,LUMON)
C
C       Main iteration loop
C       ===================
C
C:      While (expression)
17      IF(INFO.EQ.-1000)THEN
C:          Begin of Segment BVPSOG.Core
C             ----------------------------------------------------
C             2 Startup of iteration step
              IF(.NOT.(JACRFR.OR.JACRST))THEN
                LEVEL = 0
C               --------------------------------------------------
C               2.1 Scaling of variables X(NM)
                CALL BLSCLE(N,M,NM,NM1,X,XU,XW,XTHR)
                IF(ITER.NE.0)THEN
C:                Begin SetVec.Vec
                  DO 18 L1=1,NM
                    DXQA(L1)=DXQ(L1)
18                CONTINUE
C.                End SetVec.Vec
C:                FCNUM = Sum of Formula Elements (for 1,NM)
                  FCNUM = 0.0D0
                  DO 19 L1=1,NM
                    FCNUM=FCNUM+((DX(L1)/XW(L1))**2)
19                CONTINUE
C.                End MakeSum.Comp
C:                FCNMP2 = Sum of Formula Elements (for 1,NM)
                  FCNMP2 = 0.0D0
                  DO 20 L1=1,NM
                    FCNMP2=FCNMP2+((DXQ(L1)/XW(L1))**2)
20                CONTINUE
C.                End MakeSum.Comp
                  FCNUMP = FCNUM*FCNMP2
                ENDIF
                IF(ITER.NE.0)THEN
                  TH = FC-ONE
C:                FCDNM = Sum of Formula Elements (for 1,NM)
                  FCDNM = 0.0D0
                  DO 21 L1=1,NM
                    FCDNM=FCDNM+(((DXQ(L1)+TH*DX(L1))/XW(L1))**2)
21                CONTINUE
C.                End MakeSum.Comp
                  FCH = DSQRT(FCNUM/FCDNM)*FC*FC*HALF
C                 ------------------------------------------------
C                 2.1.1 Decision criterion for Jacobian updating
C                       technique:
C                       DIFAPP.EQ..TRUE. numerical
C                       differentiation,
C                       DIFAPP.EQ..FALSE. rank1 updating
                  DIFAPP = FC.LT.FCA.AND.NEW.GT.0.OR.FCH.LT.FC*
     *            SIGMA
                  FCA = FC
                  IF(NONLIN.GT.0) FC = DMIN1(FCH,ONE)
                ENDIF
C               --------------------------------------------------
C               2.2 Difference approximation of jacobian matrix A
C                   ( If Difapp.EQ..TRUE. ) or
C                   Rank-1 update of jacobian matrix A ( If Difapp
C                   .EQ..FALSE. )
                CALL BLDERA(BC,N,M,NM,XW,X1,XM,R,T2,A,B,RELDIF)
C               --------------------------------------------------
C               2.3 Determination of sparse structure of matrices
C                   A and B
                IAF = 0
                DO 22 I=1,N
                  S = ZERO
                  DO 23 K=1,N
                    TH = DABS(A(I,K))*XW(K)
                    S = S+TH
                    TH = DABS(B(I,K))*XW(K+NM1)
                    S = S+TH
23                CONTINUE
                  IF(S.LT.XTHR) S = XTHR
                  DE(I)=ONE/S
                  DO 24 K=1,N
                    IF(IA(I,K).LE.0)THEN
                      IF(A(I,K).NE.ZERO)THEN
                        IA(I,K)=1
                        IAF = 1
                      ENDIF
                    ENDIF
                    IF(IB(I,K).LE.0)THEN
                      IF(B(I,K).NE.ZERO)THEN
                        IB(I,K)=1
                        IAF = 1
                      ENDIF
                    ENDIF
24                CONTINUE
22              CONTINUE
                IF(IAF.NE.0)THEN
C                 ------------------------------------------------
C                 2.3.1 Determination of row and column
C                       permutation vectors
                  DO 25 I=1,N
                    ICOLA(I)=I
                    ICOLB(I)=I
                    IROW(I)=I
25                CONTINUE
C                 ------------------------------------------------
C                 2.3.2 Search for separable linear boundary
C                       conditions at T(1)
                  NAQ = N
                  DO 26 I=1,N
                      DO 27 K=1,N
                        IF(IB(I,K).NE.0) GOTO 9996
27                    CONTINUE
                      ISUM = 0
                      DO 28 K=1,N
                        IF(IA(I,K).NE.0)THEN
                          ISUM = ISUM+1
                          ICA = K
                        ENDIF
28                    CONTINUE
                      IF(ISUM.LE.1)THEN
                        DO 29 IS=1,N
                          IH = ICOLA(IS)
                          IF(IH.EQ.ICA) ICOLA(IS)=ICOLA(NAQ)
                          IH = IROW(IS)
                          IF(IH.EQ.I) IROW(IS)=IROW(NAQ)
29                      CONTINUE
                        ICOLA(NAQ)=ICA
                        IROW(NAQ)=I
                        NAQ = NAQ-1
                        IF(DABS(R(I)).GT.TEN*EPMACH*DABS(X(ICA)))
     *                  THEN
                          INFO = -5
                          GOTO 9998
                        ENDIF
                      ENDIF
9996                CONTINUE
26                CONTINUE
                  IF(KPRINT.GE.0.AND.NAQ.EQ.0)THEN
30                  FORMAT('0','Warning: ','attempt ','to ',
     *              'solve ','initial ','value ','problem')
                    WRITE(LUMON,30)
                  ENDIF
                  NA = N-NAQ
C                 ------------------------------------------------
C                 2.3.3 Search for separable linear boundary
C                       conditions at T(M)
                  NB = 0
                ENDIF
                IF(IAF.NE.0.AND.NAQ.NE.0)THEN
                  DO 31 I=1,NAQ
                      IR = IROW(I)
                      DO 32 K=1,N
                        IF(IA(IR,K).NE.0) GOTO 9995
32                    CONTINUE
                      ISUM = 0
                      DO 33 K=1,N
                        IF(IB(IR,K).NE.0)THEN
                          ISUM = ISUM+1
                          ICB = K
                        ENDIF
33                    CONTINUE
                      IF(ISUM.LE.1)THEN
                        NB = NB+1
                        DO 34 IS=1,N
                          IH = ICOLB(IS)
                          IF(IH.EQ.ICB) ICOLB(IS)=ICOLB(NB)
34                      CONTINUE
                        ICOLB(NB)=ICB
                        IROW(I)=IROW(NB)
                        IROW(NB)=IR
                        IF(DABS(R(IR)).GT.TEN*EPMACH*DABS(X(ICB+
     *                  NM1)))THEN
                          INFO = -5
                          GOTO 9998
                        ENDIF
                      ENDIF
9995                CONTINUE
31                CONTINUE
                  IF(KPRINT.GE.0.AND.NB.EQ.N)THEN
35                  FORMAT('0','Warning: ','attempt ','to ',
     *              'solve ','initial ','value ','problem')
                    WRITE(LUMON,35)
                  ENDIF
                  NBQ = N-NB
C                 ------------------------------------------------
C                 2.3.4 Count non-zeroes in jacobian and store
C                       their locations
                  NDIM = NM-NA-NB
                  NRS = N-NA-NB
                  INZ = 0
                  DO 36 I=1,N
                      IF(NAQ.NE.0)THEN
                        DO 37 K=1,NAQ
                          INZ = INZ+1
                          IVECT(INZ)=I
                          JVECT(INZ)=K
37                      CONTINUE
                      ENDIF
                      II = I
                      IF(M.EQ.2)THEN
                        II = 0
                        IF(NBQ.EQ.0) GOTO 9994
                        DO 38 L=1,NBQ
                          IF(ICOLB(NB+L).EQ.I) II = L
38                      CONTINUE
                        IF(II.EQ.0) GOTO 9994
                      ENDIF
                      INZ = INZ+1
                      IVECT(INZ)=I
                      JVECT(INZ)=NAQ+II
9994                CONTINUE
36                CONTINUE
                  IF(M.NE.2)THEN
                    DO 39 J=2,M1
                      DO 40 I=1,N
                          DO 41 K=1,N
                            INZ = INZ+1
                            IVECT(INZ)=(J-1)*N+I
                            JVECT(INZ)=(J-1)*N+K-NA
41                        CONTINUE
                          II = I
                          IF(J.EQ.M1)THEN
                            II = 0
                            IF(NBQ.EQ.0) GOTO 9993
                            DO 42 L=1,NBQ
                              IF(ICOLB(NB+L).EQ.I) II = L
42                          CONTINUE
                            IF(II.EQ.0) GOTO 9993
                          ENDIF
                          INZ = INZ+1
                          IVECT(INZ)=(J-1)*N+I
                          JVECT(INZ)=J*N+II-NA
9993                    CONTINUE
40                    CONTINUE
39                  CONTINUE
                  ENDIF
                  IF(NRS.NE.0)THEN
                    DO 43 I=1,NRS
                      II = I+NB
                      IF(NAQ.NE.0)THEN
                        DO 44 K=1,NAQ
                          IF(IA(IROW(II),ICOLA(K)).NE.0)THEN
                            INZ = INZ+1
                            IVECT(INZ)=NM1+I
                            JVECT(INZ)=K
                          ENDIF
44                      CONTINUE
                      ENDIF
                      IF(NBQ.NE.0)THEN
                        DO 45 KK=1,NBQ
                          K = KK+NB
                          IF(IB(IROW(II),ICOLB(K)).NE.0)THEN
                            INZ = INZ+1
                            IVECT(INZ)=NM1+I
                            JVECT(INZ)=NM1+KK-NA
                          ENDIF
45                      CONTINUE
                      ENDIF
43                  CONTINUE
                  ENDIF
                ENDIF
              ENDIF
              JACRFR = .FALSE.
              IF(.NOT.JACRST)THEN
                IF(DIFAPP)THEN
                  NEW = 0
                  KFLAG = 0
                  CALL BLDERG(FCN,N,NAQ,M,M1,NM,NM1,T,X,XU,XW,T2,
     *            TFAIL,G,ICOLA,IVPSOL,HSTART,TOL,RELDIF,KFLAG)
                  IF(KFLAG.LT.0)THEN
                    INFO = -3
                    GOTO 9998
                  ENDIF
                  IF(M.GT.2) KOUNT = KOUNT+N
                  IF(M.EQ.2) KOUNT = KOUNT+NAQ
                ELSE
                  NEW = NEW+1
                  CALL BLRK1G(N,M,M1,NM,NM1,XW,DX,HH,HHA,T1,G,FCA)
                ENDIF
              ENDIF
              JACRST = .FALSE.
C             ----------------------------------------------------
C             2.3.5 Storing of total sparse jacobian
C                   ( including row and column scaling )
              INZ = 0
              DO 46 I=1,N
                IF(NAQ.NE.0)THEN
                  DO 47 K=1,NAQ
                    INZ = INZ+1
                    E(INZ)=-G(I,ICOLA(K),1)*XW(ICOLA(K))/XW(I+N)
47                CONTINUE
                ENDIF
                II = I
                  IF(M.EQ.2)THEN
                    II = 0
                    IF(NBQ.EQ.0) GOTO 9992
                    DO 48 L=1,NBQ
                      IF(ICOLB(NB+L).EQ.I) II = L
48                  CONTINUE
                    IF(II.EQ.0) GOTO 9992
                  ENDIF
                  INZ = INZ+1
                  E(INZ)=ONE
9992            CONTINUE
46            CONTINUE
              IF(M.NE.2)THEN
                DO 49 J=2,M1
                  J0 =(J-1)*N
                  J1 = J0+N
                  DO 50 I=1,N
                    DO 51 K=1,N
                      INZ = INZ+1
                      E(INZ)=-G(I,K,J)*XW(K+J0)/XW(I+J1)
51                  CONTINUE
                    II = I
                      IF(J.EQ.M1)THEN
                        II = 0
                        IF(NBQ.EQ.0) GOTO 9991
                        DO 52 L=1,NBQ
                          IF(ICOLB(NB+L).EQ.I) II = L
52                      CONTINUE
                        IF(II.EQ.0) GOTO 9991
                      ENDIF
                      INZ = INZ+1
                      E(INZ)=ONE
9991                CONTINUE
50                CONTINUE
49              CONTINUE
              ENDIF
              IF(NRS.NE.0)THEN
                DO 53 I=1,NRS
                  II = I+NB
                  IF(NAQ.NE.0)THEN
                    DO 54 K=1,NAQ
                      IF(IA(IROW(II),ICOLA(K)).NE.0)THEN
                        INZ = INZ+1
                        E(INZ)=-A(IROW(II),ICOLA(K))*DE(IROW(II))*
     *                  XW(ICOLA(K))
                      ENDIF
54                  CONTINUE
                  ENDIF
                  IF(NBQ.NE.0)THEN
                    DO 55 KK=1,NBQ
                      K = KK+NB
                      IF(IB(IROW(II),ICOLB(K)).NE.0)THEN
                        INZ = INZ+1
                        E(INZ)=-B(IROW(II),ICOLB(K))*DE(IROW(II))*
     *                  XW(ICOLB(K)+NM1)
                      ENDIF
55                  CONTINUE
                  ENDIF
53              CONTINUE
              ENDIF
C             ----------------------------------------------------
C             2.4 Save values of R(N)and HH((M-1)*N)
C:            Begin SetVec.Vec
              DO 56 L1=1,N
                RA(L1)=R(L1)
56            CONTINUE
C.            End SetVec.Vec
C:            Begin SetVec.Vec
              DO 57 L1=1,NM1
                HHA(L1)=HH(L1)
57            CONTINUE
C.            End SetVec.Vec
              NEXT = .FALSE.
C             ----------------------------------------------------
C             3 Main-part of iteration step
C             ----------------------------------------------------
C             3.1 Solution of the linear system
C             ----------------------------------------------------
C             3.1.1 Decomposition of (N,N)-matrix A
C             ----------------------------------------------------
C             3.1.2 LU-decomposition of(NDIM,NDIM)-MATRIX E
              IF(IAF.EQ.1)THEN
                DO 58 I=1,INZ
                  IRN(I)=IVECT(I)
                  ICN(I)=JVECT(I)
58              CONTINUE
                CALL MA28AD(NDIM,INZ,E,LICN,IRN,LIRN,ICN,UQ,IKEEP,
     *          IW,WO,IFLAG)
                IAF = 0
              ELSE
                CALL MA28BD(NDIM,INZ,E,LICN,IVECT,JVECT,ICN,IKEEP,
     *          IW,WO,IFLAG)
                IF(RMIN.LT.1.0D-4) IAF = 1
                IF(IAF.EQ.1)THEN
                  JACRST = .TRUE.
                  GOTO 9998
                ENDIF
              ENDIF
              IF(IFLAG.EQ.-1.OR.IFLAG.EQ.-2)THEN
                INFO = -1
                GOTO 9998
              ENDIF
              IF(IFLAG.LE.-3)THEN
                INFO = -9
                GOTO 9998
              ENDIF
C             ----------------------------------------------------
C             3.1.3 Solution of linear (N,N)-system
              CALL BGSOLI(N,M,M1,NM,NM1,NDIM,LICN,NKEEP,NA,NAQ,NB,
     *        NBQ,NRS,ITER,LEVEL,KPRINT,EPS,REDH,TOLMIN,FC,FCA,TOL,
     *        RELDIF,EPH,EPX1H,SIGDEL,SIGDLH,COND,CORR,HH,DHH,R,A,
     *        B,G,U,DE,DU,T1,DXQ,XW,DR,RF,WO,E,IROW,ICOLA,ICOLB,
     *        ICN,IKEEP,INFO,LUMON)
              IF(INFO.NE.-1000) GOTO 9998
C             ----------------------------------------------------
C             3.2 Evaluation of scaled natural level function SUMX
C                 scaled maximum error norm CONV
C                 evaluation of (scaled) standard level function
C                 SUMF ( SUMF only, if KPRINT.GE.0 )
C                 and computation of ordinary newton corrections
C                 DX(N)
              CALL BGLVLS(N,M,NM,NM1,XW,DXQ,HH,R,DE,CONV,SUMX,SUMF,
     *        KPRINT)
C:            Begin SetVec.Vec
              DO 59 L1=1,NM
                DX(L1)=DXQ(L1)
59            CONTINUE
C.            End SetVec.Vec
C:            Begin SetVec.Vec
              DO 60 L1=1,NM
                XA(L1)=X(L1)
60            CONTINUE
C.            End SetVec.Vec
              SUMXA = SUMX
              CONVA = CONV
C             ----------------------------------------------------
C             3.3 a - priori estimate of relaxation factor FC
              JRED = 0
              IF(ITER.NE.0.AND.NONLIN.NE.0)THEN
                IF(NEW.EQ.0)THEN
C                 ------------------------------------------------
C                 3.3.1 Computation of the denominator of a-priori
C                       estimate
C:                FCDNM = Sum of Formula Elements (for 1,NM)
                  FCDNM = 0.0D0
                  DO 61 L1=1,NM
                    FCDNM=FCDNM+(((DX(L1)-DXQA(L1))/XW(L1))**2)
61                CONTINUE
C.                End MakeSum.Comp
C                 ------------------------------------------------
C                 3.3.2 New relaxation factor
                  FCDNM = FCDNM*SUMX
                  IF(FCDNM.GE.FCNUMP*FCMIN2)THEN
                    MUE = FCA*DSQRT(FCNUMP/FCDNM)
                    FC = DMIN1(MUE,ONE)
                  ELSE
                    FC = ONE
                  ENDIF
                ENDIF
                IF(FC.LT.FCMIN)THEN
                  INFO = -4
                  GOTO 9997
                ENDIF
              ENDIF
              LEVEL = 1
C             ----------------------------------------------------
C             3.4 Save natural level for later computations of
C                 corrector and print iterate
              FCNUMK = SUMX
              IF(KPRINT.GE.0)THEN
C               Print Standard - and natural level
                IF(KPRINT.GT.0)THEN
62                FORMAT('0',1X,66('*'))
                  WRITE(LUMON,62)
63                FORMAT('0',4X,'It',7X,'Levelf',10X,'Levelx',18X,
     *            'New')
                  WRITE(LUMON,63)
                ENDIF
64              FORMAT('0',4X,I2,5X,D10.3,2X,4X,D10.3,2X,13X,I2)
                WRITE(LUMON,64)ITER,SUMF,SUMXA,NEW
                IF(KPRINT.GT.0)THEN
65                FORMAT('0',1X,66('*'))
                  WRITE(LUMON,65)
                ENDIF
              ENDIF
C
C             Relaxation-factor reduction loop
C             ================================
C
C:            DO (Until)
66            CONTINUE
C               --------------------------------------------------
C               3.5 Preliminary new iterate
C:              DO (Until)
67              CONTINUE
                  FCOMPT = .FALSE.
C:                Vec X = Vec XA + Vec DX * Scalar (for 1,NM)
                  S1 = FC
                  DO 68 L1=1,NM
                    X(L1)=XA(L1)+DX(L1)*S1
68                CONTINUE
C.                End SetVec.Vec&VecxS
                  IF(ITER.GT.ITMAX)THEN
                    INFO = -2
                    GOTO 9997
                  ENDIF
C                 ------------------------------------------------
C                 3.5.1 Computation of the residual vector
                  CALL BLFCNI(IVPSOL,FCN,BC,N,M,NM,NM1,ITER,KPRINT,
     *            HSTART,FCMIN,T,X,X1,XM,T1,XU,HH,R,TOL,FC,FCOMPT,
     *            IVFAIL,KFLAG,KOUNT,INFO,LUMON)
                  IF(IVFAIL)THEN
                    INFO = -4
                    GOTO 9997
                  ENDIF
                  IF(INFO.NE.-1000) GOTO 9997
                IF(.NOT.(FCOMPT)) GOTO  67
C.              UNTIL ( expression - negated above)
C               --------------------------------------------------
C               3.5.2 Solution of linear (N,N)-system
                CALL BGSOLI(N,M,M1,NM,NM1,NDIM,LICN,NKEEP,NA,NAQ,
     *          NB,NBQ,NRS,ITER,LEVEL,KPRINT,EPS,REDH,TOLMIN,FC,
     *          FCA,TOL,RELDIF,EPH,EPX1H,SIGDEL,SIGDLH,COND,CORR,
     *          HH,DHH,R,A,B,G,U,DE,DU,T1,DXQ,XW,DR,RF,WO,E,IROW,
     *          ICOLA,ICOLB,ICN,IKEEP,INFO,LUMON)
                IF(INFO.NE.-1000) GOTO 9998
C               --------------------------------------------------
C               3.5.3 Evaluation of scaled natural level function
C                     SUMX
C                     scaled maximum error norm CONV and
C                     evaluation of (scaled) standard level
C                     function SUMF
                CALL BGLVLS(N,M,NM,NM1,XW,DXQ,HH,R,DE,CONV,SUMX,
     *          SUMF,KPRINT)
C               --------------------------------------------------
C               3.6 Convergence test
                IF(CONV.LE.EPS)THEN
                  INFO = 0
                  GOTO 9997
                ENDIF
C               --------------------------------------------------
C               3.7 Natural monotonicity test
                IF(SUMX.GT.SUMXA)THEN
C                 ------------------------------------------------
C                 3.8 Output of iterate
                  IF(KPRINT.GE.0)THEN
C                   Print Standard - and natural level, and
C                   damping factor
                    IF(KPRINT.GT.0)THEN
69                    FORMAT('0',1X,66('*'))
                      WRITE(LUMON,69)
70                    FORMAT('0',4X,'It',7X,'Levelf',10X,'Levelx',
     *                8X,'Rel.Fc.')
                      WRITE(LUMON,70)
                    ENDIF
71                  FORMAT('0',4X,I2,5X,D10.3,2X,4X,D10.3,2X,4X,F5.3)
                    WRITE(LUMON,71)ITER,SUMF,SUMX,FC
                    IF(KPRINT.GT.0)THEN
72                    FORMAT('0',1X,66('*'))
                      WRITE(LUMON,72)
                    ENDIF
                  ENDIF
                  JRED = JRED+1
                  IF(NONLIN.EQ.0)THEN
                    INFO = -4
                    GOTO 9997
                  ENDIF
C                 ------------------------------------------------
C                 3.9 Compute reduced relaxation factor
                  TH = FC-ONE
C:                FCDNM = Sum of Formula Elements (for 1,NM)
                  FCDNM = 0.0D0
                  DO 73 L1=1,NM
                    FCDNM=FCDNM+(((DXQ(L1)+TH*DX(L1))/XW(L1))**2)
73                CONTINUE
C.                End MakeSum.Comp
                  FC = DSQRT(FCNUMK/FCDNM)*FC*FC*HALF
C                 ------------------------------------------------
C                 3.10 Fail exit, if relaxation factor to small
                  JACRFR = FC.LT.FCMIN.OR.NEW.GT.0.AND.JRED.GT.1
                  IF(JACRFR.AND.NEW.EQ.0)THEN
                    INFO = -4
                    GOTO 9998
                  ENDIF
                ENDIF
              IF(.NOT.(SUMX.LE.SUMXA.OR.JACRFR)) GOTO  66
C.            UNTIL ( expression - negated above)
C
C             End of relaxation-factor reduction loop
C             =======================================
C
              IF(JACRFR)THEN
C               --------------------------------------------------
C               3.11 Restore former values for repeting iteration
C                    step
C               Restore former values
                LEVEL = 0
C:              Begin SetVec.Vec
                DO 74 L1=1,N
                  R(L1)=RA(L1)
74              CONTINUE
C.              End SetVec.Vec
C:              Begin SetVec.Vec
                DO 75 L1=1,N
                  X1(L1)=XA(L1)
75              CONTINUE
C.              End SetVec.Vec
C:              Begin SetVec.Vec
                DO 76 L1=1,N
                  XM(L1)=XA(L1+NM1)
76              CONTINUE
C.              End SetVec.Vec
C:              Begin SetVec.Vec
                DO 77 L1=1,NM
                  X(L1)=XA(L1)
77              CONTINUE
C.              End SetVec.Vec
C:              Begin SetVec.Vec&Vec
                DO 78 L1=1,NM1
                  XU(L1)=X(L1+N)+HHA(L1)
78              CONTINUE
C.              End SetVec.Vec&Vec
C:              Begin SetVec.Vec
                DO 79 L1=1,NM1
                  HH(L1)=HHA(L1)
79              CONTINUE
C.              End SetVec.Vec
                IF(KPRINT.GE.0)THEN
80                FORMAT('0',5X,I2,1X,'Not ','accepted ',
     *            'relaxation ','factor',5X,F5.3)
                  WRITE(LUMON,80)ITER,FC
                ENDIF
                IF(ITER.EQ.0)THEN
                  FC = FCMIN
                ENDIF
                DIFAPP = .TRUE.
                JACRFR = .TRUE.
                GOTO 9998
              ENDIF
C             ----------------------------------------------------
C             4 Preparations to start the following iteration step
              ITER = ITER+1
              FCA = FC
C             ----------------------------------------------------
C             4.1 Print values
              IF(KPRINT.GE.0)THEN
C               Print Standard - and natural level, and damping
C               factor
                IF(KPRINT.GT.0)THEN
81                FORMAT('0',1X,66('*'))
                  WRITE(LUMON,81)
82                FORMAT('0',4X,'It',7X,'Levelf',10X,'Levelx',8X,
     *            'Rel.Fc.')
                  WRITE(LUMON,82)
                ENDIF
83              FORMAT('0',4X,I2,5X,D10.3,2X,4X,D10.3,2X,4X,F5.3)
                WRITE(LUMON,83)ITER,SUMF,SUMX,FC
                IF(KPRINT.GT.0)THEN
84                FORMAT('0',1X,66('*'))
                  WRITE(LUMON,84)
                  DO 85 J=1,M
86                  FORMAT(D13.5,2X)
                    WRITE(LUMON,86)T(J)
87                  FORMAT((14X,3(D20.10,1X)))
                    WRITE(LUMON,87)(X(L1),L1=(J-1)*N+1,J*N)
85                CONTINUE
                ENDIF
              ENDIF
9997        CONTINUE
C.          End of Segment BVPSOG.Core
9998      CONTINUE
        GOTO 17
        ENDIF
C.      EndWhile
C
C       End of main iteration loop
C       ==========================
C
C       ----------------------------------------------------------
C       5 Exits
C       ----------------------------------------------------------
C       5.1 Solution exit
        IF(INFO.EQ.0)THEN
          ITER = ITER+1
C:        Vec X = Vec X + Vec DXQ (for 1,NM)
          DO 88 L1=1,NM
            X(L1)=X(L1)+DXQ(L1)
88        CONTINUE
C.        End SetVec.&Vec
          INFO = ITER
          IF(KPRINT.LT.0)THEN
            GOTO 9999
          ENDIF
          IF(KPRINT.GT.0)THEN
C           Print levels, damping factor of last iteration step
C           and solution info
89          FORMAT('0',1X,66('*'))
            WRITE(LUMON,89)
90          FORMAT('0',4X,'It',7X,'Levelf',10X,'Levelx',8X,
     *      'Rel.Fc.')
            WRITE(LUMON,90)
          ENDIF
91        FORMAT('0',4X,I2,5X,D10.3,2X,4X,D10.3,2X,4X,F5.3)
          WRITE(LUMON,91)ITER,SUMF,SUMX,FC
92        FORMAT('0',1X,66('*'))
          WRITE(LUMON,92)
93        FORMAT('1')
          WRITE(LUMON,93)
C         Print solution info
94        FORMAT('0','Solution ','of',1X,'boundary ','value ',
     *    'problem',' obtained',/,'0','BVPSOL',' required',I3,1X,
     *    'Iteration ','steps ','with',I4,1X,'trajectory',
     *    ' evaluations',//)
          WRITE(LUMON,94)ITER,KOUNT
95        FORMAT('0','Achieved ','relative ','accuracy',D10.3,2X)
          WRITE(LUMON,95)CONV
          IF(EPH.GT.CONV) CONV = EPH
96        FORMAT('0','Reliable ','relative ','accuracy',D10.3,2X,/)
          WRITE(LUMON,96)CONV
          S = DLOG(DBLE((M-1)*(2*N+M-1))*EPMACH)
          DO 97 J=1,M1
            S = S+DLOG(RF(J))
97        CONTINUE
          IF(S.LT.DLOG(EPS))THEN
98          FORMAT('0','This ','boundary ','value ','problem ',
     *      'can ','also ','be ','solved ','with IOPT(3)=0 set',/)
            WRITE(LUMON,98)
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       5.2 Fail exit messages
C       ----------------------------------------------------------
C       5.2.1 Gaussian decomposition failed by singular jacobian
        IF(INFO.EQ.-1.AND.KPRINT.GE.0)THEN
99        FORMAT('0','Gaussian ','elimination ','failed ','by ',
     *    'singular ','Jacobian',/)
          WRITE(LUMON,99)
        ENDIF
C       ----------------------------------------------------------
C       5.2.2 Termination after more than itmax iterations
        IF(INFO.EQ.-2.AND.KPRINT.GE.0)THEN
100       FORMAT('0','Iteration ','terminates ','after ','itmax ',
     *    '=',I3,2X,'iteration ','steps')
          WRITE(LUMON,100)ITMAX
        ENDIF
C       ----------------------------------------------------------
C       5.2.3 Integrator failed in trajectory computation
        IF(INFO.EQ.-3.AND.KPRINT.GE.0)THEN
101       FORMAT('0','Integrator ','failed ','in ','trajectory ',
     *    'computation ',/)
          WRITE(LUMON,101)
          J1 =-KFLAG
102       FORMAT('0','BVPSOL ','terminates',/,'Subinterval',I3,1X,
     *    'possibly ','insert ','new ','node',D20.11,2X,/)
          WRITE(LUMON,102)J1,TFAIL
        ENDIF
C       ----------------------------------------------------------
C       5.2.4 Convergence fail of Gauss - Newton method
        IF(INFO.EQ.-4.AND.KPRINT.GE.0)THEN
103       FORMAT('0','Gauss ','Newton ','method ','fails ','to ',
     *    'converge',/)
          WRITE(LUMON,103)
        ENDIF
C       ----------------------------------------------------------
C       5.2.5 Inconsistent initial data
        IF(INFO.EQ.-5.AND.KPRINT.GE.0)THEN
104       FORMAT('0','Error: ','initial ','data ','and ',
     *    'boundary ','conditions ','are ','inconsistent',/)
          WRITE(LUMON,104)
        ENDIF
C       ----------------------------------------------------------
C       5.2.6 Multiple shooting condition too bad -insert new
C             nodes
        IF(INFO.EQ.-6.AND.KPRINT.GE.0)THEN
105       FORMAT('0','Termination ','since ','Multiple ',
     *    'Shooting ','condition',/,' ','or ','condition ','of ',
     *    'Jacobian ','is ','too ','bad',/,'insert ','new ',
     *    'nodes',/)
          WRITE(LUMON,105)
          S = REDH/TOL
          DO 106 J=1,M1
            IF(RF(J).GT.S)THEN
107           FORMAT('0',8X,'in ','subinterval',2X,I3,/)
              WRITE(LUMON,107)J
            ENDIF
106       CONTINUE
        ENDIF
C       ----------------------------------------------------------
C       5.2.7 Insufficient error tolerance for integrator
        IF(INFO.EQ.-7.AND.KPRINT.GE.0)THEN
          TOLH = DMIN1(REDH/COND,EPS)/SIGDEL
          RELDIF = DSQRT(TOLH/SIGDEL)
108       FORMAT('0','Suggested ','integrator ','accuracy',D10.1
     *    ,2X,/,'0','Suggested ','relative ','deviation ',
     *    'parameter',D10.1,2X,/)
          WRITE(LUMON,108)TOLH,RELDIF
109       FORMAT('0','Reduce ','relative ','error ','tolerance ',
     *    'for ','integrator ','to',D10.1,2X,/,2X,'or ','insert ',
     *    'new ','nodes',/)
          WRITE(LUMON,109)TOLH
          S = REDH/TOL
          DO 110 J=1,M1
            IF(RF(J).GT.S)THEN
111           FORMAT('0',8X,'in ','subinterval',2X,I3,/)
              WRITE(LUMON,111)J
            ENDIF
110       CONTINUE
112       FORMAT('0','Reliable ','relative ','accuracy ',
     *    'greater ','than',1X,D6.1,2X,/)
          WRITE(LUMON,112)1.0D-2
        ENDIF
C       ----------------------------------------------------------
C       5.2.8 Too small storage for sparse linear system solver
        IF(INFO.EQ.-9.AND.KPRINT.GE.0)THEN
113       FORMAT('0','Too ','small ','storage ','for ','linear ',
     *    'system ','solver',/)
          WRITE(LUMON,113)
        ENDIF
C       ----------------------------------------------------------
C       5.3 Common exit
        IF(KPRINT.GE.0)THEN
C
114       FORMAT('0','Condition ','of ','Jacobian',D10.3,2X,/,'0',
     *    'Multiple ','shooting ','condition',D10.3,2X,/,'1')
          WRITE(LUMON,114)COND,SIGDLH
          IF(INFO.GT.0)THEN
115         FORMAT('0','Solution ','data',/)
            WRITE(LUMON,115)
          ENDIF
          IF(INFO.LT.0)THEN
116         FORMAT('0','Final ','data',/)
            WRITE(LUMON,116)
          ENDIF
          DO 117 J=1,M
118         FORMAT(D13.5,2X)
            WRITE(LUMON,118)T(J)
119         FORMAT((14X,3(D20.10,1X)))
            WRITE(LUMON,119)(X(L1),L1=(J-1)*N+1,J*N)
117       CONTINUE
        ENDIF
C       End of exits
C       End of subroutine BVPSOG
9999  CONTINUE
C.    End of Segment BVPSOG.Body
      RETURN
      END
      SUBROUTINE BLFCNI(IVPSOL,FCN,BC,N,M,NM,NM1,ITER,KPRINT,
     *HSTART,FCMIN,T,X,X1,XM,T1,XU,HH,R,TOL,FC,FCOMPT,REDUCT,KFLAG,
     *KOUNT,INFO,LUMON)
      IMPLICIT DOUBLEPRECISION(S)
      EXTERNAL FCN,IVPSOL,BC
      INTEGER N,M,NM,NM1,ITER,KPRINT
      DOUBLE PRECISION HSTART,FCMIN
      INTEGER KOUNT,KFLAG,INFO,LUMON
      LOGICAL REDUCT
      DOUBLE PRECISION TOL,FC
      LOGICAL FCOMPT
      DOUBLE PRECISION T(M),X(NM)
      DOUBLE PRECISION XU(NM1),HH(NM1),R(N),X1(N),XM(N),T1(N)
C:    End Parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      INTEGER J,J1,KB,KB2
      DOUBLE PRECISION HMAX,HSAVE,TJ,TJ1,H
      INTEGER L1
C:    Begin
C:    Begin of Segment FcnInt.Body
C       Computation of the trajectories (solution of M1 initial
C         value problems)
        REDUCT = .FALSE.
        KOUNT = KOUNT+1
        HSAVE = HSTART
        DO 117 J=1,M-1
          J1 = J+1
          TJ = T(J)
          TJ1 = T(J1)
          H = HSAVE
          HMAX = DABS(TJ1-TJ)
          KFLAG = 0
          KB =(J-1)*N
C:        Begin SetVec.Vec
          DO 118 L1=1,N
            T1(L1)=X(L1+KB)
118       CONTINUE
C.        End SetVec.Vec
          CALL IVPSOL(N,FCN,TJ,T1,TJ1,TOL,HMAX,H,KFLAG)
          HSAVE = H
          IF(H.EQ.ZERO)THEN
C           trajectory computation failed
            IF(ITER.EQ.0)THEN
              INFO = -3
              GOTO 9993
            ENDIF
            IF(KPRINT.GE.0)THEN
119           FORMAT('0','trajectory ','computation ','failed, ',
     *        'relaxation ','factor ','or ','pseudo-rank ','reduced',/)
              WRITE(LUMON,119)
            ENDIF
            FC = FC*HALF
            IF(FC.LT.FCMIN)THEN
              REDUCT = .TRUE.
              GOTO 9993
            ENDIF
            FCOMPT = .FALSE.
            GOTO 9993
          ENDIF
          FCOMPT = .TRUE.
C         continuity conditions
C:        Begin SetVec.Vec
          DO 120 L1=1,N
            XU(L1+KB)=T1(L1)
120       CONTINUE
C.        End SetVec.Vec
          KB2 = KB+N
C:        Begin SetVec.Vec-Vec
          DO 121 L1=1,N
            HH(L1+KB)=T1(L1)-X(L1+KB2)
121       CONTINUE
C.        End SetVec.Vec-Vec
117     CONTINUE
C       two-point boundary conditions
C:      Begin SetVec.Vec
        DO 122 L1=1,N
          XM(L1)=X(L1+NM1)
122     CONTINUE
C.      End SetVec.Vec
C:      Begin SetVec.Vec
        DO 123 L1=1,N
          X1(L1)=X(L1)
123     CONTINUE
C.      End SetVec.Vec
        CALL BC(X1,XM,R)
9993  CONTINUE
C.    End of Segment FcnInt.Body
      RETURN
      END
      SUBROUTINE BLSOLI(N,M,M1,NM,NM1,LEVEL,NE,NB,IRANK,IRANKB,
     *IREPET,NYMAX,KPRINT,EPS,REDH,TOLMIN,TOL,RELDIF,EPH,EPX1H,
     *SIGDEL,SIGDLH,E,EH,HH,DHH,R,A,B,BG,G,QE,U,QU,DE,DU,T1,T2,US,
     *DX1,D,DDX,DXQ,XW,DR,RF,IROW,ICOL,ICOLB,PIVOT,NY,INFO,LUMON)
      IMPLICIT DOUBLEPRECISION(S)
      INTEGER N,M,M1,NM,NM1,LEVEL,NE,NB,IRANK,IRANKB,IREPET,NYMAX,
     *KPRINT
      DOUBLE PRECISION EPS,REDH,TOLMIN
      DOUBLE PRECISION TOL,RELDIF,EPH,EPX1H,SIGDEL,SIGDLH
      INTEGER INFO,NY,LUMON
      INTEGER IROW(N),ICOL(N),ICOLB(N),PIVOT(N)
      DOUBLE PRECISION G(N,N,M1)
      DOUBLE PRECISION A(N,N),B(N,N),BG(N,N),E(N,N),QE(N,N)
      DOUBLE PRECISION DDX(NM),DXQ(NM),XW(NM),HH(NM1),DHH(NM1),D(N),
     *DE(N),R(N),DR(N),U(N),DU(N),QU(N),T1(N),T2(N),DX1(N),RF(M),
     *US(N)
      DOUBLE PRECISION EH(N,N)
C:    End Parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      INTEGER I,IC,I0,J,JA,JIN,JN,J1,K,K0
      DOUBLE PRECISION CORR,EPDX1,EPX1,S,TH,TOLH
      LOGICAL NOTFND
      INTEGER L1,L2
      DOUBLE PRECISION S1
C:    Begin
C:    Begin of Segment SolvIn.Body
        IF(IREPET.GE.0)THEN
C         --------------------------------------------------------
C         1 computation of condensed right - hand side U(NE)
          IF(IRANK.GT.0) CALL BLRHS1(N,NE,M1,NM1,1,HH,R,B,G,U,DE,
     *    T1,BG,IROW)
C         --------------------------------------------------------
C         2 saving of right - hand side U
          IF(IRANK.GE.NE)THEN
C:          Begin SetVec.Vec
            DO 124 L1=1,NE
              US(L1)=U(L1)
124         CONTINUE
C.          End SetVec.Vec
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       3 ( best ) constrained least squares solution of linear(NE,
C         NE) -system
C:      Vec DX1 = Scalar (for 1,N)
        S1 = ZERO
        DO 125 L1=1,N
          DX1(L1)=S1
125     CONTINUE
C.      End SetVec.S
        IF(IRANK.GT.0)THEN
          CALL BLSOLC(E,N,N,IRANKB,NE,NE,DX1,U,IRANK,D,PIVOT,
     *    IREPET,QE,T1)
        ENDIF
C       ----------------------------------------------------------
C       4 iterative refinement of DX1
        EPH = EPS
        IF(IRANK.GE.NE.AND.NE.NE.0)THEN
C:        Begin SetVec.MatxVec
          DO 126 L1=1,NE
            S1=0.0
            DO 127 L2=1,NE
              S1=S1+EH(L1,L2)*DX1(L2)
127         CONTINUE
            DU(L1)=S1
126       CONTINUE
C.        End SetVec.MatxVec
C:        Begin SetVec.Vec-Vec
          DO 128 L1=1,NE
            DU(L1)=US(L1)-DU(L1)
128       CONTINUE
C.        End SetVec.Vec-Vec
C         Solution of residual equation
          CALL BLSOLC(E,N,N,IRANKB,NE,NE,T2,DU,IRANK,D,PIVOT,
     *    IREPET,QE,T1)
C:        EPDX1 = Max.Norm of Vec T2 (for 1,NE)
          EPDX1 = 0.0
          DO 129 L1=1,NE
            S1 = DABS(T2(L1))
            IF(S1.GT.EPDX1) EPDX1=S1
129       CONTINUE
C.        End MakeMaxNorm.Vec
C:        Vec DX1 = Vec DX1 + Vec T2 (for 1,NE)
          DO 130 L1=1,NE
            DX1(L1)=DX1(L1)+T2(L1)
130       CONTINUE
C.        End SetVec.&Vec
C:        EPX1 = Max.Norm of Vec DX1 (for 1,NE)
          EPX1 = 0.0
          DO 131 L1=1,NE
            S1 = DABS(DX1(L1))
            IF(S1.GT.EPX1) EPX1=S1
131       CONTINUE
C.        End MakeMaxNorm.Vec
          EPX1H = EPDX1/EPX1
          EPH = TEN*EPDX1
          IF(EPX1H.GT.HALF)THEN
            INFO = -8
            GOTO 9992
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       5 Descaling and back - permutation of solution DX1
C:      Permutation ICOL of Vec DXQ = Vec DX1 (for 1,N)
        DO 132 L1=1,N
          L2 = ICOL(L1)
          DXQ(L2) = DX1(L1)
132     CONTINUE
C.      End SetVecByPermVec.Vec
C:      Vec DXQ = Vec DXQ * Vec XW (for 1,N)
        DO 133 L1=1,N
          DXQ(L1)=DXQ(L1)*XW(L1)
133     CONTINUE
C.      End SetVec.VecxVec
C       ----------------------------------------------------------
C       6 Recursive computation of DXQ(N,2),  ... ,  DXQ(N,M)
        CALL BLRCRS(N,M,M1,NM,NM1,1,HH,G,DXQ,T1,T2)
C       ----------------------------------------------------------
C       1 Iterative refinement sweeps NY = 1 ,  ... ,  NYMAX
        NY = 0
        SIGDEL = TEN*TEN
        SIGDLH = ZERO
        IF(EPH.LT.EPS) EPH = EPS
        IF(NYMAX.NE.0.AND.IRANK.GE.NE.AND.NE.NE.0)THEN
          IF(KPRINT.GT.0)THEN
134         FORMAT('0','Iterative ','refinement',/)
            WRITE(LUMON,134)
          ENDIF
C         --------------------------------------------------------
C         1.1 Computation of required continuity residuals DHH(N,
C             M1)
          JN = 1
          JIN = M
C         --------------------------------------------------------
C         1.2 Computation of boundary residual DR(N)
C:        DO (Until)
135       CONTINUE
C:          Begin SetVec.MatxVec
            DO 136 L1=1,N
              S1=0.0
              DO 137 L2=1,N
                S1=S1+A(L1,L2)*DXQ(L2)
137           CONTINUE
              DR(L1)=S1
136         CONTINUE
C.          End SetVec.MatxVec
C:          Begin SetVec.MatxVec
            DO 138 L1=1,N
              S1=0.0
              DO 139 L2=1,N
                S1=S1+B(L1,L2)*DXQ(L2+NM1)
139           CONTINUE
              T1(L1)=S1
138         CONTINUE
C.          End SetVec.MatxVec
C:          Vec DR = Formula (for 1,N)
            DO 140 L1=1,N
              DR(L1)=R(L1)+DR(L1)+T1(L1)
140         CONTINUE
C.          End SetVec.Comp
C           ------------------------------------------------------
C           1.3 Computation of condensed residual DU(NE)
            IF(IRANK.GT.0) CALL BLRHS1(N,NE,M1,NM1,JIN,DHH,DR,B,G,
     *      DU,DE,T1,BG,IROW)
C           ------------------------------------------------------
C           1.4 Computation of correction DDX(N)
C:          Vec DX1 = Scalar (for 1,N)
            S1 = ZERO
            DO 141 L1=1,N
              DX1(L1)=S1
141         CONTINUE
C.          End SetVec.S
            IF(IRANK.GT.0) CALL BLSOLC(E,N,N,IRANKB,NE,NE,DX1,DU,
     *      IRANK,D,PIVOT,IREPET,QE,T1)
C           ------------------------------------------------------
C           2 Descaling of DDX(N),  refinement of DXQ(N)
C:          CORR = Max.Norm of Vec DX1 (for 1,N)
            CORR = 0.0
            DO 142 L1=1,N
              S1 = DABS(DX1(L1))
              IF(S1.GT.CORR) CORR=S1
142         CONTINUE
C.          End MakeMaxNorm.Vec
C:          Permutation ICOL of Vec T1 = Vec DX1 (for 1,N)
            DO 143 L1=1,N
              L2 = ICOL(L1)
              T1(L2) = DX1(L1)
143         CONTINUE
C.          End SetVecByPermVec.Vec
C:          Vec DDX = Vec T1 * Vec XW (for 1,N)
            DO 144 L1=1,N
              DDX(L1)=T1(L1)*XW(L1)
144         CONTINUE
C.          End SetVec.VecxVec
C:          Vec DXQ = Vec DXQ + Vec DDX (for 1,N)
            DO 145 L1=1,N
              DXQ(L1)=DXQ(L1)+DDX(L1)
145         CONTINUE
C.          End SetVec.&Vec
            IF(CORR.GE.EPH)THEN
              EPH = CORR
              INFO = -8
              GOTO 9992
            ENDIF
            RF(1)=CORR
C           ------------------------------------------------------
C           3 Recursive computation of DDX(N+1),  ... ,  DDX(NM)
            CALL BLRCRS(N,M,M1,NM,NM1,JIN,DHH,G,DDX,T1,T2)
C           ------------------------------------------------------
C           3.1 Refinement of DXQ(N+1),  ... ,  DXQ(NM)
            DO 146 J=2,M
              I0 =(J-1)*N
C:            Vec DXQ = Vec DXQ + Vec DDX (for I0+1,I0+N)
              DO 147 L1=I0+1,I0+N
                DXQ(L1)=DXQ(L1)+DDX(L1)
147           CONTINUE
C.            End SetVec.&Vec
C:            CORR = Max of Formula Elements (for I0+1,I0+N)
              CORR =-7.23D75
              DO 148 L1=I0+1,I0+N
                S1=DABS(DDX(L1)/XW(L1))
                IF(S1.GT.CORR)CORR=S1
148           CONTINUE
C.            End MakeMax.Comp
              RF(J)=CORR
146         CONTINUE
C           ------------------------------------------------------
C           3.2 Determination of sweep index JN
            JA = JN
            J = 1
            NOTFND = .TRUE.
C:          While (expression)
149         IF(J.LE.M.AND.NOTFND)THEN
              IF(RF(J).GT.EPH)THEN
                NOTFND = .FALSE.
              ELSE
                JN = J
                J = J+1
              ENDIF
            GOTO 149
            ENDIF
C.          EndWhile
            NY = NY+1
            IF(KPRINT.GT.0)THEN
150           FORMAT('0','Sweep',1X,I3,1X,'starts ','at',1X,I3)
              WRITE(LUMON,150)NY,JA
151           FORMAT((1X,5(D12.3,1X)))
              WRITE(LUMON,151)(RF(L1),L1=1,M)
            ENDIF
            IF(JN.LE.JA)THEN
              INFO = -6
              GOTO 9992
            ENDIF
            IF(JN.NE.M) JIN = JN
C           ------------------------------------------------------
C           3.3 Determination and adaptation of parameters TOL AND
C               RELDIF
            IF(LEVEL.NE.0.AND.NY.LE.1)THEN
              DO 152 J=1,M1
                S = RF(J+1)/RF(J)
                IF(SIGDLH.LT.S) SIGDLH = S
                RF(J)=S
152           CONTINUE
              IF(KPRINT.GT.0)THEN
153             FORMAT('0','Norms ','of ','wronskians')
                WRITE(LUMON,153)
154             FORMAT((1X,5(D12.3,1X)))
                WRITE(LUMON,154)(RF(L1),L1=1,M1)
              ENDIF
              SIGDEL = DMAX1(SIGDLH,SIGDEL)
              TH = TOL*SIGDEL
              IF(TH.GT.REDH)THEN
                INFO = -7
                GOTO 9992
              ENDIF
              IF(TH.GT.EPH) EPH = TH
              TOLH = EPS/SIGDEL
              IF(TOLH.LT.TOLMIN) TOLH = TOLMIN
CWei;         TOL = TOLH
CWei;         RELDIF = DSQRT(TOL/SIGDEL)
              IF(KPRINT.GE.0)THEN
155             FORMAT('0','Suggested ','integrator ','accuracy',D10.1
     *          ,2X,/,'0','Suggested ','relative ','deviation ',
     *          'parameter',D10.1,2X,//,'0','Adapted ','in ',
     *          'the ','next ','iteration ','step',/)
                WRITE(LUMON,155)TOLH,RELDIF
              ENDIF
            ENDIF
            IF(JN.NE.M)THEN
              DO 156 J=JN,M1
                J1 = J+1
                DO 157 I=1,N
                  K0 =(J-1)*N
                  S = HH(I+K0)
                  DO 158 K=1,N
                    S = S+G(I,K,J)*DXQ(K+K0)
158               CONTINUE
                  DHH(I+K0)=S-DXQ(I+K0+N)
157             CONTINUE
156           CONTINUE
            ENDIF
          IF(.NOT.(JN.EQ.M)) GOTO  135
C.        UNTIL ( expression - negated above)
        ENDIF
C       End of iterative refinement sweeps
C       ----------------------------------------------------------
C       4 Projection of separated linear boundary conditions at T(
C         M)
        IF(NB.NE.0)THEN
          DO 159 K=1,NB
            IC = ICOLB(K)
            DXQ(IC+NM1)=ZERO
159       CONTINUE
        ENDIF
9992  CONTINUE
C.    End of Segment SolvIn.Body
      RETURN
      END
      SUBROUTINE BGSOLI(N,M,M1,NM,NM1,NDIM,LICN,NKEEP,NA,NAQ,NB,
     *NBQ,NRS,ITER,LEVEL,KPRINT,EPS,REDH,TOLMIN,FC,FCA,TOL,RELDIF,
     *EPH,EPX1H,SIGDEL,SIGDLH,COND,CORR,HH,DHH,R,A,B,G,U,DE,DU,T1,
     *DXQ,XW,DR,RF,WO,E,IROW,ICOLA,ICOLB,ICN,IKEEP,INFO,LUMON)
      IMPLICIT DOUBLEPRECISION(S)
      INTEGER N,M,M1,NM,NM1,NDIM,LICN,NKEEP,NA,NAQ,NB,NBQ,NRS,ITER,
     *LEVEL,KPRINT,LUMON
      DOUBLE PRECISION EPS,REDH,TOLMIN,FC,FCA
      DOUBLE PRECISION TOL,RELDIF,EPH,EPX1H,SIGDEL,SIGDLH,COND,
     *CORR
      INTEGER INFO
      INTEGER IROW(N),ICOLA(N),ICOLB(N)
      INTEGER ICN(LICN),IKEEP(NKEEP)
      DOUBLE PRECISION G(N,N,M1)
      DOUBLE PRECISION A(N,N),B(N,N)
      DOUBLE PRECISION DXQ(NM),XW(NM),HH(NM1),DHH(NM1),DE(N),R(N),
     *DR(N),U(NM),DU(NM),T1(N),RF(M),E(LICN),WO(NM)
C:    End Parameter
C:    EPMACH = relative machine precision
      DOUBLE PRECISION EPMACH,SMALL
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      INTEGER I,IC,I0,J,J1,K,MTYPE
      DOUBLE PRECISION S,TOLH
      REAL ST
      INTEGER L1,L2
      DOUBLE PRECISION S1
C:    Begin
C:    Begin of Segment SolvIn.Body
C       ----------------------------------------------------------
C       1 Computation of right - hand side U(NDIM)
        CALL ZIBCONST(EPMACH,SMALL)
        DO 127 J=1,M1
          J0 =(J-1)*N
          J1 = J0+N
          DO 128 I=1,N
            U(I+J0)=HH(I+J0)/XW(I+J1)
128       CONTINUE
127     CONTINUE
        IF(NRS.NE.0)THEN
          DO 129 I=1,NRS
            I0 = IROW(NB+I)
            U(NM1+I)=DE(I0)*R(I0)
129       CONTINUE
        ENDIF
C       ----------------------------------------------------------
C       2 Solution of linear(NDIM,NDIM) -system
        MTYPE = 1
        CALL MA28CD(NDIM,E,LICN,ICN,IKEEP,U,WO,MTYPE)
        COND = ZERO
        IF(NAQ.NE.0)THEN
          DO 130 I=1,NAQ
            IC = ICOLA(I)
            S = U(I)
            DXQ(IC)=S*XW(IC)
            COND = COND+DABS(S)
130       CONTINUE
        ENDIF
        IF(M.NE.2)THEN
          DO 131 J=2,M1
            J0 =(J-1)*N
            DO 132 I=1,N
              S = U(J0-NA+I)
              DXQ(I+J0)=S*XW(I+J0)
              COND = COND+DABS(S)
132         CONTINUE
131       CONTINUE
        ENDIF
        IF(NBQ.NE.0)THEN
          DO 133 I=1,NBQ
            IC = ICOLB(NB+I)
            S = U(NM1-NA+I)
            DXQ(IC+NM1)=S*XW(IC+NM1)
            COND = COND+DABS(S)
133       CONTINUE
        ENDIF
        IF(NB.NE.0)THEN
          DO 134 K=1,NB
            IC = ICOLB(K)
            DXQ(IC+NM1)=ZERO
134       CONTINUE
        ENDIF
        IF(NA.NE.0)THEN
          DO 135 K=1,NA
            IC = ICOLA(NAQ+K)
            DXQ(IC)=ZERO
135       CONTINUE
        ENDIF
C       ----------------------------------------------------------
C       3 Iterative refinement
        IF(KPRINT.GT.0)THEN
136       FORMAT('0','Iterative ','refinement',/)
          WRITE(LUMON,136)
        ENDIF
C       ----------------------------------------------------------
C       4 Computation of required continuity residuals DHH(NM1)
        DO 137 J=1,M1
          J0 =(J-1)*N
          J1 = J0+N
          DO 138 I=1,N
            S = HH(I+J0)
            DO 139 K=1,N
              S = S+G(I,K,J)*DXQ(K+J0)
139         CONTINUE
            DHH(I+J0)=S-DXQ(I+J1)
138       CONTINUE
137     CONTINUE
C       ----------------------------------------------------------
C       5 Computation of boundary residual DR(N)
C:      Begin SetVec.MatxVec
        DO 140 L1=1,N
          S1=0.0
          DO 141 L2=1,N
            S1=S1+A(L1,L2)*DXQ(L2)
141       CONTINUE
          DR(L1)=S1
140     CONTINUE
C.      End SetVec.MatxVec
C:      Begin SetVec.MatxVec
        DO 142 L1=1,N
          S1=0.0
          DO 143 L2=1,N
            S1=S1+B(L1,L2)*DXQ(L2+NM1)
143       CONTINUE
          T1(L1)=S1
142     CONTINUE
C.      End SetVec.MatxVec
C:      Vec DR = Formula (for 1,N)
        DO 144 L1=1,N
          DR(L1)=R(L1)+DR(L1)+T1(L1)
144     CONTINUE
C.      End SetVec.Comp
C       ----------------------------------------------------------
C       6 Computation of residual DU(NDIM)
        DO 145 J=1,M1
          J0 =(J-1)*N
          J1 = J0+N
          DO 146 I=1,N
            DU(J0+I)=DHH(I+J0)/XW(I+J1)
146       CONTINUE
145     CONTINUE
        IF(NRS.NE.0)THEN
          DO 147 I=1,NRS
            I0 = IROW(NB+I)
            DU(NM1+I)=DE(I0)*DR(I0)
147       CONTINUE
        ENDIF
C       ----------------------------------------------------------
C       7 Computation of correction
        MTYPE = 1
        CALL MA28CD(NDIM,E,LICN,ICN,IKEEP,DU,WO,MTYPE)
C       ----------------------------------------------------------
C       8 Refinement of DXQ
        CORR = ZERO
        IF(NAQ.NE.0)THEN
          DO 148 I=1,NAQ
            IC = ICOLA(I)
            S = DU(I)
            DXQ(IC)=DXQ(IC)+S*XW(IC)
            CORR = CORR+DABS(S)
148       CONTINUE
        ENDIF
        IF(M.NE.2)THEN
          DO 149 J=2,M1
            J0 =(J-1)*N
            DO 150 I=1,N
              S = DU(J0-NA+I)
              DXQ(I+J0)=DXQ(I+J0)+S*XW(I+J0)
              CORR = CORR+DABS(S)
150         CONTINUE
149       CONTINUE
        ENDIF
        IF(NBQ.NE.0)THEN
          DO 151 I=1,NBQ
            IC = ICOLB(NB+I)
            S = DU(NM1-NA+I)
            DXQ(IC+NM1)=DXQ(IC+NM1)+S*XW(IC+NM1)
            CORR = CORR+DABS(S)
151       CONTINUE
        ENDIF
        COND = CORR/(COND*EPMACH)
        IF(KPRINT.GT.0)THEN
152       FORMAT('0','Norm ','of ','residual',D12.3,2X)
          WRITE(LUMON,152)CORR
        ENDIF
C       End of iterative refinement
C       ----------------------------------------------------------
C       9 Determination and adaptation of parameters TOL and
C         RELDIF
        IF(LEVEL.NE.0)THEN
          DO 153 J=1,M1
            J0 =(J-1)*N
            J1 = J0+N
            SIGDEL = ZERO
            ST = ZERO
            DO 154 I=1,N
              II = N
              IF(J.EQ.1) II = NAQ
              S = DABS(DXQ(I+J0))/XW(I+J0)
              IF(ST.LT.S) ST = S
              S = ZERO
              DO 155 K=1,II
                KK = K
                IF(J.EQ.1) KK = ICOLA(KK)
                S = S+G(I,KK,J)*DXQ(KK+J0)/XW(I+J1)
155           CONTINUE
              S = DABS(S)
              IF(SIGDEL.LT.S) SIGDEL = S
154         CONTINUE
            RF(J)=SIGDEL/ST+ONE
153       CONTINUE
          IF(KPRINT.GT.0)THEN
156         FORMAT('0','Norms ','of ','wronskians')
            WRITE(LUMON,156)
157         FORMAT((1X,5(D12.3,1X)))
            WRITE(LUMON,157)(RF(L1),L1=1,M1)
          ENDIF
          SIGDLH = ZERO
          DO 158 J=1,M1
            IF(SIGDLH.LT.RF(J)) SIGDLH = RF(J)
158       CONTINUE
          SIGDEL = SIGDLH
          IF(FC.EQ.ONE.AND.FCA.EQ.ONE.AND.ITER.GT.0) SIGDEL =
     *    SIGDLH*COND
          SIGDEL = DMAX1(SIGDEL,TEN)
          EPH = TOL*SIGDEL
          IF(EPH.GT.REDH)THEN
            INFO = -7
            GOTO 9989
          ENDIF
          TOLH = EPS/SIGDEL
          IF(TOLH.LT.TOLMIN) TOLH = TOLMIN
CWEI;     TOL = TOLH
CWEI;     RELDIF = DSQRT(TOL/SIGDEL)
          IF(KPRINT.GE.0)THEN
159         FORMAT('0','Suggested ','integrator ','accuracy',D10.1
     *      ,2X,/,'0','Suggested ','relative ','deviation ',
     *      'parameter',D10.1,2X,/)
            WRITE(LUMON,159)TOLH,RELDIF
160         FORMAT('0','Adapted ','in ','the ','next ',
     *      'iteration ','step',/)
            WRITE(LUMON,160)
          ENDIF
        ENDIF
9989  CONTINUE
C.    End of Segment SolvIn.Body
      RETURN
      END
      SUBROUTINE BLPRCD(LUMON,COND,SENS,SMALIN,J,IRANK)
      IMPLICIT DOUBLEPRECISION(S)
      DOUBLE PRECISION COND,SENS,SMALIN
      INTEGER J,IRANK,LUMON
C:    End Parameter
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION SENSP
C:    Begin
      IF(SENS.LT.ONE)THEN
        SENSP = SENS*SMALIN
160     FORMAT('0','Subcondition (',I2,',',I2,') ',D10.3,2X,/,'0',
     *  'Sensitivity  (',I2,',',I2,') ',D10.3,2X,/)
        WRITE(LUMON,160)J,IRANK,COND,J,IRANK,SENSP
      ELSE
161     FORMAT('0','Subcondition ','(',I2,',',I2,') ',D10.3,2X,/,
     *  '0','Sensitivity ','(',I2,',',I2,') ',D10.3,2X,' *',D7.0
     *  ,2X,/)
        WRITE(LUMON,161)J,IRANK,COND,J,IRANK,SENS,SMALIN
      ENDIF
      END
      SUBROUTINE BLPRCV(LUMON,CONV,EPH)
      IMPLICIT DOUBLEPRECISION(S)
      DOUBLE PRECISION CONV,EPH
      INTEGER LUMON
C:    End Parameter
C:    Begin
162   FORMAT('0','Achieved ','relative ','accuracy',D10.3,2X)
      WRITE(LUMON,162)CONV
      IF(EPH.GT.CONV) CONV = EPH
163   FORMAT('0','Reliable ','relative ','accuracy',D10.3,2X,/)
      WRITE(LUMON,163)CONV
      END
      SUBROUTINE BLSCLE(N,M,NM,NM1,X,XU,XW,XTHR)
      IMPLICIT DOUBLEPRECISION(S)
      INTEGER N
      INTEGER M
      INTEGER NM
      INTEGER NM1
      DOUBLE PRECISION X(NM)
      DOUBLE PRECISION XW(NM)
      DOUBLE PRECISION XU(NM1)
      DOUBLE PRECISION XTHR
C:    End Parameter
C     PROVIDES SCALING XW(NM)OF VARIABLES X(NM)
C:    EPMACH = relative machine precision
      DOUBLE PRECISION EPMACH,SMALL
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      DOUBLE PRECISION RED
      PARAMETER (RED=1.0D-2)
      INTEGER I,J,J0,J1
      DOUBLE PRECISION XMAX
      INTEGER L1
C:    Begin
C:    Vec XW = Formula (for 1,N)
      CALL ZIBCONST(EPMACH,SMALL)
      DO 164 L1=1,N
        XW(L1)=DABS(X(L1))
164   CONTINUE
C.    End SetVec.Comp
C     ------------------------------------------------------------
C     1 Arithmetic mean for XW(2*N)... XW(M*N)
      DO 165 J=1,M-1
        J0 =(J-1)*N
        J1 = J0+N
        DO 166 I=1,N
          XW(I+J1)=(DABS(X(I+J1))+DABS(XU(I+J0)))*HALF
166     CONTINUE
165   CONTINUE
C     ------------------------------------------------------------
C     2 Threshold
      DO 167 I=1,N
        XMAX = ZERO
        DO 168 J=0,NM1,N
          IF(XMAX.LT.XW(I+J)) XMAX = XW(I+J)
168     CONTINUE
        XMAX = XMAX*RED
        IF(XMAX.LT.XTHR) XMAX = XTHR
        DO 169 J=0,NM1,N
          IF(XW(I+J).LT.XMAX) XW(I+J)=XMAX
169     CONTINUE
167   CONTINUE
      RETURN
C     End of subroutine BLSCLE
      END
      SUBROUTINE BLLVLS(N,M,NM,NM1,XW,DXQ,HH,R,DE,CONV,SUMX,SUMF,
     *KPRINT)
      IMPLICIT DOUBLEPRECISION(S)
C
      INTEGER N,M,NM,NM1,KPRINT
      DOUBLE PRECISION XW(NM),DXQ(NM),HH(NM1),R(N),DE(N)
      DOUBLE PRECISION CONV,SUMX,SUMF
C:    End Parameter
C:    SMALL = squareroot of "smallest positive machine number
C     divided by relative machine precision"
      DOUBLE PRECISION SMALL
      PARAMETER (SMALL=4.94D-32)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER I,J,J0,J1
      DOUBLE PRECISION S
      INTEGER L1
C:    Begin
C     ------------------------------------------------------------
C     1 Evaluation of scaled natural level function SUMX and
C       scaled maximum error norm CONV
      CONV = ZERO
      SUMX = ZERO
      DO 170 J=1,NM
        S = DABS(DXQ(J)/XW(J))
        IF(CONV.LT.S) CONV = S
        SUMX = SUMX+S*S
170   CONTINUE
C     ------------------------------------------------------------
C     2 Evaluation of (scaled) standard level function sumfs (only
C       if needed for print)
C:    SUMF = Sum of Formula Elements (for 1,N)
      SUMF = 0.0D0
      DO 171 L1=1,N
        SUMF=SUMF+((R(L1)*DE(L1)/SMALL)**2)
171   CONTINUE
C.    End MakeSum.Comp
      DO 172 J=1,M-1
        J0 =(J-1)*N
        J1 = J0+N
        DO 173 I=1,N
          SUMF = SUMF+(HH(I+J0)/XW(I+J1))**2
173     CONTINUE
172   CONTINUE
C     End of subroutine BLLVLS
      RETURN
      END
      SUBROUTINE BGLVLS(N,M,NM,NM1,XW,DXQ,HH,R,DE,CONV,SUMX,SUMF,
     *KPRINT)
      IMPLICIT DOUBLEPRECISION(S)
C
      INTEGER N,M,NM,NM1,KPRINT
      DOUBLE PRECISION XW(NM),DXQ(NM),HH(NM1),R(N),DE(N)
      DOUBLE PRECISION CONV,SUMX,SUMF
C:    End Parameter
C:    SMALL = squareroot of "smallest positive machine number
C     divided by relative machine precision"
      DOUBLE PRECISION SMALL
      PARAMETER (SMALL=4.94D-32)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER I,J,J0,J1
      DOUBLE PRECISION S
      INTEGER L1
C:    Begin
C     ------------------------------------------------------------
C     1 Evaluation of scaled natural level function SUMX and
C       scaled maximum error norm CONV
      CONV = ZERO
      SUMX = ZERO
      DO 167 J=1,NM
        S = DABS(DXQ(J)/XW(J))
        IF(CONV.LT.S) CONV = S
        SUMX = SUMX+S*S
167   CONTINUE
C     ------------------------------------------------------------
C     2 Evaluation of (scaled) standard level function sumfs (only
C       if needed for print)
C:    SUMF = Sum of Formula Elements (for 1,N)
      SUMF = 0.0D0
      DO 168 L1=1,N
        SUMF=SUMF+((R(L1)*DE(L1))**2)
168   CONTINUE
C.    End MakeSum.Comp
      DO 169 J=1,M-1
        J0 =(J-1)*N
        J1 = J0+N
        DO 170 I=1,N
          SUMF = SUMF+(HH(I+J0)/XW(I+J1))**2
170     CONTINUE
169   CONTINUE
C     End of subroutine BGLVLS
      RETURN
      END
      SUBROUTINE BLRHS1(N,NE,M1,NM1,JIN,HH,R,B,G,U,DE,V,BG,IROW)
      IMPLICIT DOUBLEPRECISION(S)
      INTEGER N,NE,M1,NM1,JIN
      DOUBLE PRECISION HH(NM1),R(N)
      DOUBLE PRECISION B(N,N)
      DOUBLE PRECISION G(N,N,M1)
      DOUBLE PRECISION U(N),DE(N),V(N)
      DOUBLE PRECISION BG(N,N)
      INTEGER IROW(N)
C:    End Parameter
C     Computation of condensed right-hand side U(NE)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER I,IR,J,JJ,J1,K,K0,L,M2
      DOUBLE PRECISION S,TH
C:    Begin
      DO 174 I=1,NE
        IR = IROW(I)
        U(I)=DE(IR)*R(IR)
174   CONTINUE
      IF(JIN.GT.M1)THEN
        RETURN
      ENDIF
      DO 175 I=1,NE
        IR = IROW(I)
        S = U(I)
        K0 = NM1-N
        DO 176 K=1,N
          TH = DE(IR)*B(IR,K)
          BG(I,K)=TH
          S = S+TH*HH(K+K0)
176     CONTINUE
        U(I)=S
175   CONTINUE
      IF(M1.EQ.1.OR.JIN.EQ.M1) RETURN
      M2 = M1-1
      DO 177 JJ=JIN,M2
        J = M2+JIN-JJ
        J1 = J+1
        DO 178 I=1,NE
          DO 179 K=1,N
            S = ZERO
            DO 180 L=1,N
              S = S+BG(I,L)*G(L,K,J1)
180         CONTINUE
            V(K)=S
179       CONTINUE
          S = U(I)
          K0 =(J-1)*N
          DO 181 K=1,N
            S = S+V(K)*HH(K+K0)
            BG(I,K)=V(K)
181       CONTINUE
          U(I)=S
178     CONTINUE
177   CONTINUE
C     End of subroutine BLRHS1
      RETURN
      END
      SUBROUTINE BLRCRS(N,M,M1,NM,NM1,JIN,HH,G,DX,U,V)
      IMPLICIT DOUBLEPRECISION(S)
      INTEGER N,M,M1,NM,NM1,JIN
      DOUBLE PRECISION HH(NM1)
      DOUBLE PRECISION G(N,N,M1)
      DOUBLE PRECISION DX(NM),U(N),V(N)
C:    End Parameter
C     Recursive solution of m1 linear(N,N)-systems
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER I,J,J0,J1,K
      INTEGER L1
      DOUBLE PRECISION S
C:    Begin
C:    Begin SetVec.Vec
      DO 182 L1=1,N
        U(L1)=DX(L1)
182   CONTINUE
C.    End SetVec.Vec
      DO 183 J=1,M1
        J0 =(J-1)*N
        J1 = J0+N
        DO 184 I=1,N
          IF(J.GE.JIN)THEN
            S = HH(I+J0)
          ELSE
            S = ZERO
          ENDIF
          DO 185 K=1,N
            S = S+G(I,K,J)*U(K)
185       CONTINUE
          V(I)=S
          DX(I+J1)=S
184     CONTINUE
C:      Begin SetVec.Vec
        DO 186 L1=1,N
          U(L1)=V(L1)
186     CONTINUE
C.      End SetVec.Vec
183   CONTINUE
C     End of subroutine BLRCRS
      RETURN
      END
      SUBROUTINE BLPRJC(N,NE,IRANK,DEL,U,D,V,QE,PIVOT)
      IMPLICIT DOUBLEPRECISION(S)
C
      INTEGER IRANK,N,NE
      INTEGER PIVOT(N)
      DOUBLE PRECISION DEL
      DOUBLE PRECISION QE(N,N)
      DOUBLE PRECISION U(N),D(N),V(N)
C:    End Parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1,L2
      INTEGER I,IRK1
      DOUBLE PRECISION S,SH
C:    Begin
      DO 187 I=1,NE
        V(I)=U(PIVOT(I))
187   CONTINUE
      IRK1 = IRANK+1
      DEL = ZERO
      DO 188 I=IRK1,NE
C:      SH = Col I of QE * Vec V (for 1,I-1)
        SH = 0.0
        DO 189 L1=1,I-1
          SH = SH+QE(L1,I)*V(L1)
189     CONTINUE
C.      End MakeSProd.ColxVec
        S =(V(I)-SH)/D(I)
        DEL = DEL-S*S
        V(I)=S
188   CONTINUE
      DO 190 I=IRK1,NE
        K = NE+IRK1-I
        S = V(K)
        IF(K.NE.NE)THEN
C:        SH = Row K of QE * Vec V (for K+1,NE)
          SH = 0.0
          DO 191 L1=K+1,NE
            SH=SH+QE(K,L1)*V(L1)
191       CONTINUE
C.        End MakeSProd.RowxVec
          S = S-SH
        ENDIF
        S = S/D(K)
        V(K)=S
190   CONTINUE
      DO 192 I=1,IRANK
C:      S = Row I of QE * Vec V (for IRK1,NE)
        S = 0.0
        DO 193 L1=IRK1,NE
          S=S+QE(I,L1)*V(L1)
193     CONTINUE
C.      End MakeSProd.RowxVec
        V(I)=-S
192   CONTINUE
C:    Permutation PIVOT of Vec U = Vec V (for 1,NE)
      DO 194 L1=1,NE
        L2 = PIVOT(L1)
        U(L2) = V(L1)
194   CONTINUE
C.    End SetVecByPermVec.Vec
C     End of subroutine BLPRJC
      RETURN
      END
      SUBROUTINE BLDERA(BC,N,M,NM,XW,X1,XM,R,RH,A,B,RELDIF)
      IMPLICIT DOUBLEPRECISION(S)
      EXTERNAL BC
      INTEGER N,M,NM
      DOUBLE PRECISION XW(NM),X1(N),XM(N),R(N),RH(N)
      DOUBLE PRECISION A(N,N),B(N,N)
      DOUBLE PRECISION RELDIF
C:    End Parameter
C     Difference approx. of boundary derivative matrices A(N,N)and
C       B(N,N)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER I,K,NM1
      DOUBLE PRECISION S,XH
C:    Begin
      NM1 = N*(M-1)
      DO 195 K=1,N
        XH = X1(K)
        S = RELDIF*XW(K)
        IF(XH.LT.ZERO) S =-S
        X1(K)=XH+S
        CALL BC(X1,XM,RH)
        X1(K)=XH
        S = ONE/S
        DO 196 I=1,N
          A(I,K)=(RH(I)-R(I))*S
196     CONTINUE
        XH = XM(K)
        S = RELDIF*XW(K+NM1)
        IF(XH.LT.ZERO) S =-S
        XM(K)=XH+S
        CALL BC(X1,XM,RH)
        XM(K)=XH
        S = ONE/S
        DO 197 I=1,N
          B(I,K)=(RH(I)-R(I))*S
197     CONTINUE
195   CONTINUE
C     END SUBROUTINE BLDERA
      RETURN
      END
      SUBROUTINE BLDERG(FCN,N,NE,M,M1,NM,NM1,T,X,XU,XW,XJ,TJ,G,
     *ICOL,IVPSOL,HSTART,TOL,RELDIF,KFLAG)
      IMPLICIT DOUBLEPRECISION(S)
      INTEGER N
      INTEGER NE
      INTEGER M
      INTEGER M1
      INTEGER NM
      INTEGER NM1
      DOUBLE PRECISION T(M)
      DOUBLE PRECISION X(NM)
      DOUBLE PRECISION XU(NM1)
      DOUBLE PRECISION XW(NM)
      DOUBLE PRECISION XJ(N)
      DOUBLE PRECISION TJ
      DOUBLE PRECISION G(N,N,M1)
      INTEGER ICOL(N)
      EXTERNAL IVPSOL
      DOUBLE PRECISION HSTART
      DOUBLE PRECISION TOL
      DOUBLE PRECISION RELDIF
      INTEGER KFLAG
C:    End Parameter
C     Difference approximation of Wronskian Matrices G(1),.., G(M1)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER I,IK,J,J0,J1,K
      DOUBLE PRECISION HMAX,H,HSAVE,S,TJ1,TJA,TH
      EXTERNAL FCN
C:    Begin
      HSAVE = HSTART
      DO 198 J=1,M-1
        J0 =(J-1)*N
        J1 = J+1
        TJA = T(J)
        TJ1 = T(J1)
        HMAX = DABS(TJ1-TJA)
        DO 199 IK=1,N
          I = ICOL(IK)
          H = HSAVE
          IF(J.NE.1.OR.IK.LE.NE)THEN
            TJ = TJA
            KFLAG = 0
            DO 200 K=1,N
              XJ(K)=X(K+J0)
200         CONTINUE
            TH = XJ(I)
            S = RELDIF*XW(I+J0)
            IF(TH.LT.ZERO) S =-S
            XJ(I)=TH+S
            S = ONE/S
            CALL IVPSOL(N,FCN,TJ,XJ,TJ1,TOL,HMAX,H,KFLAG)
            IF(H.EQ.ZERO)THEN
              KFLAG =-J
              RETURN
            ENDIF
            DO 201 K=1,N
              G(K,I,J)=S*(XJ(K)-XU(K+J0))
201         CONTINUE
          ENDIF
199     CONTINUE
        HSAVE = H
198   CONTINUE
      KFLAG = 0
C     END OF SUBROUTINE BLDERG
      RETURN
      END
      SUBROUTINE BLRK1G(N,M,M1,NM,NM1,XW,DX,HH,HHA,DXJ,G,FCA)
      IMPLICIT DOUBLEPRECISION(S)
      INTEGER N,M,M1,NM,NM1
      DOUBLE PRECISION XW(NM),DX(NM),HH(NM1),HHA(NM1),DXJ(N)
      DOUBLE PRECISION G(N,N,M1)
      DOUBLE PRECISION FCA
C:    End Parameter
C     RANK-1 UPDATES OF WRONSKIAN MATRICES G(1),..., G(M1)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER I,J,J0,K
      DOUBLE PRECISION DNM,FCH,S,T
C:    Begin
      FCH = FCA-ONE
      DO 202 J=1,M1
        J0 =(J-1)*N
        DNM = ZERO
        DO 203 I=1,N
          T = DX(I+J0)/XW(I+J0)
          DXJ(I)=T/XW(I+J0)
          DNM = DNM+T*T
203     CONTINUE
        DNM = DNM*FCA
        IF(DNM.NE.ZERO)THEN
          DO 204 K=1,N
            T = DXJ(K)/DNM
            DO 205 I=1,N
              S = G(I,K,J)
              IF(S.NE.ZERO) G(I,K,J)=S+T*(HH(I+J0)+FCH*HHA(I+J0))
205         CONTINUE
204       CONTINUE
        ENDIF
202   CONTINUE
C     END OF SUBROUTINE BLRK1G
      RETURN
      END
C
C*    Group  Initial value problem solver (Code DIFEX1)   
C
      SUBROUTINE BLDFX1 (N,FCN,T,Y,TEND,TOL,HMAX,H,KFLAG)
C
C* Begin prologue DIFEX1
C
C  ---------------------------------------------------------------------
C
C* Title
C
C    Explicit extrapolation integrator for non-stiff systems of
C    ordinary first-order differential equations.
C
C* Written by        P. Deuflhard, U. Nowak, U. Poehle
C                    Adapted by L. Weimann for use with BVPSOL
C* Purpose           Solution of systems of initial value problems
C* Method            Explicit mid-point rule discretization with
C                    h**2-extrapolation
C* Category          i1a1c1. - System of nonstiff first order
C                              differential equations
C* Keywords          extrapolation, ODE, explicit mid-point rule,
C                    nonstiff
C* Version           1.2 , August 1991
C* Latest Change     January 2004
C* Library           CodeLib
C* Code              Fortran 77
C                    Double Precision
C* Environment       Standard version for FORTRAN77 environments on
C                    PCs, workstations, and hosts
C* Copyright     (c) Konrad-Zuse-Zentrum fuer Informationstechnik
C                    Berlin (ZIB)
C                    Takustrasse 7, D-14195 Berlin-Dahlem
C                    phone : + 49/30/84185-0
C                    fax   : + 49/30/84185-125
C* Contact           Rainald Ehrig
C                    ZIB, Numerical Analysis and Modelling
C                    phone : + 49/30/84185-282
C                    fax   : + 49/30/84185-107
C                    e-mail: ehrig@zib.de
C
C  ---------------------------------------------------------------------
C
C* Licence
C  -------
C
C  You may use or modify this code for your own non-commercial
C  purposes for an unlimited time. 
C  In any case you should not deliver this code without a special 
C  permission of ZIB.
C  In case you intend to use the code commercially, we oblige you
C  to sign an according licence agreement with ZIB.
C
C
C* Warranty
C  --------
C 
C  This code has been tested up to a certain level. Defects and
C  weaknesses, which may be included in the code, do not establish
C  any warranties by ZIB. ZIB does not take over any liabilities
C  which may follow from acquisition or application of this code.
C
C
C* Software status 
C  ---------------
C
C  This code is under care of ZIB and belongs to ZIB software
C  class II.
C
C
C  ---------------------------------------------------------------------
C
C* References:
C
C /1/ W. B. Gragg:
C     On Extrapolation Algorithms for Ordinary
C     Initial Value Problems
C     SIAM J. Numer. Anal. 2, 384-404 (1965)
C
C /2/ R. Bulirsch, J. Stoer:
C     Numerical Treatment of Ordinary Differential Equations
C     by Extrapolation Methods
C     Num. Math. 8, 1-13 (1966)
C
C /3/ P. Deuflhard:
C     Order and Stepsize Control in Extrapolation Methods
C     Numer. Math. 41, 399-422 (1983)
C
C
C* External subroutine: (to be supplied by the user)
C
C    FCN           EXT  Subroutine FCN (N,T,Y,DY)
C                       Right-hand side of first-order
C                       differential equations
C                       N      Number of first-order ODEs
C                       T      Actual position
C                       Y(N)   Values at T
C                       DY(N)  Derivatives at T
C
C
C* Parameters: (* marks input/output parameters)
C
C    N         I   IN   Number of first-order ODEs
C  * T         D   IN   Starting point of integration
C                  OUT  Achieved final point of integration
C  * Y         D   IN   Array of initial values Y(1),...,Y(N)
C                  OUT  Array of final values
C    TEND      D   IN   Prescribed final point of integration
C    TOL       D   IN   Prescribed relative precision (.GT.0)
C    HMAX      D   IN   Maximum permitted stepsize
C  * H         D   IN   Initial stepsize guess
C                  OUT  Stepsize proposal for next integration step
C                       (H .EQ. 0. ,if DIFEX1 fails to proceed)
C  * KFLAG     I   IN   Print parameter
C                         0   no output
C                         1   Integration monitor
C                         2   Intermediate solution points  T,Y(I),I=1,N
C                         3   Integration monitor and solution points
C                       +10   if KFLAG is augmented by 10, a time
C                             monitor is printed additionally (this may
C                             be very expensive in terms of cpu-time)
C                  OUT  Error flag
C                       .GE. 0  Successful integration
C                               (KFLAG not altered internally)
C                       -2   More than NSTMAX basic integration steps
C                            per interval have been performed
C                       -3   More than JRMAX stepsize reductions
C                            occurred per basic integration step
C                       -4   Stepsize proposal for next basic
C                            integration step too small
C
C
C* End prologue
C  ------------
C
C
C    COMMON /STAT/ NFCN, NSTEP, NACCPT, NREJCT, NDEC, NSOL
C                       Internally initialized, for statistical
C                       purposes
C    NFCN               Number of FCN-evaluations
C    NSTEP              Number of integration steps
C    NACCPT             Number of steps accepted
C    NREJCT             Number of steps rejected
C    NDEC               Number of decompositions
C    NSOL               Number of substitutions
C
C* Type declaration
C
      INTEGER I, J, JK, JL, JM, JMACT, JOPT, JRED, JRMAX, J1, K, KFIN,
     $     KFLAG, KM, KMACT, KOPT, K1, L, LOUT, M, MAXODE, MDT, M1, N,
     $     NACCPT, NDEC, NFCN, NJ, NREJCT, NSOL, NSTEP, NSTMAX
C
      DOUBLE PRECISION ALPHA, AWK, D, DMAX1, DMIN1, DT, DUMMY, DY, DYM,
     $     DZ, BLDFER, EPMACH, EPSAFE, ERR, FC, FCK, FCM, FCO, FMIN,
     $     FNJ, FN1, FN, H, HALF, HJ,  HJ2, HMAX, HMAXU, HN, HREST,
     $     HRTRN, OMJ, OMJO, ONE, PCT101, PCT90, QUART, RED, RMAX, RO,
     $     SAFE, SMALL, T, TEND, TN, TOL, TOLH, TOLMIN, U, Y, YK, YM,
     $     YMAX, YWGT, ZERO
C
      LOGICAL QFIRST, QKONV, QLAST, QPRMON, QPRSOL, QRED
C
      CHARACTER CHGDAT*20, PRODCT*8
C
      EXTERNAL FCN
C
C
C* Constants problem oriented: (to be supplied by the user)
C
C    MAXODE    I   K    Maximal number of first-order ODEs
C
      PARAMETER ( MAXODE = 1024          )
C
C* Constants machine oriented: (to be verified by the user)
C
C    EPMACH    D   K    Relative machine precision
C    LOUT      I   K    Output is written on logical unit LOUT
C    SMALL     D   K    Square-root of smallest positive machine number
C
CCRY  (adapted to Cray X-MP)
CCRY  PARAMETER ( EPMACH = 7.106D-15     ,
CCRY $            EPSAFE = EPMACH*10.0D0 ,
CCRY $            LOUT   = 6             ,
CCRY $            SMALL  = 6.771D-1234   )
C
CIBM  (adapted to Siemens 7.865, IBM 370-compatible)
CIBM  PARAMETER ( EPMACH = 2.22D-16      ,
CIBM $            EPSAFE = EPMACH*10.0D0 ,
CIBM $            LOUT   = 6             ,
CIBM $            SMALL  = 7.35D-40      )
C
CSUN  (adapted to sun)
C      PARAMETER ( EPMACH = 0.1085D-18    ,
C     $            EPSAFE = EPMACH*10.0D0 ,
C     $            LOUT   = 6             ,
C     $            SMALL  = 0.2223D-161   )
C
C* Other Constants:
C
C    HALF      D   K    1/2
C    ONE       D   K    1
C    PCT101    D   K    101 Percent
C    PCT90     D   K    90 Percent
C    QUART     D   K    1/4
C    ZERO      D   K    0
C
      PARAMETER ( HALF   = 0.5  D0       ,
     $            ONE    = 1.0  D0       ,
     $            PCT101 = 1.01 D0       ,
     $            PCT90  = 0.9  D0       ,
     $            QUART  = 0.25 D0       ,
     $            ZERO   = 0.0  D0       )
C
C* Control parameters: (to be supplied by the user)
C  standard values fixed below
C
C    NSTMAX    I   K    Maximum permitted number of integration steps
C                       per interval  =  10000
C    JRMAX     I   K    Maximum permitted number of stepsize reductions
C    KM        I   K    Prescribed maximum column number
C    JM        I   K    Associated maximum row number
C                       (JM = KM + 1)
C    MDT       I   K    Associated dimension of DT
C    BLDFSQ         EXT  Subroutine BLDFSQ(JM,NJ)
C                       Generate stepsize sequence with respect to /1/
C                       JM     Maximum row number
C                       NJ     Array(JM) of stepsize sequence
C    BLDFSC        EXT  Subroutine BLDFSC (MODE, Y, N, YOLD, YWGT,
C                                          YMAX, THREL, THABS)
C                       Scaling for DIFEX1
C                       MODE   ='INITIAL '    Initial scaling
C                              ='INTERNAL'    Scaling during 
C                                             discretization
C                              ='ACCEPTED'    Rescaling if step accepted
C                              Else           Error
C                       Y      Array of values Y(1),...,Y(N)
C                       N      Length of vectors Y, YOLD, YWGT, and YMAX
C                       YOLD   Array of old values
C                       YWGT   Array of scaled values
C                       YMAX   Array of maximum values
C                       THREL  Relative threshold value
C                       THABS  Absolute threshold value
C    BLDFER        EXT  Double Precision function BLDFER(Y, N, YWGT)
C                       Scaled root mean square error
C                       Y      Array of values Y(1),...,Y(N)
C                       N      Length of vectors Y and YWGT
C                       YWGT   Array of scaled values
C
      PARAMETER ( NSTMAX = 10000         ,
     $            JRMAX  = 10            ,
     $            KM     = 8             ,
     $            JM     = KM + 1        ,
     $            MDT    = MAXODE*JM     )
C
C* Internal parameters: (modification not recommended)
C
C
      PARAMETER ( FMIN   = 1.0  D-3      ,
     $            RMAX   = 0.75 D0       ,
     $            RO     = QUART         ,
     $            SAFE   = 0.7  D0       )
C
C
C* Local variables: (workspace)
C
C
C
C    QFIRST    L   V    First integration step
C    QKONV     L   V    Convergence detected
C    QLAST     L   V    Last integration step
C    QPRMON    L   V    Print integration monitor
C    QPRSOL    L   V    Print intermediate solution points
C
C* Dimensions:
C
      DIMENSION ALPHA(JM,JM), AWK(JM), D(JM,JM), DT(MAXODE,JM),
     $     DY(MAXODE), DYM(MAXODE), DZ(MAXODE), FCK(KM), NJ(JM),
     $     Y(MAXODE), YK(MAXODE), YM(MAXODE), YMAX(MAXODE),
     $     YWGT(MAXODE)
C
      COMMON /DXSTAT/ NFCN, NSTEP, NACCPT, NREJCT, NDEC, NSOL
C
C*******  Revision 1 *******  Latest change:
      DATA      CHGDAT      /'August 27, 1991'/
      DATA      PRODCT      /'BLDFX1'/
C***************************
C
C
C* Modification history
C  --------------------
C
C
C  1.0       Feb  9, 1988    First release at ZIB
C  1.1       Mar 27, 1991    Vectorize extrapolation loop,
C                            Time monitor
C  1.2       Aug 27, 1991    Allow reverse integration direction
C
C
      DATA  DT/MDT*0.D0/
C
C---1. Initial preparations
      CALL ZIBCONST(EPMACH,SMALL)
      EPSAFE = EPMACH*10.0D0
      LOUT   = 6
      QPRMON = (KFLAG .EQ. 1 .OR. KFLAG .EQ. 3)
      QPRSOL = (KFLAG .GE. 2)
      HRTRN = H
C
      DO 1001 I = 1, N
         YMAX(I) = ZERO
 1001 CONTINUE
C
      HREST = TEND - T
      H = SIGN (DMIN1 (DABS(H), DABS(HREST)), HREST)
      HMAX = DABS(HMAX)
      HMAXU = HMAX
      FCM = DMAX1 (DABS(H)/HMAX, FMIN)
      KMACT = KM
      JMACT = JM
      CALL BLDFSQ (JM, NJ)
      FN = DBLE (N)
      FN1 = DBLE (NJ(1))
      TOLH = RO*TOL
      TOLMIN = EPSAFE*FN
      IF (TOL .LT. TOLMIN) THEN
         WRITE (LOUT, 10002) PRODCT, TOL, TOLMIN
         TOL = TOLMIN
      ENDIF
C
C---  Compute amount of work per row of extrapolation tableau
      AWK(1) = FN1 + ONE
      DO 101 J=2,JM
         J1 = J - 1
         FNJ = DBLE (NJ(J))
         AWK(J) = AWK(J1) + FNJ
         DO 1011 K=1,J1
            D(J,K) = (FNJ / DBLE (NJ(K)))*(FNJ / DBLE (NJ(K)))
 1011    CONTINUE
C
         IF (J .NE. 2) THEN
            DO 1012 K1=2,J1
               K = K1 - 1
               ALPHA(J1,K) = TOLH**((AWK(K1) - AWK(J)) /
     $              ((AWK(J) - AWK(1) + ONE)*DBLE(K + K1)))
 1012       CONTINUE
C
         ENDIF
 101  CONTINUE
C
C---1.2 Determination of maximum column number in extrapolation
C---    tableau (information theoretic concept, ref./3/)
      KOPT = 1
      JOPT = 2
 121  CONTINUE
C     DO WHILE (JOPT .LT. KMACT .AND.
C               AWK(JOPT+1)*PCT101 .LE. AWK(JOPT)*ALPHA(JOPT,KOPT))
         IF (JOPT .GE. KMACT .OR.
     $     AWK(JOPT+1)*PCT101 .GT. AWK(JOPT)*ALPHA(JOPT,KOPT)) GOTO 122
C                                                              Exit 121
         KOPT = JOPT
         JOPT = JOPT + 1
         GOTO  121
C     ENDDO
 122  KMACT = KOPT + 1
      JMACT = JOPT
      IF (QPRMON) WRITE (LOUT, 11221)
     $     PRODCT, CHGDAT, TOL, KMACT, NJ
C
      IF (QPRSOL) WRITE (LOUT, 11222)
      NSTEP = 0
      QFIRST = .TRUE.
      QLAST = .FALSE.
C      NFCN = 0
      KFIN = 0
      OMJO = ZERO
      CALL BLDFSC ('INITIAL ', Y, N, DUMMY, YWGT, YMAX, TOL, ONE)
C
C---2. Basic integration step
 2    CONTINUE
C     DO WHILE (T .NE. TEND)
         IF (QPRMON) WRITE (LOUT, 12001) NSTEP,NFCN,T,H,KFIN,KOPT
         IF (QPRSOL) WRITE (LOUT, 12002) NSTEP,NFCN,T,H,(Y(I),I=1,N)
         JRED = 0
C
C---     Explicit euler starting step
         CALL FCN (N, T, Y, DZ)
         NFCN = NFCN + 1
C
C---3.   Basic discretization step
 3       CONTINUE
C        DO WHILE (JRED .LE. JRMAX .AND. .NOT. QKONV)
            IF (QLAST) THEN
               TN = TEND
            ELSE
               TN = T + H
            ENDIF
            IF (TN .EQ. T) THEN
C              Error 4
               IF (QPRMON) WRITE (LOUT, 13001) PRODCT
               KFLAG = -4
               GOTO  9
C              Exit to Return
            ENDIF
C
C---3.1     Internal discretization
            DO 31 J=1,JMACT
               M = NJ(J)
               M1 = M - 1
               KFIN = J - 1
               FNJ = DBLE (M)
               HJ = H / FNJ
               HJ2 = HJ + HJ
               DO 3101 I=1,N
                  YK(I) = Y(I)
                  YM(I) = Y(I) + HJ*DZ(I)
 3101          CONTINUE
C
C---3.1.3      Explicit mid-point rule
               DO 313 K=1,M1
                  CALL FCN (N, T + HJ*DBLE (K), YM, DY)
                  NFCN = NFCN + 1
                  DO 3135 I=1,N
                     U = YK(I) + HJ2*DY(I)
                     YK(I) = YM(I)
                     YM(I) = U
 3135             CONTINUE
 313           CONTINUE
C
C---3.1.4      Smoothing final step
               CALL FCN (N, TN, YM, DY)
               NFCN = NFCN + 1
               DO 3141 I = 1,N
                  YM(I) = (YM(I) + YK(I) + HJ*DY(I))*HALF
 3141          CONTINUE
C
C
C---3.1.5      Extrapolation
               DO 3153 I=1,N
                  DY(I) = YM(I)
                  YK(I) = DT(I,1)
                  DT(I,1) = DY(I)
 3153          CONTINUE
C
               DO 3158 K=2,J
                  JK = J - K + 1
C
                  DO 3155 I=1,N
                     DYM(I) = (DY(I) - YK(I)) / (D(J,JK) - ONE)
                     DY(I) = D(J,JK)*DYM(I)
                     YK(I) = DT(I,K)
                     DT(I,K) = DYM(I)
                     YM(I) = DYM(I) + YM(I)
 3155             CONTINUE
C
 3158          CONTINUE
C
               IF (J .NE. 1) THEN
C
C---3.1.6         Convergence monitor
                  CALL BLDFSC ('INTERNAL',YM,N,Y,YWGT,YMAX,TOL,ONE)
                  ERR = BLDFER (DYM, N, YWGT)
                  QKONV = ERR .LE. TOL
                  ERR = ERR / TOLH
C
C---              Order control
                  K = J - 1
                  FC = ERR**(ONE / DBLE(K + J))
                  FCK(K) = FC
C
C---              Order window
                  IF (J .GE. KOPT .OR. QFIRST .OR. QLAST) THEN
                     IF (QKONV) GOTO 25
C                                Exit 3 for next basic integration step
C
C---                 Check for possible stepsize reduction
                     RED = ONE / FC
                     QRED = .FALSE.
                     IF (K .EQ. KMACT .OR. K .EQ. JOPT) THEN
                        RED = RED*SAFE
                        QRED = .TRUE.
                     ELSE
                        IF (K .EQ. KOPT) THEN
                           RED = RED*ALPHA(JOPT,KOPT)
                           IF (RED .LT. ONE) THEN
                              RED = ONE / FC
                              QRED = .TRUE.
                           ENDIF
                        ELSE
                           IF (KOPT .EQ. KMACT) THEN
                              RED = RED*ALPHA(KMACT,K)
                              IF (RED .LT. ONE) THEN
                                 RED = RED * SAFE
                                 QRED = .TRUE.
                              ENDIF
                           ELSE
                              RED = RED*ALPHA(JOPT,K)
                              IF (RED .LT. ONE) THEN
                                 RED = ALPHA(KOPT,K) / FC
                                 QRED = .TRUE.
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                     IF (QRED) GOTO 32
C                              Exit 3.1 to stepsize reduction
                  ENDIF
               ENDIF
 31         CONTINUE
C
C---3.2     Prepare stepsize reduction
 32         CONTINUE
C
C---3.5     Stepsize reduction
            RED = DMIN1 (RED, RMAX)
            H = H*RED
            IF (NSTEP .GT. 0) QLAST = .FALSE.
            JRED = JRED + 1
            IF (QPRMON) WRITE (LOUT, 13501) JRED,RED,
     $           KFIN,KOPT,KMACT
            IF (JRED .GT. JRMAX) THEN
C              Error 3
               IF (QPRMON) WRITE (LOUT, 13502) JRMAX
               KFLAG = -3
               GOTO  9
C              Exit to Return
            ENDIF
            GOTO  3
C        ENDDO
C
C        ************************************************
C---2.5  Preparations for next basic integration step
 25      NSTEP = NSTEP + 1
         QFIRST = .FALSE.
         IF (NSTEP .GT. NSTMAX) THEN
C           Error 2
C           Emergency exit, if too many steps taken
            IF (QPRMON) WRITE (LOUT, 12501) PRODCT, NSTMAX
            KFLAG = -2
            GOTO  9
C           Exit to return
         ENDIF
C
C---     Restoring
         DO 251 I=1, N
            Y(I) = YM(I)
 251     CONTINUE
C
         T = TN
         IF (T .EQ. TEND) GOTO 9
C                         Exit to return
         CALL BLDFSC ('ACCEPTED', Y, N, DUMMY, YWGT, YMAX, TOL, ONE)
C
C---2.7  Order and stepsize selection
C
C---2.7.1 Stepsize restrictions
         HMAX = DMIN1(HMAXU, DABS(H)/FMIN)
         FCM = DABS(H) / HMAX
C
C---2.7.2 Optimal order determination
         KOPT = 1
         JOPT = 2
         FCO = DMAX1 (FCK(1), FCM)
         OMJO = FCO*AWK(2)
         IF (KFIN .GE. 2) THEN
            DO 272 L=2,KFIN
               JL = L + 1
               FC = DMAX1 (FCK(L), FCM)
               OMJ = FC*AWK(JL)
               IF (OMJ*PCT101 .LE. OMJO .AND. L .LT. KMACT) THEN
                  KOPT = L
                  JOPT = JL
                  OMJO = OMJ
                  FCO = FC
               ENDIF
 272        CONTINUE
         ENDIF
         HREST = TEND - T
         HN = H / FCO
C
C---2.7.3 Possible increase of order
         IF (DABS(HN) .LT. DABS(HREST)) THEN
            IF ((JRED .EQ. 0 .OR. NSTEP .EQ. 0) .AND.
     $           KOPT .GE. KFIN .AND. KOPT .NE. KMACT) THEN
               FC = DMAX1 (FCO/ALPHA(JOPT,KOPT), FCM)
               JL = JOPT + 1
               IF (AWK(JL)*FC*PCT101 .LE. OMJO .AND.
     $              JOPT .LT. KMACT) THEN
                  FCO = FC
                  HN = H / FCO
                  KOPT = JOPT
                  JOPT = JOPT + 1
               ENDIF
            ENDIF
         ENDIF
C
C---2.7.4 Stepsize selection
         H = HN
         HRTRN = H
         IF (DABS(H) .GT. DABS(HREST)*PCT90) THEN
            H = HREST
            QLAST = .TRUE.
         ENDIF
         GOTO  2
C     ENDDO
C
C---9. Exit
 9    HMAX = HMAXU
      IF (KFLAG .LT. 0) THEN
C        Fail exit
         H = ZERO
      ELSE
C        Solution exit
         H = HRTRN
         IF (QPRMON) WRITE (LOUT, 12001) NSTEP,NFCN,T,H,KFIN,KOPT
         IF (QPRSOL) WRITE (LOUT, 12002) NSTEP,NFCN,T,H,(Y(I),I=1,N)
      ENDIF
      RETURN
C
C
10001 FORMAT(//,' ',A8,'  - Error -  '
     $      ,   ' Direction if integration is reverse to convention.')
10002 FORMAT(//,' ',A8,'  - Warning -'
     $      ,   ' Desired tolerance ', D10.3, ' too small.', /,
     $      22X,' tolerance set to  ', D10.3, '.')
C
11221 FORMAT(1H0,A8,' - ',A20,/,
     $       1H0,' Rel.prec. TOL ',D10.3,' max.col. ',I3, /,
     $       ' sequence ',(1H ,13I4))
11222 FORMAT(//,5X,4HStep,3X,7HF-calls,8X,1HT,25X,1HH,5X,7HY1(T)..,//)
12001 FORMAT(1H ,2I9,D20.11,D12.4,I9,I6)
12002 FORMAT(1H ,2I9,D20.11,D12.4,4D20.11,/,(1H ,50X,4D20.11))
12501 FORMAT(//,' ',A8,'  - Error -  '
     $      ,18H more than NSTMAX=,I3,18H integration steps,//)
13001 FORMAT(//,' ',A8,'  - Error -  '
     $      ,40H stepsize reduction failed to succeed  ,//)
13501 FORMAT(1H ,I3,27H Stepsize reduction factor ,D10.3,
     $      ' KFIN',I3,' KOPT',I3,' KMAX',I3)
13502 FORMAT(//,' ',A8,'  - Error -  '
     $      ,17H more than JRMAX=,I3,29H stepsize reductions per step,/)
C
C
C End BLDFX1
C
      END
      SUBROUTINE BLDFSQ(M,NJ)
      INTEGER I, M, NJ
      DIMENSION NJ(M)
C
C  Set stepsize sequence for DIFEX1
C
      NJ(1) = 2
      DO 10 I=2,M
         NJ(I) = NJ(I-1) + 2
 10   CONTINUE
C
      RETURN
      END
      SUBROUTINE BLDFSC (MODE, Y, N, YOLD, YWGT, YMAX, THREL, THABS)
C
C     Scaling for DIFEX1
C
C       (May be altered for real life applications
C        by the skillful user)
C
C
C* Parameters:
C
C    MODE      C*8 IN   ='INITIAL '    Initial scaling
C                       ='INTERNAL'    Scaling during discretization
C                       ='ACCEPTED'    Rescaling if step accepted
C                       Else           Error
C    Y         D   IN   Array of values Y(1),...,Y(N)
C    N         I   IN   Length of vectors Y, YOLD, YWGT, and YMAX
C    YOLD      D   IN   Array of old values
C    YWGT      D   OUT  Array of scaled values new
C    YMAX      D   IN   Array of maximum values old
C                  OUT  Array of maximum values new
C    THREL     D   IN   Relative threshold value
C    THABS     D   IN   Absolute threshold value
C
C* Local variables:
C
C    YUSER     D   V    User defined array of maximum values
C
C* Type declaration
C
      INTEGER I, LOUT, MAXODE, N
C
      DOUBLE PRECISION DABS, DMAX1, EPMACH, ONE, THABS, THREL, U, Y,
     $     YMAX, YOLD, YUSER, YWGT, ZERO, SMALL
C
      CHARACTER MODE*8
C
C* Constants:
C
C    EPMACH    D   K    Relative machine precision
C    LOUT      I   K    Output is written on logical unit LOUT
C    MAXODE    I   K    Maximal number of first-order ODEs
C    ONE       D   K    1.0
C    ZERO      D   K    0.0
C
CCRY  (adapted to Cray X-MP)
CCRY  PARAMETER ( EPMACH = 7.106D-15     ,
C
CIBM  (adapted to Siemens 7.865, IBM 370-compatible)
CIBM  PARAMETER ( EPMACH = 2.22D-16      ,
C
CSUN  (adapted to sun)
      PARAMETER (
     $            LOUT   = 6             ,
     $            MAXODE = 1024          ,
     $            ONE    = 1.0  D0       ,
     $            ZERO   = 0.0  D0       )
C
      DIMENSION Y(N), YOLD(N), YWGT(N), YMAX(N), YUSER(MAXODE)
      SAVE YUSER
      CALL ZIBCONST(EPMACH,SMALL)
      IF (MODE .EQ.          'INITIAL '         ) THEN
C                             --------
         DO 100 I=1,N
            YUSER(I) = DABS (YMAX(I))
            U = DABS (Y(I))
            IF (U .LT. EPMACH) U = ONE
            YMAX(I) = DMAX1 (U, YUSER(I), THABS)
            YWGT(I) = YMAX(I)
 100     CONTINUE
C
      ELSE IF (MODE .EQ.     'INTERNAL'         ) THEN
C                             --------
         DO 200 I=1,N
            YWGT(I) = DMAX1 (YMAX(I)*THREL, DABS(Y(I)),
     $                       DABS(YOLD(I)), YUSER(I), THABS)
 200     CONTINUE
C
      ELSE IF (MODE .EQ.     'ACCEPTED'         ) THEN
C                             --------
         DO 300 I=1,N
            YMAX(I) = DMAX1 (YMAX(I), DABS(Y(I)))
 300     CONTINUE
C
      ELSE
         WRITE (LOUT, '(//,A,/)')
     $      ' D1SCAL    - ERROR -   Illegal mode'
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION BLDFER(Y, N, YWGT)
C* Title:
C
C  Scaled root mean square error
C
C
C* Parameters:
C
C    Y         D   IN   Array of values Y(1),...,Y(N)
C    N         I   IN   Length of Vectors Y and YWGT
C    YWGT      D   IN   Array of scaled values
C
C* Type declaration
C
      INTEGER I, N
C
      DOUBLE PRECISION DBLE, DSQRT, SUM, Y, YWGT, ZERO
C
C* Constants:
C
C    ZERO      D   K    0
C
      PARAMETER ( ZERO   = 0.0  D0       )
C
      DIMENSION Y(N), YWGT(N)
C
      SUM = ZERO
      DO 100 I=1,N
         SUM = SUM + (Y(I) / YWGT(I)) * (Y(I) / YWGT(I))
 100  CONTINUE
C
      BLDFER = DSQRT(SUM / DBLE(N))
      RETURN
      END
