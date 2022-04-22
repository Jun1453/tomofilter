      PROGRAM RCONST
      IMPLICIT NONE
*
*     This routine reads singular matrices and reconstructs
*     matrix G for testing
*     Modified from the Example Program PDGESV of ScaLAPACK
*
*     .. Parameters ..
      INTEGER            DLEN_, IA, JA, IU, JU, M, N, MB, NB, RSRC,
     $                   CSRC, MXLLDA, MXLLDB, NRHS, NBRHS, NOUT,
     $                   MXLOCR, MXLOCC, MXRHSC, GB,
     $                   IVT, JVT, ISIG, JSIG
      PARAMETER          ( DLEN_ = 9, IA = 1, JA = 1, IU = 1, JU = 1,
     $                   IVT = 1, JVT = 1, ISIG = 1, JSIG = 1,
     $                   M = 251617, N = 36092, MB = 200, NB = 200,
     $                   RSRC = 0, CSRC = 0, MXLLDA = 251617,
     $                   NRHS = 1, NBRHS = 1, NOUT = 6,
     $                   MXLOCR = 9000, MXLOCC = 3092 )
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      CHARACTER          JOBU, JOBVT
      PARAMETER          ( JOBU = 'V', JOBVT = 'V' )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICTXT, INFO, MYCOL, MYROW, NPCOL, NPROW,
     $                   LDA, LDU, LDVT, NQ, WPDGESVD,
     $                   PTRA, PTRAC, PTRD, PTRWORK, PTRS, PTRSC, PTRU,
     $                   PTRR, PTRUC, PTRVT, PTRVTC, SETHET, SIZE, SIZEQ
      DOUBLE PRECISION   ANORM, BNORM, EPS, RESID, XNORM, SIGMAX, SIGK
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCU( DLEN_ ),
     $                   DESCVT( DLEN_ ), DESCSIG( DLEN_ ),
     $                   DESCR( DLEN_ ), DESCS( DLEN_ )
!      $                   IPIV( MXLOCR+NB )
!       DOUBLE PRECISION   A( MXLLDA, MXLOCC ), A0( MXLLDA, MXLOCC ),
!      $                   B( MXLLDB, MXRHSC ), B0( MXLLDB, MXRHSC ),
      DOUBLE PRECISION   WORK( 400000000 )
!      $                   A( MXLLDA, MXLOCC ), U( MXLLDA, MXLOCC ),
!      $                   D( MXLOCC ),    E( MXLOCR ),
!      $                   TAUQ( MXLOCC ), TAUP( MXLOCR ),
!      $                   VT( MXLLDA, MXLOCC ), S( MXLLDA ),          
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      DOUBLE PRECISION   PDLAMCH, PDLANGE
      EXTERNAL           PDLAMCH, PDLANGE, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GRIDEXIT, BLACS_GRIDINFO,
     $                   DESCINIT, MATINIT, PDGEMM, PDGESV, PDLACPY,
     $                   SL_INIT, PDGESVD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN
*     ..
*     .. Data statements ..
      DATA               NPROW / 28 / , NPCOL / 12 /
*     ..
*     .. My Variables ..
      CHARACTER*40 FILENAM
*     ..
*     .. Executable Statements ..
*
*     INITIALIZE THE PROCESS GRID
*
      CALL SL_INIT( ICTXT, NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     If I'm not in the process grid, go to the end of the program
*
      IF( MYROW.EQ.-1 )
     $   GO TO 10
*
*     DISTRIBUTE THE MATRIX ON THE PROCESS GRID
*     Initialize the array descriptors for the matrices A and B
*
      SIZE = MIN( M, N )
      LDA = NUMROC( M, NB, MYROW, 0, NPROW )
      LDA = MAX( 1, LDA )
      NQ = NUMROC( N, NB, MYCOL, 0, NPCOL )
      LDU = LDA
      SIZEQ = NUMROC( SIZE, NB, MYCOL, 0, NPCOL )
      LDVT = NUMROC( SIZE, NB, MYROW, 0, NPROW )
      LDVT = MAX( 1, LDVT )
      CALL DESCINIT( DESCA, M, N, NB, NB, 0, 0, ICTXT, LDA, INFO )
      CALL DESCINIT( DESCU, M, SIZE, NB, NB, 0, 0, ICTXT, LDU, INFO )
      CALL DESCINIT( DESCVT, SIZE, N, NB, NB, 0, 0, ICTXT, LDVT, INFO )
      CALL DESCINIT( DESCS, SIZE, SIZE, NB, NB, 0, 0, ICTXT, LDVT, INFO)
      CALL DESCINIT( DESCR, SIZE, SIZE, NB, NB, 0, 0, ICTXT, LDVT, INFO)

*
*     Set some pointers to work array in order to do "dummy" calls.
*
      PTRA = 2
      PTRU = PTRA + LDA*NQ
      PTRR = PTRU + LDU*SIZEQ
      PTRS = PTRR + LDVT*NQ
      PTRVT = PTRS + LDVT*NQ
      
      PTRWORK = PTRVT + LDVT*NQ
*
*     "Dummy" calls -- return required workspace in work(1) without
*     any calculation. 
*
      !   CALL PDGESVD( 'V', 'V', M, N, WORK( PTRA ), IA, JA, DESCA, 
      !  $              WORK( PTRS ), WORK( PTRU ), IU, JU, DESCU, 
      !  $              WORK( PTRVT ), IVT, JVT, DESCVT,
      !  $              WORK( PTRWORK ), -1, INFO )     !   WPDGESVD = INT( WORK( PTRWORK ) )
      !   WPDGESVD = INT( WORK( PTRWORK ) )
      ! WRITE(*,*) WPDGESVD
*
*     Generate matrices A and B and distribute to the process grid
*
      ! CALL MATINIT( A, DESCA, B, DESCB )
      CALL PDLAREAD( 'U.out', WORK( PTRU ), DESCU, 0, 0,
     $                WORK(PTRWORK) )
      CALL PDLAREAD( 'V.out', WORK( PTRVT ), DESCVT, 0, 0,
     $                WORK(PTRWORK) )
      CALL PDLAREAD( 'S.out', WORK( PTRS ), DESCS, 0, 0,
     $                WORK(PTRWORK) )
      ! CALL PDLAREAD( 'obs_example', B, DESCB, 0, 0, WORK )
     
*     Invert diagonal matrix S to Z
      !CALL PDELSET( WORK( PTRA ), 1, 1, DESCA, DBLE(114514.1919) )
      CALL PDELGET( 'N', 'I', SIGMAX, WORK( PTRA ), 1, 1, DESCA )
      !WRITE(*,*) MYROW,MYCOL,SIGMAX, 0.001*SIGMAX
*     Update singular values and set zeros for elements under lowcut
      DO 30 I = 1, SIZE
         CALL PDELGET( 'N', 'I', SIGK, WORK( PTRS ), I, I, DESCS )
         IF ( SIGK .LT. 0.001*SIGMAX .OR. I .GT. 6000 ) THEN
            !WRITE(*,*) MYROW,MYCOL,I,SIGK, SIGMAX
            CALL PDELSET( WORK( PTRWORK ), I, I, DESCS, ZERO )
         ELSE
            CALL PDELSET( WORK( PTRS ), I, I, DESCS, DBLE(SIGK) )
         ENDIF
   30 CONTINUE
*
*     Make a copy of A and B for checking purposes
*
      ! CALL PDLACPY( 'All', N, N, A, 1, 1, DESCA, A0, 1, 1, DESCA )
      ! CALL PDLACPY( 'All', N, NRHS, B, 1, 1, DESCB, B0, 1, 1, DESCB )
*
*     CALL THE SCALAPACK ROUTINE
*     Solve the linear system A * X = B
*
*     R = S * VT
      CALL PDGEMM( 'N', 'N', N, N, N, ONE,
     $             WORK( PTRS ), 1, 1, DESCS,
     $             WORK( PTRVT ), 1, 1, DESCVT,
     $             ZERO, WORK( PTRR ), 1, 1, DESCR )
*     A = U * R
      CALL PDGEMM( 'N', 'N', M, N, N, ONE,
     $             WORK( PTRU ), 1, 1, DESCU,
     $             WORK( PTRR ), 1, 1, DESCR,
     $             ZERO, WORK( PTRA ), 1, 1, DESCA )
*     V = R * Z
*      CALL PDGEMM( 'N', 'N', N, N, N, ONE,
*     $             WORK( PTRR ), 1, 1, DESCR,
*     $             WORK( PTRS ), 1, 1, DESCS,
*     $             ZERO, WORK( PTRVT ), 1, 1, DESCVT )
    !   CALL PDGEMM( 'N', 'N', N, NRHS, N, ONE, A0, 1, 1, DESCA, B,
    !  $             1, 1, DESCB, -ONE, B0, 1, 1, DESCB )
      ! WRITE( *, * )IA, JA
      ! PTRVT = PTRU + LDU*SIZEQ
      ! PTRWORK = PTRVT + LDVT*NQ
   !    CALL PDGESVD( JOBU, JOBVT, M, N,
   !   $              WORK( PTRA ), IA, JA, DESCA, 
   !   $              WORK( PTRS ),
   !   $              WORK( PTRU ), IU, JU, DESCU,
   !   $              WORK( PTRVT ), IVT, JVT, DESCVT,
   !   $              WORK( PTRWORK ), WPDGESVD, INFO )
!       CALL PDGESVD('V','V',M,N,A,IA,JA,DESCA,S,U,IU,JU,DESCU,
!      $             VT,IVT,JVT,DESCVT,WORK,-1,INFO)
!       ( N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB,
!      $             INFO )
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9999 )
         WRITE( NOUT, FMT = 9998 )M, N, NB
         WRITE( NOUT, FMT = 9997 )NPROW*NPCOL, NPROW, NPCOL
         WRITE( NOUT, FMT = 9996 )INFO
      END IF
! *
! *     Compute residual ||A * X  - B|| / ( ||X|| * ||A|| * eps * N )
! *
!       EPS = PDLAMCH( ICTXT, 'Epsilon' )
!       ANORM = PDLANGE( 'I', N, N, A, 1, 1, DESCA, WORK )
!       BNORM = PDLANGE( 'I', N, NRHS, B, 1, 1, DESCB, WORK )
!       CALL PDGEMM( 'N', 'N', N, NRHS, N, ONE, A0, 1, 1, DESCA, B, 1, 1,
!      $             DESCB, -ONE, B0, 1, 1, DESCB )
!       XNORM = PDLANGE( 'I', N, NRHS, B0, 1, 1, DESCB, WORK )
!       RESID = XNORM / ( ANORM*BNORM*EPS*DBLE( N ) )
! *
!       IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
!          IF( RESID.LT.10.0D+0 ) THEN
!             WRITE( NOUT, FMT = 9995 )
!             WRITE( NOUT, FMT = 9993 )RESID
!          ELSE
!             WRITE( NOUT, FMT = 9994 )
!             WRITE( NOUT, FMT = 9993 )RESID
!          END IF
!       END IF
*
*     My code
*
      ! WRITE( FILENAM, FMT = 8999 )'d', MYROW, MYCOL
      CALL PDLAWRITE( 'A.out', M, N, WORK( PTRA ), IA, JA, DESCA,
     $                0, 0, WORK(PTRWORK) )
   !    CALL PDLAWRITE( 'U.out', M, SIZE, WORK( PTRU ), IU, JU, DESCU,
   !   $                0, 0, WORK(PTRWORK) )
   !    CALL PDLAWRITE( 'test.out', M, N, WORK( PTRA ), IA, JA,
   !   $                DESCA, 0, 0, WORK(PTRWORK) )
   !    CALL PDLAWRITE( 'R.out', SIZE, SIZE, WORK( PTRVT ), IVT, JVT,
   !   $                DESCVT, 0, 0, WORK(PTRWORK) )
!       CALL PDLAWRITE( 'test.out.obs', M, NRHS, B, IB, JB, DESCB, 0, 0,
!      $                WORK )
      ! WRITE( *, * ) D(2), E(2)
*
*              Diagonal Matrix with specified singular values.
*
!       CALL PDLASET( 'All', SIZE, SIZE, ZERO, ZERO, WORK( PTRWORK ), 
!      $              ISIG, JSIG, DESCSIG )
! *
!       DO 40 I = 1, SIZE
!          CALL PDELSET( WORK( PTRWORK ), I, I, DESCSIG,
!      $                 WORK( PTRS+I-1 ) )
!    40 CONTINUE

!       CALL PDLAWRITE( 'S.out', SIZE, N, WORK( PTRWORK ), ISIG, JSIG,
!      $                 DESCSIG, 0, 0, WORK(PTRWORK) )
*
*     RELEASE THE PROCESS GRID
*     Free the BLACS context
*
      CALL BLACS_GRIDEXIT( ICTXT )
   10 CONTINUE
*
*     Exit the BLACS
*
      CALL BLACS_EXIT( 0 )
*
 9999 FORMAT( / 'ScaLAPACK Example Program #1 -- May 1, 1997' )
 9998 FORMAT( / 'Solving Ax=b where A is a ', I3, ' by ', I3,
     $      ' matrix with a block size of ', I3 )
 9997 FORMAT( 'Running on ', I3, ' processes, where the process grid',
     $      ' is ', I3, ' by ', I3 )
 9996 FORMAT( / 'INFO code returned by PDGESV = ', I3 )
 9995 FORMAT( /
     $   'According to the normalized residual the solution is correct.'
     $       )
 9994 FORMAT( /
     $ 'According to the normalized residual the solution is incorrect.'
     $       )
 9993 FORMAT( / '||A*x - b|| / ( ||x||*||A||*eps*N ) = ', 1P, E16.8 )
 8999 FORMAT( A, 2I3.3)
      STOP
      END
      SUBROUTINE MATINIT( AA, DESCA, B, DESCB )
*
*     MATINIT generates and distributes matrices A and B (depicted in
*     Figures 2.5 and 2.6) to a 2 x 3 process grid
*
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   AA( * ), B( * )
*     ..
*     .. Parameters ..
      INTEGER            CTXT_, LLD_
      PARAMETER          ( CTXT_ = 2, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICTXT, MXLLDA, MYCOL, MYROW, NPCOL, NPROW
      DOUBLE PRECISION   A, C, K, L, P, S
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Executable Statements ..
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      S = 19.0D0
      C = 3.0D0
      A = 1.0D0
      L = 12.0D0
      P = 16.0D0
      K = 11.0D0
*
      MXLLDA = DESCA( LLD_ )
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         AA( 1 ) = S
         AA( 2 ) = -S
         AA( 3 ) = -S
         AA( 4 ) = -S
         AA( 5 ) = -S
         AA( 1+MXLLDA ) = C
         AA( 2+MXLLDA ) = C
         AA( 3+MXLLDA ) = -C
         AA( 4+MXLLDA ) = -C
         AA( 5+MXLLDA ) = -C
         AA( 1+2*MXLLDA ) = A
         AA( 2+2*MXLLDA ) = A
         AA( 3+2*MXLLDA ) = A
         AA( 4+2*MXLLDA ) = A
         AA( 5+2*MXLLDA ) = -A
         AA( 1+3*MXLLDA ) = C
         AA( 2+3*MXLLDA ) = C
         AA( 3+3*MXLLDA ) = C
         AA( 4+3*MXLLDA ) = C
         AA( 5+3*MXLLDA ) = -C
         B( 1 ) = 0.0D0
         B( 2 ) = 0.0D0
         B( 3 ) = 0.0D0
         B( 4 ) = 0.0D0
         B( 5 ) = 0.0D0
      ELSE IF( MYROW.EQ.0 .AND. MYCOL.EQ.1 ) THEN
         AA( 1 ) = A
         AA( 2 ) = A
         AA( 3 ) = -A
         AA( 4 ) = -A
         AA( 5 ) = -A
         AA( 1+MXLLDA ) = L
         AA( 2+MXLLDA ) = L
         AA( 3+MXLLDA ) = -L
         AA( 4+MXLLDA ) = -L
         AA( 5+MXLLDA ) = -L
         AA( 1+2*MXLLDA ) = K
         AA( 2+2*MXLLDA ) = K
         AA( 3+2*MXLLDA ) = K
         AA( 4+2*MXLLDA ) = K
         AA( 5+2*MXLLDA ) = K
      ELSE IF( MYROW.EQ.0 .AND. MYCOL.EQ.2 ) THEN
         AA( 1 ) = A
         AA( 2 ) = A
         AA( 3 ) = A
         AA( 4 ) = -A
         AA( 5 ) = -A
         AA( 1+MXLLDA ) = P
         AA( 2+MXLLDA ) = P
         AA( 3+MXLLDA ) = P
         AA( 4+MXLLDA ) = P
         AA( 5+MXLLDA ) = -P
      ELSE IF( MYROW.EQ.1 .AND. MYCOL.EQ.0 ) THEN
         AA( 1 ) = -S
         AA( 2 ) = -S
         AA( 3 ) = -S
         AA( 4 ) = -S
         AA( 1+MXLLDA ) = -C
         AA( 2+MXLLDA ) = -C
         AA( 3+MXLLDA ) = -C
         AA( 4+MXLLDA ) = C
         AA( 1+2*MXLLDA ) = A
         AA( 2+2*MXLLDA ) = A
         AA( 3+2*MXLLDA ) = A
         AA( 4+2*MXLLDA ) = -A
         AA( 1+3*MXLLDA ) = C
         AA( 2+3*MXLLDA ) = C
         AA( 3+3*MXLLDA ) = C
         AA( 4+3*MXLLDA ) = C
         B( 1 ) = 1.0D0
         B( 2 ) = 0.0D0
         B( 3 ) = 0.0D0
         B( 4 ) = 0.0D0
      ELSE IF( MYROW.EQ.1 .AND. MYCOL.EQ.1 ) THEN
         AA( 1 ) = A
         AA( 2 ) = -A
         AA( 3 ) = -A
         AA( 4 ) = -A
         AA( 1+MXLLDA ) = L
         AA( 2+MXLLDA ) = L
         AA( 3+MXLLDA ) = -L
         AA( 4+MXLLDA ) = -L
         AA( 1+2*MXLLDA ) = K
         AA( 2+2*MXLLDA ) = K
         AA( 3+2*MXLLDA ) = K
         AA( 4+2*MXLLDA ) = K
      ELSE IF( MYROW.EQ.1 .AND. MYCOL.EQ.2 ) THEN
         AA( 1 ) = A
         AA( 2 ) = A
         AA( 3 ) = -A
         AA( 4 ) = -A
         AA( 1+MXLLDA ) = P
         AA( 2+MXLLDA ) = P
         AA( 3+MXLLDA ) = -P
         AA( 4+MXLLDA ) = -P
      END IF
      RETURN
      END
      SUBROUTINE SL_INIT( ICTXT, NPROW, NPCOL )
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, NPCOL, NPROW
*     ..
*
*  Purpose
*  =======
*
*  SL_INIT initializes an NPROW x NPCOL process grid using a row-major
*  ordering  of  the  processes. This routine retrieves a default system
*  context  which  will  include all available processes. In addition it
*  spawns the processes if needed.
*
*  Arguments
*  =========
*
*  ICTXT   (global output) INTEGER
*          ICTXT specifies the BLACS context handle identifying the
*          created process grid.  The context itself is global.
*
*  NPROW   (global input) INTEGER
*          NPROW specifies the number of process rows in the grid
*          to be created.
*
*  NPCOL   (global input) INTEGER
*          NPCOL specifies the number of process columns in the grid
*          to be created.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            IAM, NPROCS
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GET, BLACS_GRIDINIT, BLACS_PINFO,
     $                   BLACS_SETUP
*     ..
*     .. Executable Statements ..
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
*
*     If machine needs additional set up, do it now
*
      IF( NPROCS.LT.1 ) THEN
         IF( IAM.EQ.0 )
     $      NPROCS = NPROW*NPCOL
         CALL BLACS_SETUP( IAM, NPROCS )
      END IF
*
*     Define process grid
*
      CALL BLACS_GET( -1, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
*
      RETURN
*
*     End of SL_INIT
*
      END
