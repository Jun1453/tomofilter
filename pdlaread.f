      SUBROUTINE PDLAREAD( FILENAM, A, DESCA, IRREAD, ICREAD,
     $                     WORK )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     August 12, 2001 
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    FILENAM
      INTEGER            ICREAD, IRREAD
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), WORK( * )
*     ..
*
* Purpose
* =======
*
*  PDLAREAD reads from a file named FILNAM a matrix and distribute
*  it to the process grid.
*
*  Only the process of coordinates {IRREAD, ICREAD} read the file.
*
*  WORK must be of size >= MB_ = DESCA( MB_ ).
*
* Further Details
* ===============
*
* Contributed by Song Jin, University of Tennessee, 1996.
*
* =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            NIN
      PARAMETER          ( NIN = 11 )
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ISIOPROCESSOR
      INTEGER            CSRC, I, ICTXT, IEND, ISIZE, ISTART, J, JEND,
     $                   JSIZE, JSTART, LDD, LWORK, M, MB, MM, MYCOL,
     $                   MYROW, N, NB, NN, NPCOL, NPROW, RSRC
      DOUBLE PRECISION   TMP
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Local Arrays ..
      INTEGER            DESCWORK( DLEN_ ), IWORK( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, IGEBR2D, IGEBS2D,
     $                   PDGEADD, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
*     ..
*     .. Executable Statements ..
      LWORK = DESCA( MB_ )
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      ISIOPROCESSOR = ( ( MYROW.EQ.IRREAD ) .AND. ( MYCOL.EQ.ICREAD ) )
*
      IF( ISIOPROCESSOR ) THEN
        ! TODO: read binary file
         OPEN( NIN, FILE = FILENAM, STATUS = 'OLD',
     $        FORM = 'UNFORMATTED', ACCESS = 'STREAM', ERR = 50 )
         REWIND ( NIN )
         READ( NIN, ERR = 50 ) IWORK( 1 )
         READ( NIN, ERR = 50 ) IWORK( 2 )
         CALL IGEBS2D( ICTXT, 'All', ' ', 2, 1, IWORK, 2 )
      ELSE
         CALL IGEBR2D( ICTXT, 'All', ' ', 2, 1, IWORK, 2, IRREAD,
     $                 ICREAD )
      END IF
      M = IWORK( 1 )
      N = IWORK( 2 )
*
      MM = MAX( 1, MIN( M, LWORK ) )
      NN = MAX( 1, INT( LWORK / MM ) )
      MB = MM
      NB = NN
      RSRC = IRREAD
      CSRC = ICREAD
      LDD = MAX( 1, MM )
      CALL DESCSET( DESCWORK, MM, NN, MB, NB, RSRC, CSRC, ICTXT, LDD )
*
      DO 40 JSTART = 1, N, NN
         JEND = MIN( N, JSTART+NN-1 )
         JSIZE = JEND - JSTART + 1
         DO 30 ISTART = 1, M, MM
            IEND = MIN( M, ISTART+MM-1 )
            ISIZE = IEND - ISTART + 1
            ALPHA = ONE
            BETA = ZERO
            IF( ISIOPROCESSOR ) THEN
               DO 20 J = 1, JSIZE
                  DO 10 I = 1, ISIZE
                     READ( NIN, ERR = 50 ) TMP
                     WORK( I+( J-1 )*LDD ) = DBLE(TMP)
   10             CONTINUE
   20          CONTINUE
            END IF
*
            CALL PDGEADD( 'NoTrans', ISIZE, JSIZE, ALPHA, WORK, 1, 1,
     $                    DESCWORK, BETA, A, ISTART, JSTART, DESCA )
   30    CONTINUE
   40 CONTINUE
      IF( ISIOPROCESSOR ) THEN
         CLOSE ( NIN, ERR = 50 )
      END IF
      WORK( 1 ) = DESCA( MB_ )
      RETURN
   50 CONTINUE
      CALL PXERBLA( DESCA( CTXT_ ), 'PLAWRITE', 1 )
      WORK( 1 ) = DESCA( MB_ )
      RETURN
*
*     End of PDLAREAD
*
      END
