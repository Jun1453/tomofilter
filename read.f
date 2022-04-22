      PROGRAM READ
      INTEGER*4 I, J, M, N
      REAL*8 TMP(500000)
      OPEN( 9999, FILE='kernel_input', STATUS='UNKNOWN',
     $      FORM='UNFORMATTED' )
      READ( 9999 ) M
      READ( 9999 ) N

      OPEN( 8888, FILE='double_input', STATUS='REPLACE',
     $      ACCESS='STREAM', FORM='UNFORMATTED' )
!     $       FORM='FORMATTED' )
      WRITE( 8888 ) M
      WRITE( 8888 ) N

      DO 20 I = 1, N
         READ( 9999 ) ( TMP( J ) , J = 1, M)
         DO 30 J = 1, M
            WRITE( 8888 ) TMP( J )
30    CONTINUE
20    CONTINUE
!     DO 30 I=1,N
!        
      END
