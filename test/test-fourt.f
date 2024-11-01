      PROGRAM TESTFOURT

      DIMENSION DATA(8),WORK(50),NN(1)
      COMPLEX DATA
      DATA NN/8/
      
C     INPUT DATA https://rosettacode.org/wiki/Fast_Fourier_transform      
      DATA(1)=COMPLEX(1.,0.)
      DATA(2)=COMPLEX(1.,0.)
      DATA(3)=COMPLEX(1.,0.)
      DATA(4)=COMPLEX(1.,0.)
      DATA(5)=COMPLEX(0.,0.)
      DATA(6)=COMPLEX(0.,0.)
      DATA(7)=COMPLEX(0.,0.)
      DATA(8)=COMPLEX(0.,0.)

      WRITE  (*,*) "DATA"
      DO 1 I=1,8
      PRINT 10,DATA(I)
10    FORMAT (2H (,F6.3,3H , ,F6.3,1H))  
1     CONTINUE

      CALL FOURT(DATA,NN,1,-1,1,WORK)

C     EXPECT OUTPUT
C     4,0
C     1,-2.414
C     0,0
C     1,-0.414
C     0,0
C     1,0.414
C     0,0
C     1,2.414      
      WRITE (*,*) "FFT (FOURT.F) OF DATA"
      DO 2 I=1,8
      PRINT 10,DATA(I)
2     CONTINUE
      
      CALL FOURT(DATA,NN,1,1,1,WORK)

      WRITE (*,*) "TRANSFORM BACK"
      DO 3 I=1,8
      PRINT 10,DATA(I)
3     CONTINUE
      STOP
      END
  
