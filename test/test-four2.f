      PROGRAM TESTFOUR2

      DIMENSION DATA(8),NN(1)
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
      WRITE (*,*) DATA(I)
1     CONTINUE

      CALL FOUR2(DATA,NN,1,-1)

C     EXPECT OUTPUT
C     4,0
C     1,-2.414
C     0,0
C     1,-0.414
C     0,0
C     1,0.414
C     0,0
C     1,2.414      
      WRITE (*,*) "FFT (FOUR2.F) OF DATA"
      DO 2 I=1,8
      WRITE (*,*) DATA(I)
2     CONTINUE
      
      CALL FOUR2(DATA,NN,1,1)

      WRITE (*,*) "TRANSFORM BACK"
      DO 3 I=1,8
      WRITE (*,*) DATA(I)
3     CONTINUE
      STOP
      END