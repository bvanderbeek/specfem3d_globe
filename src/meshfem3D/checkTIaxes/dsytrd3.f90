      SUBROUTINE DSYTRD3(A, Q, D, E)
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION D(3)
      DOUBLE PRECISION E(2)

      INTEGER          N
      PARAMETER        ( N = 3 )

      DOUBLE PRECISION U(N), P(N)
      DOUBLE PRECISION OMEGA, F
      DOUBLE PRECISION K, H, G
      INTEGER          I, J

      DO 10 I = 1, N
        Q(I,I) = 1.0D0
        DO 11, J = 1, I-1
          Q(I, J) = 0.0D0
          Q(J, I) = 0.0D0
   11   CONTINUE
   10 CONTINUE

      H = A(1,2)**2 + A(1,3)**2
      IF (A(1,2) .GT. 0.0D0) THEN
        G = -SQRT(H)
      ELSE
        G = SQRT(H)
      END IF
      E(1)  = G
      F     = G * A(1,2)
      U(2)  = A(1,2) - G
      U(3)  = A(1,3)

      OMEGA = H - F
      IF (OMEGA > 0.0D0) THEN
        OMEGA = 1.0D0 / OMEGA
        K     = 0.0D0
        DO 20 I = 2, N
          F    = A(2,I)*U(2) + A(I,3)*U(3)
          P(I) = OMEGA * F
          K    = K + U(I) * F
  20    CONTINUE
        K = 0.5D0 * K * OMEGA**2

        DO 30 I = 2, N
          P(I) = P(I) - K * U(I)
  30    CONTINUE

        D(1) = A(1,1)
        D(2) = A(2,2) - 2.0D0 * P(2) * U(2)
        D(3) = A(3,3) - 2.0D0 * P(3) * U(3)

        DO 40, J = 2, N
          F = OMEGA * U(J)
          DO 41 I = 2, N
            Q(I,J) = Q(I,J) - F * U(I)
   41     CONTINUE
   40   CONTINUE
            
        E(2) = A(2, 3) - P(2) * U(3) - U(2) * P(3)
      ELSE
        DO 50 I = 1, N
          D(I) = A(I, I)
  50    CONTINUE
        E(2) = A(2, 3)
      END IF
      
      END SUBROUTINE

