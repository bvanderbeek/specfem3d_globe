      SUBROUTINE DSYEVQ3(A, Q, W)
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

      INTEGER          N
      PARAMETER        ( N = 3 )

      DOUBLE PRECISION E(3)
      DOUBLE PRECISION G, R, P, F, B, S, C, T
      INTEGER          NITER
      INTEGER          L, M, I, J, K

      EXTERNAL         DSYTRD3
      
      CALL DSYTRD3(A, Q, W, E)

      DO 10 L = 1, N-1
        NITER = 0

        DO 11 I = 1, 50
          DO 20 M = L, N-1
            G = ABS(W(M)) + ABS(W(M+1))
            IF (ABS(E(M)) + G .EQ. G) THEN
              GO TO 30
            END IF
   20     CONTINUE
   30     IF (M .EQ. L) THEN
            GO TO 10
          END IF

          NITER = NITER + 1
          IF (NITER >= 30) THEN
            PRINT *, 'DSYEVQ3: No convergence.'
            RETURN
          END IF

          G = (W(L+1) - W(L)) / (2.0D0 * E(L))
          R = SQRT(1.0D0 + G**2)
          IF (G .GE. 0.0D0) THEN
            G = W(M) - W(L) + E(L)/(G + R)
          ELSE
            G = W(M) - W(L) + E(L)/(G - R)
          END IF

          S = 1.0D0
          C = 1.0D0
          P = 0.0D0
          DO 40 J = M - 1, L, -1
            F = S * E(J)
            B = C * E(J)
            IF (ABS(F) .GT. ABS(G)) THEN
              C      = G / F
              R      = SQRT(1.0D0 + C**2)
              E(J+1) = F * R
              S      = 1.0D0 / R
              C      = C * S
            ELSE
              S      = F / G
              R      = SQRT(1.0D0 + S**2)
              E(J+1) = G * R
              C      = 1.0D0 / R
              S      = S * C
            END IF

            G      = W(J+1) - P
            R      = (W(J) - G) * S + 2.0D0 * C * B
            P      = S * R
            W(J+1) = G + P
            G      = C * R - B

            DO 50 K = 1, N
              T         = Q(K, J+1)
              Q(K, J+1) = S * Q(K, J) + C * T
              Q(K, J)   = C * Q(K, J) - S * T
   50       CONTINUE
   40     CONTINUE
          W(L) = W(L) - P
          E(L) = G
          E(M) = 0.0D0
   11   CONTINUE
   10 CONTINUE
  
      END SUBROUTINE

