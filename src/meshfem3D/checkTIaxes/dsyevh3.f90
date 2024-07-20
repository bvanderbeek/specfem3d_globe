      SUBROUTINE DSYEVH3(A, Q, W)
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

      DOUBLE PRECISION EPS
      PARAMETER        ( EPS = 2.2204460492503131D-16 )

      DOUBLE PRECISION NORM
      DOUBLE PRECISION ERROR
      DOUBLE PRECISION T, U
      INTEGER          J

      EXTERNAL         DSYEVC3, DSYEVQ3

      CALL DSYEVC3(A, W)


      T       = MAX(ABS(W(1)), ABS(W(2)), ABS(W(3)))
      U       = MAX(T, T**2)
      ERROR   = 256.0D0 * EPS * U**2
      Q(1, 2) = A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2)
      Q(2, 2) = A(1, 3) * A(1, 2) - A(2, 3) * A(1, 1)
      Q(3, 2) = A(1, 2)**2

      Q(1, 1) = Q(1, 2) + A(1, 3) * W(1)
      Q(2, 1) = Q(2, 2) + A(2, 3) * W(1)
      Q(3, 1) = (A(1,1) - W(1)) * (A(2,2) - W(1)) - Q(3,2)
      NORM    = Q(1, 1)**2 + Q(2, 1)**2 + Q(3, 1)**2

      IF (NORM .LE. ERROR) THEN
        CALL DSYEVQ3(A, Q, W)
        RETURN
      ELSE
        NORM = SQRT(1.0D0 / NORM)
        DO 20, J = 1, 3
          Q(J, 1) = Q(J, 1) * NORM
   20   CONTINUE
      END IF
 
      Q(1, 2) = Q(1, 2) + A(1, 3) * W(2)
      Q(2, 2) = Q(2, 2) + A(2, 3) * W(2)
      Q(3, 2) = (A(1,1) - W(2)) * (A(2,2) - W(2)) - Q(3, 2)
      NORM    = Q(1, 2)**2 + Q(2, 2)**2 + Q(3, 2)**2
      IF (NORM .LE. ERROR) THEN
        CALL DSYEVQ3(A, Q, W)
        RETURN
      ELSE
        NORM = SQRT(1.0D0 / NORM)
        DO 40, J = 1, 3 
          Q(J, 2) = Q(J, 2) * NORM
   40   CONTINUE
      END IF

      Q(1, 3) = Q(2, 1) * Q(3, 2) - Q(3, 1) * Q(2, 2)
      Q(2, 3) = Q(3, 1) * Q(1, 2) - Q(1, 1) * Q(3, 2)
      Q(3, 3) = Q(1, 1) * Q(2, 2) - Q(2, 1) * Q(1, 2)

      END SUBROUTINE

