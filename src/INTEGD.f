c version: Feb 23 1992
C INTEG. THIS SUBROUTINES DOES INVERSE DIFFERENCING OR INTEGRATION
C   OF TIME SERIES.
C
C-----------------PROGRAM AUTHOR: A.I. MCLEOD-------------------------
C-----------------DATE: SEPTEMBER 1977--------------------------------
C-----------------COPYRIGHT BY A.I. MCLEOD----------------------------
C
C**********************************************************************
      SUBROUTINE INTEGD(Z, NCAP, N, ID, IDS, ISEA, ZINIT, IDCAP)
C**********************************************************************
C
C INPUT PARAMETERS:
C   Z - VECTOR OF DIMENSION NCAP. Z(1),...,Z(N) CONTAINS THE TIME
C       SERIES TO BE INTEGRATED.
C   NCAP - LENGTH OF INTEGRATED TIME SERIES. NCAP=N+IDCAP=N+ID+ISEA*IDS
C   N - LENGTH OF DIFFERENCED TIME SERIES
C   ID - DEGREE OF NONSEASONAL DIFFERENCING OPERATOR
C   IDS - DEGREE OF SEASONAL DIFFERENCING OPERATOR
C   ISEA - SEASONAL PERIOD
C   ZINIT - VECTOR OF LENGTH IDCAP CONTAINING THE INITIAL VALUES
C          OF THE INTEGRATED TIME SERIES.
C   IDCAP - DIMENSION OF ZINIT. IDCAP=ID+ISEA*IDS
C
C OUTPUT PARAMETERS
C   Z - VECTOR OF DIMENSION NCAP CONTAINING THE INTEGRATED TIME
C       SERIES. Z(1),...,Z(IDCAP) IS IDENTICAL TO ZINIT(1),...,ZINIT(IDC
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DIMENSION Z(NCAP), ZINIT(IDCAP)
      IF( IDCAP .EQ. 0 ) RETURN
      DO I = 1, N
        II = N + 1 - I
        III = NCAP + 1 - I
        Z(III) = Z(II)
      ENDDO
      DO I = 1, IDCAP
        Z(I) = ZINIT(I)
      ENDDO
      IDP1 = IDCAP + 1
      IDP2 = IDCAP + 2
      IF( ID .NE. 0 ) THEN
        DO IN = 1, ID
          IDMIN = ID - IN
          IA = IDCAP
          DO I = 1, IDCAP
            ZINIT(I) = Z(I)
          ENDDO
          IF( IDMIN .NE. 0 ) THEN
            DO I = 1, IDMIN
              IA = IA - 1
              DO J = 1, IA
                JP1 = J + 1
                ZINIT(J) = ZINIT(JP1) - ZINIT(J)
              ENDDO
            ENDDO
          ENDIF
          IF( .NOT.((ISEA .EQ. 0) .OR. (IDS .EQ. 0)) ) THEN
            IA = IDCAP
            DO I = 1, IDS
              IA = IA - ISEA
              DO J = 1, IA
                INPSEA = IN + ISEA + J - 1
                INPJ = IN + J - 1
                ZINIT(INPJ) = ZINIT(INPSEA) - ZINIT(INPJ)
              ENDDO
            ENDDO
          ENDIF
          Z(IDP1) = Z(IDP1) + ZINIT(IN)
          DO I = IDP2, NCAP
            IM1 = I - 1
            Z(I) = Z(I) + Z(IM1)
          ENDDO
        ENDDO
      ENDIF
      IF( (ISEA .EQ. 0) .OR. (IDS .EQ. 0) ) RETURN
      ISPS = IDP1 + ISEA
      IWS = ID + 1
      DO IN = 1, IDS
        IDMIN = IDS - IN
        IA = IDS*ISEA
        DO I = 1, IDCAP
          ZINIT(I) = Z(I)
        ENDDO
        IF( IDMIN .NE. 0 ) THEN
          DO I = 1, IDMIN
            IA = IA - ISEA
            DO J = 1, IA
              JPSEA = IWS + J + ISEA - 1
              JPIWS = IWS + J - 1
              ZINIT(JPIWS) = ZINIT(JPSEA) - ZINIT(JPIWS)
            ENDDO
          ENDDO
        ENDIF
        DO I = 1, ISEA
          IWSPI = IWS + I - 1
          IDPI = IDP1 + I - 1
          Z(IDPI) = Z(IDPI) + ZINIT(IWSPI)
        ENDDO
        DO I = ISPS, NCAP
          IMSEA = I - ISEA
          Z(I) = Z(I) + Z(IMSEA)
        ENDDO
        IWS = IWS + ISEA
      ENDDO
      RETURN
      END
