PROGRAM HEISENM
C
C     MULTI-CASE HEISENBERG UNCERTAINTY DEMO
C     WITH FOURIER-SPACE MOMENTUM AND FINITE-DIFFERENCE
C     MOMENTUM-OPERATOR CROSS-CHECK
C
C     MODEL
C     -----
C     PSI(X) ~ EXP(-(X-X0)**2/(4*SIGMA**2)) * EXP(I*K0*X)
C
C     THEN |PSI|**2 IS A GAUSSIAN WITH STD DEV SIGMA.
C
C     EXACT IDEAL CONTINUOUS RESULTS
C     ------------------------------
C     DELTA X    = SIGMA
C     <P>        = HBAR*K0
C     DELTA P    = HBAR/(2*SIGMA)
C     DELX*DELP  = HBAR/2
C
C     DATA DICTIONARY
C     ---------------
C     N        NUMBER OF X GRID POINTS
C     NCASE    NUMBER OF TEST CASES
C     SIGMAS   CASE LIST OF POSITION STD DEV
C     K0S      CASE LIST OF CENTRAL WAVENUMBER
C     X0S      CASE LIST OF POSITION CENTRE
C
      IMPLICIT NONE

      INTEGER N, NCASE
      PARAMETER (N=256)
      PARAMETER (NCASE=4)

      REAL SIGMAS(NCASE), K0S(NCASE), X0S(NCASE)
      INTEGER IC

      DATA SIGMAS /0.5E0, 1.0E0, 2.0E0, 1.0E0/
      DATA K0S    /0.0E0, 0.0E0, 0.0E0, 2.0E0/
      DATA X0S    /0.0E0, 0.0E0, 0.0E0, 0.0E0/

      WRITE(*,*) ' '
      WRITE(*,*) 'MULTI-CASE HEISENBERG UNCERTAINTY DEMO'
      WRITE(*,*) 'FFT AND FINITE-DIFFERENCE CROSS-CHECK'
      WRITE(*,*) ' '

      DO 100 IC = 1, NCASE
         CALL RUNCASE(N, IC, SIGMAS(IC), K0S(IC), X0S(IC))
  100 CONTINUE

      STOP
      END


      SUBROUTINE RUNCASE(N, IC, SIGMA, K0, X0)
      IMPLICIT NONE

      INTEGER N, IC
      INTEGER I, J

      REAL SIGMA, K0, X0
      REAL PI, HBAR, L, DX, DK, EPS
      REAL X(512), P(512)
      REAL REPSI(512), IMPSI(512)
      REAL REPHI(512), IMPHI(512)

      REAL AMP, EDGE
      REAL PROBX, PROBK
      REAL NORMX, NORMK
      REAL XMEAN, X2MEAN, VARX, DELX

      REAL PMFFT, P2FFT, VARFFT, DELPFFT
      REAL PMFD,  P2FD,  VARFD,  DELPFD

      REAL KVAL, THETA, CTH, STH
      REAL SUMRE, SUMIM
      REAL DARE, DAIM, D2ARE, D2AIM

      REAL EXX, EXPM, EXDP, EXPROD
      REAL ERRX, ERRPMF, ERRPFF
      REAL ERRDPF, ERRDPD, ERRUF, ERRUD

      PI   = 3.14159265358979E0
      HBAR = 1.0E0
      L    = 40.0E0
      DX   = L / FLOAT(N)
      DK   = 2.0E0 * PI / L
      EPS  = 1.0E-6

      IF (N .GT. 512) THEN
         WRITE(*,*) 'ERROR: N EXCEEDS ARRAY LIMIT 512'
         STOP
      ENDIF

C     BUILD X GRID AND UNNORMALIZED PSI(X)
      DO 100 I = 1, N
         X(I) = -0.5E0*L + FLOAT(I-1)*DX
         AMP = EXP( -((X(I)-X0)**2) / (4.0E0*SIGMA*SIGMA) )
         REPSI(I) = AMP * COS(K0*X(I))
         IMPSI(I) = AMP * SIN(K0*X(I))
  100 CONTINUE

C     RAW X NORM
      NORMX = 0.0E0
      DO 110 I = 1, N
         PROBX = REPSI(I)*REPSI(I) + IMPSI(I)*IMPSI(I)
         NORMX = NORMX + PROBX*DX
  110 CONTINUE

      IF (NORMX .LE. 0.0E0) THEN
         WRITE(*,*) 'ERROR: NONPOSITIVE X NORM'
         STOP
      ENDIF

      NORMX = SQRT(NORMX)

      DO 120 I = 1, N
         REPSI(I) = REPSI(I) / NORMX
         IMPSI(I) = IMPSI(I) / NORMX
  120 CONTINUE

C     BOUNDARY DIAGNOSTIC
      EDGE = 0.0E0
      AMP = SQRT(REPSI(1)*REPSI(1) + IMPSI(1)*IMPSI(1))
      IF (AMP .GT. EDGE) EDGE = AMP
      AMP = SQRT(REPSI(N)*REPSI(N) + IMPSI(N)*IMPSI(N))
      IF (AMP .GT. EDGE) EDGE = AMP

C     X MOMENTS
      XMEAN  = 0.0E0
      X2MEAN = 0.0E0
      DO 130 I = 1, N
         PROBX  = REPSI(I)*REPSI(I) + IMPSI(I)*IMPSI(I)
         XMEAN  = XMEAN  + X(I)      * PROBX * DX
         X2MEAN = X2MEAN + X(I)*X(I) * PROBX * DX
  130 CONTINUE

      VARX = X2MEAN - XMEAN*XMEAN
      IF (VARX .LT. 0.0E0 .AND. ABS(VARX) .LT. EPS) VARX = 0.0E0
      IF (VARX .LT. 0.0E0) THEN
         WRITE(*,*) 'ERROR: NEGATIVE X VARIANCE', VARX
         STOP
      ENDIF
      DELX = SQRT(VARX)

C     FOURIER MOMENTUM CALCULATION
      DO 200 J = 1, N
         KVAL = DK * FLOAT(J-1-N/2)
         P(J) = HBAR * KVAL
         SUMRE = 0.0E0
         SUMIM = 0.0E0

         DO 180 I = 1, N
            THETA = -KVAL * X(I)
            CTH   = COS(THETA)
            STH   = SIN(THETA)
            SUMRE = SUMRE + (REPSI(I)*CTH
     1             - IMPSI(I)*STH) * DX
            SUMIM = SUMIM + (REPSI(I)*STH
     1             + IMPSI(I)*CTH) * DX
  180    CONTINUE

         REPHI(J) = SUMRE / SQRT(2.0E0*PI)
         IMPHI(J) = SUMIM / SQRT(2.0E0*PI)
  200 CONTINUE

      NORMK = 0.0E0
      DO 210 J = 1, N
         PROBK = REPHI(J)*REPHI(J) + IMPHI(J)*IMPHI(J)
         NORMK = NORMK + PROBK*DK
  210 CONTINUE

      IF (NORMK .LE. 0.0E0) THEN
         WRITE(*,*) 'ERROR: NONPOSITIVE K NORM'
         STOP
      ENDIF

      NORMK = SQRT(NORMK)

      DO 220 J = 1, N
         REPHI(J) = REPHI(J) / NORMK
         IMPHI(J) = IMPHI(J) / NORMK
  220 CONTINUE

      PMFFT = 0.0E0
      P2FFT = 0.0E0
      DO 230 J = 1, N
         PROBK = REPHI(J)*REPHI(J) + IMPHI(J)*IMPHI(J)
         PMFFT = PMFFT + P(J)      * PROBK * DK
         P2FFT = P2FFT + P(J)*P(J) * PROBK * DK
  230 CONTINUE

      VARFFT = P2FFT - PMFFT*PMFFT
      IF (VARFFT .LT. 0.0E0 .AND. ABS(VARFFT) .LT. EPS)
     1    VARFFT = 0.0E0
      IF (VARFFT .LT. 0.0E0) THEN
         WRITE(*,*) 'ERROR: NEGATIVE FFT P VARIANCE', VARFFT
         STOP
      ENDIF
      DELPFFT = SQRT(VARFFT)

C     FINITE-DIFFERENCE MOMENTUM OPERATOR
C     <P>  = HBAR * INTEGRAL (A*DB/DX - B*DA/DX) DX
C     <P2> = -HBAR**2 * INTEGRAL (A*D2A/DX2 + B*D2B/DX2) DX
      PMFD = 0.0E0
      P2FD = 0.0E0
      DO 300 I = 2, N-1
         DARE = (REPSI(I+1) - REPSI(I-1)) / (2.0E0*DX)
         DAIM = (IMPSI(I+1) - IMPSI(I-1)) / (2.0E0*DX)

         D2ARE = (REPSI(I+1) - 2.0E0*REPSI(I)
     1          + REPSI(I-1)) / (DX*DX)
         D2AIM = (IMPSI(I+1) - 2.0E0*IMPSI(I)
     1          + IMPSI(I-1)) / (DX*DX)

         PMFD = PMFD + HBAR * (REPSI(I)*DAIM
     1        - IMPSI(I)*DARE) * DX

         P2FD = P2FD - HBAR*HBAR * (REPSI(I)*D2ARE
     1        + IMPSI(I)*D2AIM) * DX
  300 CONTINUE

      VARFD = P2FD - PMFD*PMFD
      IF (VARFD .LT. 0.0E0 .AND. ABS(VARFD) .LT. EPS)
     1    VARFD = 0.0E0
      IF (VARFD .LT. 0.0E0) THEN
         WRITE(*,*) 'ERROR: NEGATIVE FD P VARIANCE', VARFD
         STOP
      ENDIF
      DELPFD = SQRT(VARFD)

C     EXACT VALUES
      EXX    = SIGMA
      EXPM   = HBAR * K0
      EXDP   = HBAR / (2.0E0*SIGMA)
      EXPROD = 0.5E0 * HBAR

      ERRX   = ABS(DELX - EXX)
      ERRPMF = ABS(PMFFT - EXPM)
      ERRPFF = ABS(PMFD  - EXPM)
      ERRDPF = ABS(DELPFFT - EXDP)
      ERRDPD = ABS(DELPFD  - EXDP)
      ERRUF  = ABS(DELX*DELPFFT - EXPROD)
      ERRUD  = ABS(DELX*DELPFD  - EXPROD)

      WRITE(*,*) ' '
      WRITE(*,*) '----------------------------------------'
      WRITE(*,*) 'CASE ', IC
      WRITE(*,*) 'SIGMA                 = ', SIGMA
      WRITE(*,*) 'K0                    = ', K0
      WRITE(*,*) 'X0                    = ', X0
      WRITE(*,*) 'L                     = ', L
      WRITE(*,*) 'N                     = ', N
      WRITE(*,*) 'DX                    = ', DX
      WRITE(*,*) 'DK                    = ', DK
      WRITE(*,*) 'MAX BOUNDARY |PSI|    = ', EDGE

      WRITE(*,*) ' '
      WRITE(*,*) 'POSITION'
      WRITE(*,*) '<X>                   = ', XMEAN
      WRITE(*,*) 'DELTA X               = ', DELX
      WRITE(*,*) 'EXACT DELTA X         = ', EXX
      WRITE(*,*) 'ABS ERR DELTA X       = ', ERRX

      WRITE(*,*) ' '
      WRITE(*,*) 'MOMENTUM FROM FOURIER SPACE'
      WRITE(*,*) '<P> FFT               = ', PMFFT
      WRITE(*,*) 'EXACT <P>             = ', EXPM
      WRITE(*,*) 'ABS ERR <P> FFT       = ', ERRPMF
      WRITE(*,*) 'DELTA P FFT           = ', DELPFFT
      WRITE(*,*) 'EXACT DELTA P         = ', EXDP
      WRITE(*,*) 'ABS ERR DELTA P FFT   = ', ERRDPF
      WRITE(*,*) 'DELX*DELP FFT         = ', DELX*DELPFFT
      WRITE(*,*) 'ABS ERR PRODUCT FFT   = ', ERRUF

      WRITE(*,*) ' '
      WRITE(*,*) 'MOMENTUM FROM FINITE DIFFERENCES'
      WRITE(*,*) '<P> FD                = ', PMFD
      WRITE(*,*) 'EXACT <P>             = ', EXPM
      WRITE(*,*) 'ABS ERR <P> FD        = ', ERRPFF
      WRITE(*,*) 'DELTA P FD            = ', DELPFD
      WRITE(*,*) 'EXACT DELTA P         = ', EXDP
      WRITE(*,*) 'ABS ERR DELTA P FD    = ', ERRDPD
      WRITE(*,*) 'DELX*DELP FD          = ', DELX*DELPFD
      WRITE(*,*) 'ABS ERR PRODUCT FD    = ', ERRUD

      WRITE(*,*) ' '
      IF (DELX*DELPFFT .GE. 0.5E0*HBAR) THEN
         WRITE(*,*) 'FFT CHECK: DELX*DELP >= HBAR/2'
      ELSE
         WRITE(*,*) 'FFT CHECK: BELOW HBAR/2, NUMERICAL ERROR'
      ENDIF

      IF (DELX*DELPFD .GE. 0.5E0*HBAR) THEN
         WRITE(*,*) 'FD  CHECK: DELX*DELP >= HBAR/2'
      ELSE
         WRITE(*,*) 'FD  CHECK: BELOW HBAR/2, NUMERICAL ERROR'
      ENDIF

      RETURN
      END
~/downloads $
