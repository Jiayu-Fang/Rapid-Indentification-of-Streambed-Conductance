PROGRAM MAIN
	IMPLICIT NONE
	REAL (KIND = 8) :: KA,SS,D,K,C,DA,L
	REAL (KIND = 8) :: ALPHA,A,DX,T,DT,F
	REAL (KIND = 8) :: X(50001)
	REAL (KIND = 8) :: H1(50001),H(50001)
	REAL (KIND = 8) :: AP(50001),APE(50001),APW(50001),SOUR(50001)
	REAL (KIND = 8) :: AP_N(50001),APW_N(50001)
	REAL (KIND = 8) :: RESL
	INTEGER :: N,IK

	L = 500.0; KA = 1.0; SS = 1.0E-4; D = 1.0; K = 0.001; C = 0.1
	DA = 50.0; ALPHA = KA*DA/(SS*DA); A = D*KA/K;
	DX = 0.01; X(1) = 0.0
	DO IK = 2,50001
		X(IK) = X(IK-1)+DX;
	END DO 
	H1 = 0.0; H = 0.0
	T = 0.0; DT = 0.01
	N = 1
	DO WHILE (T<2.01)
		DO IK = 2,50000
			AP(IK) = (2.0*KA/DX)+1.0/DT*DX*SS
			APW(IK) = -KA/DX
			APE(IK) = -KA/DX
			SOUR(IK) = H1(IK)/DT*DX*SS;
		END DO
		IK = 50001
		AP(IK) = -KA/DX+2.0*KA/DX+DX/DT*SS
		APW(IK) = -KA/DX
		APE(IK) = 0.0
		SOUR(IK) = H1(IK)/DT*DX*SS
		IK = 1
		F = 1.0
		AP(IK) = K/D+KA/DX+DX/DT*SS
		APW(IK) = 0.0
		APE(IK) = -KA/DX
		SOUR(IK) = K/D*F+H1(IK)/DT*DX*SS
		T = T+DT
	! LU DECOMPOSITION
		AP_N(1) = AP(1)
		DO IK = 2,50001
			APW_N(IK) = APW(IK)/AP_N(IK-1)
			AP_N(IK) = AP(IK) - APW_N(IK)*APE(IK-1)
		END DO
		IK = 1
		H(IK) = SOUR(IK)
		DO IK = 2,50001
			H(IK) = SOUR(IK)-APW_N(IK)*H(IK-1)
		END DO
		IK = 50001
		H(IK) = H(IK)/AP_N(IK)
		DO IK = 50000,1,-1
			H(IK) = (H(IK) - H(IK+1)*APE(IK))/AP_N(IK)
		END DO 
		DO IK = 1,50001
			H1(IK) = H(IK)
		END DO
		RESL = 0.0
		DO IK = 2,50000
			RESL = MAX(ABS(SOUR(IK) - APE(IK)*H(IK+1) - AP(IK)*H(IK) - APW(IK)*H(IK-1)),RESL)
		END DO
		WRITE(*,*) 'RESL = ', RESL
		IK = 1
		WRITE(*,*) ABS(SOUR(IK)-AP(IK)*H(IK)-APE(IK)*H(IK+1))
		IK = 3
		WRITE (*,*) ABS(SOUR(IK) - APE(IK)*H(IK+1) - AP(IK)*H(IK) - APW(IK)*H(IK-1))
	END DO
	OPEN (1, FILE = 'RESULTS.DAT')
	DO IK = 1,50001
		WRITE (1,'(10F20.10)') H(IK),APW(IK),AP(IK),APE(IK),SOUR(IK),APW_N(IK),AP_N(IK),APE(IK)
	END DO
	CLOSE (1)
END PROGRAM
