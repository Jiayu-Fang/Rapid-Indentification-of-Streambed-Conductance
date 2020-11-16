PROGRAM MAIN
	IMPLICIT NONE 
	INTEGER :: S_NUM
	REAL (KIND = 8) :: ALPHA,A,KA,DA,S,DISTANCE,DD
	REAL (KIND = 8) :: POW,POW_S
	REAL (KIND = 8) :: A1,A2
	REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) :: HEAD,T_HEAD,HEAD_GW
	REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) :: K_HEAD
	REAL (KIND = 8) :: HEAD_CAL1,HEAD_CAL2
	REAL (KIND = 8) :: DT1,DT2
	REAL (KIND = 8) :: DT,X
	REAL (KIND = 8) :: T1,T2
	REAL (KIND = 8) :: T_TOL
	CHARACTER (LEN = 256) :: PATH,NAME1,NAME2,NAME3
    INTEGER :: I,J,K,JJ
    INTEGER :: N,K1,K2
    INTEGER :: J1,J2
	INTEGER :: I_STEP,IK
    REAL (KIND = 8) :: T_1,T_2,H1,H2
    REAL (KIND = 8) :: T_11,T_22
    REAL (KIND = 8) :: FK1,FK2,PN1,PN2
	REAL (KIND = 8) :: SSE1,SSE2,SSE
	REAL (KIND = 8) :: SS1,SS2,SS3
	REAL (KIND = 8) :: STEP_SIZE,STEP_TOLERANCE1,STEP_TOLERANCE2
	REAL (KIND = 8) :: STEP_SIZE1,STEP_SIZE_DEF
	REAL (KIND = 8) :: XX,X_DX,DF
	REAL (KIND = 8) :: DFF
	REAL (KIND = 8) :: XX1,XX2,XX3
	REAL (KIND = 8) :: X_DX1,DF1
	REAL (KIND = 8) :: Y1,Y2,X1,X2
	REAL (KIND = 8) :: C,CT
	REAL (KIND = 8) :: CT1,CT2,CT3
	REAL (KIND = 8) :: SSEX1,SSEX2

	CALL getcwd(PATH)
	NAME1 = trim(PATH)//'\INPUT.DAT'
	NAME2 = trim(PATH)//'\MEASUREMENT.DAT'
	NAME3 = trim(PATH)//'\CALCULATION.DAT'
	OPEN (1, FILE = trim(NAME1))
	OPEN (2, FILE = trim(NAME2))
	OPEN (3, FILE = trim(NAME3))

	READ (1,*) DISTANCE,KA,DA,S,DD ! Distance, Aquifer Hydraulic Conductivity,Aquifer Thickness, Storativity
	READ (1,*) DT
	READ (1,*) POW,POW_S
	READ (1,*) STEP_SIZE_DEF,STEP_TOLERANCE1,STEP_TOLERANCE2
	READ (1,*) STEP_SIZE1,C
	READ (1,*) XX
	STEP_SIZE = STEP_SIZE_DEF
	CLOSE (1)

	ALPHA = KA*DA/S

	READ (2,*) S_NUM
	ALLOCATE (T_HEAD(S_NUM),HEAD(S_NUM),HEAD_GW(S_NUM))
	DO I = 1,S_NUM
		READ (2,*) T_HEAD(I),HEAD(I),HEAD_GW(I)
	END DO

	CLOSE (2)
	
	JJ = FLOOR((T_HEAD(S_NUM)-POW)*1.0/POW_S+0.5)
	ALLOCATE (K_HEAD(JJ))
	K_HEAD(1) = XX
	K1 = 1
	K2 = 1
	T_TOL = 0.5*POW_S
	DO I = 1,JJ
		T1 = (I-1)*POW_S
		T2 = T1+POW
		DO WHILE (T_HEAD(K1).LT.T1)
			K1 = K1+1
		END DO
		K2 = K1
		DO WHILE (K2.LE.S_NUM.AND.T_HEAD(K2).LT.T2)
			K2 = K2+1
		END DO
		IF (I.NE.1) THEN 
			K_HEAD(I) = K_HEAD(I-1)
		END IF
		STEP_SIZE = STEP_SIZE_DEF
		X_DX = 1.0
		X = K_HEAD(I)
		I_STEP = 0
		SSE1 = 1.0
		DF = 1.0
		CT = 1.0
		DFF = 1.0
		DO WHILE (ABS(CT).GE.STEP_TOLERANCE2.AND.ABS(DFF).GT.STEP_TOLERANCE2.AND.I_STEP.LE.20)
			SS1 = 0.0; SS2 = 0.0; SS3 = 0.0
			DO IK = 1,3
                X1 = X
                X2 = X+STEP_SIZE1
                SSE1 = 0.0
                SSE2 = 0.0
                SSE = 0.0
                DF1 = 0.0
                CT = 0.0
                DO K = K1,K2
                    H1 = 0.0
                    H2 = 0.0
                    A1 = DD*KA/X1
                    A2 = DD*KA/X2
                    T_2 = T_HEAD(K)
                    T_1 = T_HEAD(K) - DT
                    J1 = K
                    J2 = K
                    DO WHILE ((T_1.GE.T_HEAD(1)).AND.((DISTANCE/2.0/SQRT(ALPHA*T_1)).LT.8.0))
                        J = J1
                        DO WHILE (T_2.LT.T_HEAD(J).AND.J.GT.1)
                            J = J - 1
                        END DO
                        J1 = J
                        FK2 = (HEAD(J+1)-HEAD(J))/(T_HEAD(J+1)-T_HEAD(J))
                        J = J2
                        DO WHILE (T_1.LT.T_HEAD(J).AND.J.GT.1)
                            J = J - 1
                        END DO
                        J2 = J
                        FK1 = (HEAD(J+1)-HEAD(J))/(T_HEAD(J+1)-T_HEAD(J))
                        T_11 = T_HEAD(K) - T_1
                        T_22 = T_HEAD(K) - T_2
                        CALL H_CAL (PN1,DISTANCE,ALPHA,A1,T_11,K)
                        CALL H_CAL (PN2,DISTANCE,ALPHA,A1,T_22,K)
                        H1 = H1+DT/6.0*(FK1*PN1+FK2*PN2+2.0*FK1*PN2+2.0*FK2*PN1)
                        CALL H_CAL (PN1,DISTANCE,ALPHA,A2,T_11,K)
                        CALL H_CAL (PN2,DISTANCE,ALPHA,A2,T_22,K)
                        H2 = H2+DT/6.0*(FK1*PN1+FK2*PN2+2.0*FK1*PN2+2.0*FK2*PN1)
                        T_1 = T_1 - DT
                        T_2 = T_2 - DT
                    END DO
                    HEAD_CAL1 = HEAD(1) + H1
                    HEAD_CAL2 = HEAD(1) + H2
                    DF1 = (HEAD_CAL2 - HEAD_CAL1)/(X2 - X1)
                    SSE1 = SSE1 + DF1*DF1
                    SSE2 = SSE2 + DF1*(HEAD_GW(K) - HEAD_CAL1)
                    CT = CT + ((HEAD_CAL2 - HEAD_GW(K))**2 - (HEAD_CAL1 - HEAD_GW(K))**2)/(X2 - X1)
                    SSE = SSE + (HEAD_GW(K) - HEAD_CAL1) * (HEAD_GW(K) - HEAD_CAL1)
                END DO
                IF (IK.EQ.1) THEN
                    SS1 = SSE
                    CT1 = CT/(K2 - K1)
                    XX1 = X
                    X = X + SSE2/(SSE1 + SSE1*STEP_SIZE)
                    X = MAX(X,1.0E-9)
                    SSEX1 = SSE1; SSEX2 = SSE2
                END IF
                IF (IK.EQ.2) THEN
                    SS2 = SSE
                    CT2 = CT/(K2 - K1)
                    XX2 = X
                    X = X + SSE2/(SSE1 + SSE1*STEP_SIZE/C)
                    X = MAX(X,1.0E-9)
                END IF
                IF (IK.EQ.3) THEN 
                    SS3 = SSE
                    CT3 = CT/(K2 - K1)
                    XX3 = X
                END IF
            END DO
!			WRITE (*,*) SS1,SS2,SS3
!			WRITE (*,*) XX1,XX2,XX3
			IF (SS1.GE.SS2.AND.SS2.GE.SS3) THEN 
                STEP_SIZE = STEP_SIZE/C
                CT = CT3
                X = XX3
				DFF = ABS(XX3 - XX1)
                I_STEP = I_STEP + 1
				SSE = SS3
            END IF
            IF (SS1.GE.SS2.AND.SS3.GE.SS2) THEN 
                STEP_SIZE = STEP_SIZE
                CT = CT2
                X = XX2
				DFF = ABS(XX2 - XX1)
                I_STEP = I_STEP + 1
				SSE = SS2
            END IF
            IF (SS1.LE.SS2.AND.SS1.LE.SS3) THEN 
                STEP_SIZE = STEP_SIZE * C
                X = XX1 + SSEX2/(SSEX1 + SSEX1*STEP_SIZE)
				DFF = 1.0
                DO WHILE (SS2.GE.SS1.AND.ABS(DFF).GT.STEP_TOLERANCE2)
                    X1 = X
                    X2 = X + STEP_SIZE1
                    SSE1 = 0.0; SSE2 = 0.0
                    SSE = 0.0
                    SS2 = 0.0
                    DF1 = 0.0
                    CT = 0.0
                    DO K = K1,K2
                        H1 = 0.0; H2 = 0.0
                        A1 = DD*KA/X1; A2 = DD*KA/X2
                        T_2 = T_HEAD(K)
                        T_1 = T_HEAD(K) - DT
                        J1 = K
                        J2 = K
                        DO WHILE ((T_1.GE.T_HEAD(1)).AND.((DISTANCE/2.0/SQRT(ALPHA*T_1)).LT.8.0))
                            J = J1
                            DO WHILE (T_2.LT.T_HEAD(J).AND.J.GT.1)
                                J = J - 1
                            END DO
                            J1 = J
                            FK2 = (HEAD(J+1)-HEAD(J))/(T_HEAD(J+1)-T_HEAD(J))
                            J = J2
                            DO WHILE (T_1.LT.T_HEAD(J).AND.J.GT.1)
                                J = J - 1
                            END DO
                            J2 = J
                            FK1 = (HEAD(J+1)-HEAD(J))/(T_HEAD(J+1)-T_HEAD(J))
                            T_11 = T_HEAD(K) - T_1
                            T_22 = T_HEAD(K) - T_2
                            CALL H_CAL (PN1,DISTANCE,ALPHA,A1,T_11,K)
                            CALL H_CAL (PN2,DISTANCE,ALPHA,A1,T_22,K)
                            H1 = H1+DT/6.0*(FK1*PN1+FK2*PN2+2.0*FK1*PN2+2.0*FK2*PN1)
                            CALL H_CAL (PN1,DISTANCE,ALPHA,A2,T_11,K)
                            CALL H_CAL (PN2,DISTANCE,ALPHA,A2,T_22,K)
                            H2 = H2+DT/6.0*(FK1*PN1+FK2*PN2+2.0*FK1*PN2+2.0*FK2*PN1)
                            T_1 = T_1 - DT
                            T_2 = T_2 - DT
                        END DO
                        HEAD_CAL1 = HEAD(1) + H1
                        HEAD_CAL2 = HEAD(1) + H2
                        DF1 = (HEAD_CAL2 - HEAD_CAL1)/(X2 - X1)
                        SSE1 = SSE1 + DF1*DF1
                        SSE2 = SSE2 + DF1*(HEAD_GW(K) - HEAD_CAL1)
                        CT = CT + ((HEAD_CAL2 - HEAD_GW(K))**2 - (HEAD_CAL1 - HEAD_GW(K))**2)/(X2 - X1)
                        SSE = SSE + (HEAD_GW(K) - HEAD_CAL1) * (HEAD_GW(K) - HEAD_CAL1)
                    END DO
                    SS2 = SSE
                    STEP_SIZE = STEP_SIZE * C
                    XX2 = X + SSE2/(SSE1 + SSE1*STEP_SIZE)
					DFF = ABS(X - XX2)
                    X = MAX(XX2,1.0E-9)
                    CT2 = CT/(K2 - K1)
                END DO
                X  = X - SSE2/(SSE1 + SSE1*STEP_SIZE)
                STEP_SIZE = STEP_SIZE/C
                X = MAX(X,1.0E-9)
                CT = CT2
                I_STEP = I_STEP + 1
            END IF
!			WRITE (*,*) I_STEP
!			WRITE (*,*) X,STEP_SIZE
!			WRITE (*,*) SS1,SS2,SS3
!			X_DX = DF*(X2-X1)
!			WRITE (*,*) 'I = ', I
!			WRITE (*,*) X1, X2
!			WRITE (*,*) SSE1, SSE2
!			WRITE (*,*) X_DX
		END DO	
		K_HEAD(I) = X
		WRITE (3,*) T_TOL,K_HEAD(I)
		WRITE (*,*) T_TOL,K_HEAD(I),SSE
!		WRITE (*,*) CT,I_STEP
		T_TOL = T_TOL+POW_S
	END DO
	CLOSE (3)		

END PROGRAM 

SUBROUTINE H_CAL(PN,DISTANCE,ALPHA,A,T,K)
    IMPLICIT NONE
    REAL (KIND = 8) :: PN,DISTANCE,ALPHA,A,T
    REAL (KIND = 8) :: X1
	REAL (KIND = 8) :: X
	INTEGER :: K

    X1 = DISTANCE/2.0/SQRT(ALPHA*T)
	
 !   PN = ERFC(X1)-MIN(1.7E+38,EXP(DISTANCE/A+ALPHA*T/A/A))*ERFC(X1+SQRT(ALPHA*T)/A)
	IF ((DISTANCE/A+ALPHA*T/A/A).LT.709.0) THEN 
		PN = ERFC(X1)-EXP(DISTANCE/A+ALPHA*T/A/A)*ERFC(X1+SQRT(ALPHA*T)/A)
	ELSE 
		PN = ERFC(X1)
	END IF

!	PN = ERFC(X1)-EXP(MIN(DISTANCE/A+ALPHA*T/A/A,MAXEXPONENT(X)))*ERFC(X1+SQRT(ALPHA*T)/A)

END SUBROUTINE
