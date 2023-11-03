PROGRAM loadmodl
  USE mar_data
  USE calc_mod
  INTEGER I
  REAL (OP) :: TEST,t1,t2

  TEST=1.0_OP
  ETA=EPSILON(TEST)
  TOL=TINY(TEST)/ETA
  NLAY=0
  NDAT=0
  TNT=0
  TX_DESIG=''
  MAXNT=0
  DO I=1,100
     UN_NO(I)=.FALSE.
  ENDDO

 DIPOLE_PARTS=20

  CALL LOAD(2)
! calc is embedded in likelihood computation code (see below)
  call cpu_time(t1)
  call calc
  call cpu_time(t2)
	   write(6,*) 'computation time [s]', t2-t1
  call save(2)


 IF (ALLOCATED(TNESS)) DEALLOCATE(TNESS,RES,ANIS)
  IF (ALLOCATED(SYS_D)) DEALLOCATE(SYS_D,SYS_T,SYS_TYPE,SYS_PER,SYS_N)
  IF (ALLOCATED(NT)) THEN
     DEALLOCATE(NT,X,Y,Z,SYS_RESP,SYS_NAME,LFCAL,DELAY,ANG1,IMP)
     DEALLOCATE(CALF,ANG2,TYPE,CHI,TIME,VOLT,ERROR,CVOLT,STEP_OFF,DERIVATE)
  ENDIF
  IF (ALLOCATED(WTS)) DEALLOCATE(WTS)
  IF (ALLOCATED(LFIXR)) DEALLOCATE(LFIXR,LFIXA,LFIXT)

END PROGRAM loadmodl

!=========================================================================================

      SUBROUTINE forward2(model,dout,dobs,sdev,Eout,npar,ndata,nrec,CFR,dc)

!=========================================================================================

      USE mar_data
	  integer npar,ndata,nrec,i,k,j,nsdev
      REAL (OP) :: model(npar),dout(ndata),dobs(ndata),dcdout(nrec),CFR(nrec),dc(nrec),sdev(ndata)
      REAL (OP) Eout,Eout1

      do i=1,nlay-1
         tness(i) = model((i-1)*2+1)
          res(i)   = 10**model(i*2)
      enddo
      res(nlay) = 10**model(2*nlay-1)  
      
      
	     if (.not. lfixwr) w_res = model(2*nlay) 
		 	
! --------------------------------------
! now the time is shifted about dt and new data is calculated for that time and compared with the observed data that has been calculated for a different time
         if (.not. lfixwd) then
		     do j=1,nrec
               do i=1,nt(j)
                  time(j,i)=time(j,i)+model(npar-nrec+j)
	!			  write(6,*) 'i,j,time:',i,j,time(j,i)
               enddo
             enddo
	     endif
!	write(6,*) 'before calc'
		 ! ----------
      call calc
	     ! ----------
  !  write(6,*) 'after calc'
      k = 1
      do i=1,nrec
         do j=1,nt(i)
            dout(k) = cvolt(i,j)
            k = k+1
         enddo
      enddo
	  
	  k = 0
	     do i=1,nrec
		    dcdout(i) = 0.
			CFR(i) = 0.
	      do j=1,nt(i)
	         dcdout(i) = dcdout(i)+(dout(k+j)/sdev(k+j))**2
          enddo

          do j=1,nt(i)
             CFR(i) = CFR(i)+(dobs(k+j)*dout(k+j))/(sdev(k+j)**2*dcdout(i))
          enddo
		  k=sum(nt(1:i))
	     enddo

      if (LFCAL(1)) then
      Eout = 0.
	  k=1
       do i=1,nrec
         do j=1,nt(i)
            Eout=Eout+(volt(i,j)-cvolt(i,j))**2/(2.*sdev(k)**2)
            k = k+1
         enddo
      enddo
	  if ((k-1)/Eout .gt. 100) then
!	  write(6,*) Eout, k-1, (k-1)/Eout
	  	  dc(1:nrec)=0.
	  else 
	  dc(1:nrec)=(k-1)/Eout
	  endif
 
	! 
      else
	  k=0 
	  Eout = 0.
	     do i=1,nrec
		 Eout1 =0.
           do j=1,nt(i)

           Eout1 = Eout1+((dobs(k+j)-CFR(i)*dout(k+j))/sdev(k+j))**2/2

           enddo
             dc(i)=sqrt(2*Eout1/nt(i))
			 CHI(i)=nt(i)/2*log(Eout1)
             Eout=Eout+nt(i)/2*log(Eout1)  		
	          k=sum(nt(1:i))
	      enddo 
TCHI=Eout;
		  
      endif
     END SUBROUTINE FORWARD2
!=========================================================================================
