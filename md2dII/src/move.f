      subroutine move()

c     procedure to advance by one time step 
c     declaration of variables
c     kinetic energy calculation is moved to the 
c     andersen thermostat

      include 'md2d.h'
      data iflg /0/
      save

c     initialize constants

      dt = time_step
      dt2=time_step*time_step
      dt2i=1.0/(2.0*time_step)
      dti=1.0/time_step

c     verlet with position and velocity
c     write(*,*) 'verlet with positions and velocity'
c     stop

      if(iflg.eq.0.or.iscale_freq.eq.0) then
c     if(iflg.eq.1) stop

      do i=1,natoms
         temp_mi = sqrtmass(i)*sqrtmass(i)
         do j=1,DIM


            temp=pos(j,i)+vel(j,i)*dt + 0.50*force(j,i)*dt2*temp_mi

c           vel(j,i) = vel(j,i)+force(j,i)*dt*temp_mi

c           andersen thermostat moving at half step, the previous half step is
c           done in thermostat.f
            vel(j,i) = vel(j,i)+0.50*force(j,i)*dt*temp_mi

c           for scale velocity type


C           temp = 2.0*pos(j,i) - pos_t(j,i) + force(j,i)*dt2*temp_mi
c           vel(j,i) = (temp-pos(j,i))*dti
            pos_t(j,i)=pos(j,i)
            pos(j,i)=temp
            temp = (pos(j,i)-pos_t(j,i))
            ave_disp2 = ave_disp2+abs(temp)
         enddo
      enddo

C        if andersen thermostat, iflg = 0 or iscale_freq=0 always
         iflg = 1
      

      else
c     stop
      do i=1,natoms
         temp_mi = sqrtmass(i)*sqrtmass(i)
         do j=1,DIM
C           temp=pos(j,i)+vel(j,i)*dt+0.50*force(j,i)*dt2*temp_mi
C           vel(j,i) = vel(j,i) + 
C    1               0.50*force(j,i)*time_step*temp_mi

            temp = 2.0*pos(j,i) - pos_t(j,i) + force(j,i)*dt2*temp_mi
            vel(j,i) = (temp-pos_t(j,i))*dt2i
            pos_t(j,i)=pos(j,i)
            pos(j,i)=temp
            temp = (pos(j,i)-pos_t(j,i))
            ave_disp2 = ave_disp2+abs(temp)
         enddo
      enddo
      endif

      return

      end
