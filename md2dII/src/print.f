      subroutine print_atoms()
      include 'md2d.h'
c     print*, '#', time, ene_pot/natoms, 
c    1 (0.5*ene_kinetic/natoms),
c    2 pressure/2/temperature_in/boxlength(XX)/boxlength(YY)
      call print_ene()
C     ene_pot=ZERO
C     ene_kinetic=ZERO
   32 format('i',i4,8f10.5)
      do i = 1, natoms
          tempposx = pos(XX,i)
     1             -nint(pos(XX,i)/boxlength(XX))*boxlength(XX)
          tempposy = pos(YY,i)
     1             -nint(pos(YY,i)/boxlength(YY))*boxlength(YY)

          temp_m = 1.0/sqrtmass(i)
          temp_m = temp_m*temp_m
          write(*,32) i, size(i), temp_m, tempposx, tempposy, 
     1            vel(XX,i), 
     2            vel(YY,i), 
     3            force(XX,i), 
     4            force(YY,i) 
      enddo

      end

      subroutine print_ene()
      include 'md2d.h'
   32 format('#',8f12.7)
      write(*,32) (ene_pot+0.50d0*ene_kinetic)/natoms, ene_pot/natoms, 
     1 (0.50d0*ene_kinetic/natoms),
     2 pressure, ave_ene/interval_output, ave_ene2/interval_output,
     3 ave_pres/interval_output, ave_pres2/interval_output

      ave_ene = ZERO
      ave_ene2 = ZERO
      ave_pres = ZERO
      ave_pres2 = ZERO
      energy=ZERO
      energy2=ZERO
      end

      subroutine lammps_out()

      include 'md2d.h'
      data iopenQ /NULL/
      save

      if(iopenQ.eq.NULL) then
         open(unit=config_out,file=config,status='unknown')
         iopenQ=NOTNULL
      endif

      print*, 'start'
c     do i = 1, natoms
         write(config_out,*) 'test success'
         stop
c     enddo
      end
