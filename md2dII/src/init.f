      subroutine init(parser_flag)
      include 'md2d.h'
      integer parser_flag
      parser_flag=0
      
      
      pi = 4.0*atan(1.0)
      twopi = pi*2.0
      sqrtpi = sqrt(pi)
      sqrttwo = sqrt(2.0)

      do i = 1, DIM
         boxlength(i) = 10.0d0
      enddo

      natoms = 0
      ave_disp2=ZERO
      time = ZERO
      myseed = 12345
      stop_time = 250.0
      time_step = 0.015
      friction = 1.0
      temperature_in = 1.2
      temperature = 1.2
      sqrttmp = sqrt(temperature_in)
      ncollision = 10
      iscale_freq=0
      iscale_int=300
      iscale_max=100
      ave_vel = ZERO
      ave_vel2 = ZERO
      interval_output=4000
      ngrow=NULL
      irestartQ=NULL
      ave_ene = ZERO
      ave_ene2 = ZERO
      ave_pres = ZERO
      ave_pres2 = ZERO
      energy=ZERO
      energy2=ZERO

      config = 'testing/config.out'
      thermo = 'testing/thermo.out'
      restart = 'testing/restart.out'


      read(*,*) string_buffer

      if(string_buffer.eq.'@!') then
         print*, string_buffer
         call parser(parser_flag)
      else
         write(*,*) "the calling sequence is './runme < input > output'"
         stop
c        call md2d_data_in()
      endif

      do i = 1, natoms
         sqrtmass(i) = 1.0/sqrt(sqrtmass(i))
      enddo

      if(irestartQ.eq.NULL) then
         call velocity_scaling()
      endif

      end

C-----------------------------------------------------------
C    DIRECT MODE INPUT

      subroutine md2d_data_in()
      include 'md2d.h'

      natoms = 1
      open(unit=10, file='md2d.data', status = 'old')

 191  read(10,*,end=192) size(natoms),sqrtmass(natoms),
     1                  pos(XX,natoms),pos(YY,natoms)

c        print*, size(natoms), sqrtmass(natoms), pos(XX,natoms),
c    1           pos(YY,natoms)
         natoms = natoms + 1
         
         goto 191

 192  continue

      close(10)

      natoms = natoms - 1 

      end
