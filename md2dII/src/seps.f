      program separations

      include 'md2d.h'
      integer parser_flag
      parser_flag=1

C     minimal initialization
      call init0()

      read(*,*) string_buffer

      if(string_buffer.eq.'@!') then
c        print*, string_buffer
         do i=0, MAXOUTPUT

C     READ IN ONE FRAME
            call parser(parser_flag)

C     DO ALL THE ANALYSIS IN CALCULATE
            call calculate()


C
C     IF THERE ARE NO MORE FRAMES THEN STOP
C
            if(parser_flag.le.0) then
               goto 191
            endif
         enddo
 191     continue
      else
         write(*,*) "the calling sequence is './runme < input > output'"
         stop
      endif
      

      end

C---------------------------------------------------------------
      subroutine calculate()
      include 'md2d.h'
      dimension boxi(DIM)
      data l_counter /0/
      save
   32 format('z', f5.2,f10.5)

      l_counter = l_counter + 1

      print*, '#',  l_counter

      do k=1, DIM
         boxi(k) = 1/boxlength(k)
      enddo

      do i=1, natoms-1
         do j=i+1, natoms
            r2=ZERO
            do k=1, DIM
               dist = pos(k,i)-pos(k,j)
               dist = dist - nint(dist*boxi(k))*boxlength(k)
               r2 = r2 + dist*dist
            enddo
            write(*,32) abs(size(i)-size(j)), sqrt(r2)
c calculate pressure
c-------------------------------------------------
         enddo
      enddo
      end

      subroutine init0()
      include 'md2d.h'
      pi = 4.0*atan(1.0)
      twopi = pi*2.0
      sqrtpi = sqrt(pi)
      sqrttwo = sqrt(2.0)

      do i = 1, DIM
         boxlength(i) = 10.0d0
      enddo

C---------------------------------------------------------------
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
      ave_vel = ZERO
      ave_vel2 = ZERO
      interval_output=4000
      ngrow=NULL
      irestartQ=NULL

      config = 'testing/config.out'
      thermo = 'testing/thermo.out'
      restart = 'testing/restart.out'

      end
