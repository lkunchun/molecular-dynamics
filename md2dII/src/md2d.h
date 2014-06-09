C     MD2D version 1.0 (01.08.02)  
C     ---use for educational purpose only
C     ---program to simulation a system of 2-d Lennard-Jones spheres
C     ---Ref: David Chandler, Lecture Notes, Stat. Mech. of Liquids.
C      
C     Contact Person: Kun-Chun Lee
C                     klee@chem.ucla.edu
C
C     Type 'man md2d' for more info
C
C---------------------------------------------------------------------

      implicit double precision(a-h, o-z)
      include 'parser_f.h'
      integer MAXCYCLE, ENE_OUTPUT, MAXOUTPUT, MAXCHARS, MINCHARS, MAXATOM,MAX_NEIGHBOR 
      integer DIM, NULL, NOTNULL, XX,YY, config_out, restart_out1, restart_out2, 
      integer thermo_out, read_line
      parameter(MAXCYCLE=2000000)
      parameter(ENE_OUTPUT=8)
      parameter(MAXOUTPUT=8000)
      parameter (MAXCHARS=80)
      parameter (MINCHARS=20)
      parameter (MAXATOM=256)
      parameter (MAX_NEIGHBOR=64, MAX_RCUT=4.5, DISPMAX=5e-2,
     1           NEIGHBOR_CUT=4.0d0)
      parameter(NEUPDATE=10)
      parameter(DIM=2)

      parameter (ZERO=0.0d0)
      parameter (NULL=-1, NOTNULL = 1)
      parameter (XX=1, YY=2)
      parameter (config_out=51)
      parameter (restart_out1=52, restart_out2=53, thermo_out=54)
      parameter (read_line = 31)

C---------------------------------------------------------------------

      dimension pos(DIM,MAXATOM), pos_t(DIM,MAXATOM),
     1          force(DIM,MAXATOM), vel(DIM,MAXATOM), 
     2          size(MAXATOM), sqrtmass(MAXATOM)
      dimension neighbor_list(MAX_NEIGHBOR,MAXATOM)
      dimension boxlength(DIM)

      character*20 config
      character*20 thermo
      character*20 restart

C----------------------------------------------------------------------

      common /env/ time, stop_time, time_step, boxlength 
      common /env_counters/ myseed, natoms, interval_output
      common /env_flags/ irestartQ
      common /atom/ pos, force, size, 
     1              neighbor_list, pos_t, vel, sqrtmass
      common /thermo/ energy, ene_kinetic, ene_pot, energy2, 
     1                pressure, pressure2, ave_ene, ave_ene2,
     2                ave_pres, ave_pres2

C    termperature_in is specified by the user, temperature  is used to 
C    measure the kinetic energy of the system

C  Michael Lewis: change ncollision to be a double to store lower freq 
      common /heatbath/ friction, temperature_in, 
     1                  temperature, sqrttmp, ncollision
	  common /scalebath/ iscale_freq, iscale_max, iscale_int
      common /dynamic/ ave_disp2, ave_vel, ave_vel2
      common /constants/ pi, twopi, sqrtpi, sqrttwo
      common /file_name/ config, thermo, restart
      common /output_frequency/ nconfig, nthermo, nrestart

C----------------------------------------------------------------------

   31 format(a80)

