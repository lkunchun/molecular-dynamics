      program genme
c     parameter(lattice_const=1)
      parameter(lattice_const=1.1)
      integer lattice, j,i,n

      read(*,*) lattice
      do i=0, lattice - 1
         do j=0, lattice - 1
            print*, 'atom 1 1',float(i)*lattice_const, 
     1                         float(j)*lattice_const
         enddo
      enddo
      end
