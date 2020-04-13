     integer function ibravais(spacegroup,lattice)
	 ! Set the Bravais lattice
     implicit none
	 integer,intent(in) :: spacegroup
	 character*1,intent(in) :: lattice

!(1)  "triclinic   primitive      -1     C_i " 
!(2)  "monoclinic  primitive      2/m    C_2h"
!(3)  "monoclinic  base centered  2/m    C_2h"
!(4)  "orthorombic primitive      mmm    D_2h"
!(5)  "orthorombic base-centered  mmm    D_2h"
!(6)  "orthorombic body-centered  mmm    D_2h"
!(7)  "orthorombic face-centered  mmm    D_2h"
!(8)  "tetragonal  primitive      4/mmm  D_4h"
!(9)  "tetragonal  body-centered  4/mmm  D_4h"
!(10) "trigonal    primitive      -3m    D_3d"
!(11) "hexagonal   primitive      6/mmm  D_6h"
!(12) "cubic       primitive      m3m    O_h "
!(13) "cubic       face-centered  m3m    O_h "
!(14) "cubic       body-centered  m3m    O_h "

          ibravais=0
		  write(*,*) 'sgroup,lattice ',spacegroup,lattice
		  
                if (spacegroup .lt. 1  .or. spacegroup .gt. 230) then
				    call wlog('Invalid spacegroup in ibravais')
					ibravais = -1
					return
				endif

                if (spacegroup .le. 2) then 
                     ibravais = 1
                else
                    if (spacegroup .le. 15) then 
                        if (lattice.eq.'P') then 
                             ibravais = 2
                        else
                             ibravais = 3
                        endif
                    else
                        if (spacegroup .le. 74) then
						    ibravais=5
						    if(lattice.eq.'P') ibravais=4
							if(lattice.eq.'I') ibravais=6
							if(lattice.eq.'F') ibravais=7
						else
                            if (spacegroup .le. 142) then
                                if ( lattice.eq.'P' ) then 
                                     ibravais = 8
                                else
                                     ibravais = 9
                                endif
                            else
                                if (spacegroup .le. 167) then
                                     ibravais = 10
                                else
                                    if (spacegroup .le. 194) then
                                         ibravais = 11
                                    else
						                if(lattice.eq.'P') ibravais=12
							            if(lattice.eq.'I') ibravais=14
							            if(lattice.eq.'F') ibravais=13
                                    endif
                                endif
                            endif
                        endif
                    endif
                endif

           if (ibravais.lt.1 .or. ibravais.gt.14) stop 'Invalid result in ibravais'
		   
			
	return
	end 
			   
