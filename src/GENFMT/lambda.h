!     common /lambda/  
!    4   mlam(lamtot), 	!mu for each lambda
!    5   nlam(lamtot),	!nu for each lambda
!    1   lamx, 		!max lambda in problem
!    2   laml0x, 	!max lambda for vectors involving absorbing atom
!    3   mmaxp1, nmax 	!max mu in problem + 1, max nu in problem
      common /lambda/ mlam(lamtot), nlam(lamtot), lamx, laml0x,         &
     &                mmaxp1, nmax
