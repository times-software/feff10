! Written by J Kas - 7/29/2014
REAL(8) FUNCTION IntDoubleLorentz(rem1,rem2,gamch,gam,a,b,omega,iinf)
  REAL(8), INTENT(IN) :: rem1, rem2, gamch, gam, a, b, omega
  REAL(8),PARAMETER :: pi = 3.1415926535897932385
  IF(iinf.GE.0) THEN
     IntDoubleLorentz = (2.d0*gam*(a*(gam**2 - gamch**2 + (rem1 - rem2)**2) + &
     &       b*(gam**2*rem1 + gamch**2*(rem1 - 2*rem2) + rem1*(rem1 - rem2)**2))* &
     &     ATAN((omega - rem1)/gamch) + &
     &    gamch*(2.d0*(a*(-gam**2 + gamch**2 + (rem1 - rem2)**2) + &
     &          b*((gamch**2 + (rem1 - rem2)**2)*rem2 + gam**2*(-2.d0*rem1 + rem2)))* &
     &        ATAN((omega - rem2)/gam) + &
     &       gam*(2.d0*a*(-rem1 + rem2) + b*(gam**2 - gamch**2 - rem1**2 + rem2**2))* &
     &        (LOG(gamch**2 + (omega - rem1)**2) - LOG(gam**2 + (omega - rem2)**2))))/ &
     &  (2.d0*gam*gamch*((gam - gamch)**2 + (rem1 - rem2)**2)* &
     &    ((gam + gamch)**2 + (rem1 - rem2)**2))
  ELSE
     !IntDoubleLorentz =         ((a*(2*(gam + gamch)*pi + rem1 - rem2))/ &
     !&     ((gam + gamch)**2 + (rem1 - rem2)**2) + & 
     !&    (a*(-rem1 + rem2))/((gam - gamch)**2 + (rem1 - rem2)**2))/(4.*gam*gamch)
     IntDoubleLorentz = a*(gam + gamch)*pi/((rem1-rem2)**2+(gam+gamch)**2)/(2.d0*gam*gamch)
  END IF
  IntDoubleLorentz = IntDoubleLorentz*gam/pi
END FUNCTION IntDoubleLorentz
