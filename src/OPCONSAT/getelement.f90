FUNCTION GetElement(iz)
  INTEGER iz
  CHARACTER(2) GetElement
  CHARACTER(2) Elements(100)


  DATA Elements(1) /'H '/
  DATA Elements(2) /'He'/
  DATA Elements(3) /'Li'/
  DATA Elements(4) /'Be'/
  DATA Elements(5) /'B '/
  DATA Elements(6) /'C '/
  DATA Elements(7) /'N '/
  DATA Elements(8) /'O '/
  DATA Elements(9) /'F '/
  DATA Elements(10) /'Ne'/
  DATA Elements(11) /'Na'/
  DATA Elements(12) /'Mg'/
  DATA Elements(13) /'Al'/
  DATA Elements(14) /'Si'/
  DATA Elements(15) /'P '/
  DATA Elements(16) /'S '/
  DATA Elements(17) /'Cl'/
  DATA Elements(18) /'Ar'/
  DATA Elements(19) /'K '/
  DATA Elements(20) /'Ca'/
  DATA Elements(21) /'Sc'/
  DATA Elements(22) /'Ti'/
  DATA Elements(23) /'V '/
  DATA Elements(24) /'Cr'/
  DATA Elements(25) /'Mn'/
  DATA Elements(26) /'Fe'/
  DATA Elements(27) /'Co'/
  DATA Elements(28) /'Ni'/
  DATA Elements(29) /'Cu'/
  DATA Elements(30) /'Zn'/
  DATA Elements(31) /'Ga'/
  DATA Elements(32) /'Ge'/
  DATA Elements(33) /'As'/
  DATA Elements(34) /'Se'/
  DATA Elements(35) /'Br'/
  DATA Elements(36) /'Kr'/
  DATA Elements(37) /'Rb'/
  DATA Elements(38) /'Sr'/
  DATA Elements(39) /'Y '/
  DATA Elements(40) /'Zr'/
  DATA Elements(41) /'Nb'/
  DATA Elements(42) /'Mo'/
  DATA Elements(43) /'Tc'/
  DATA Elements(44) /'Ru'/
  DATA Elements(45) /'Rh'/
  DATA Elements(46) /'Pd'/
  DATA Elements(47) /'Ag'/
  DATA Elements(48) /'Cd'/
  DATA Elements(49) /'In'/
  DATA Elements(50) /'Sn'/
  DATA Elements(51) /'Sb'/
  DATA Elements(52) /'Te'/
  DATA Elements(53) /'I '/
  DATA Elements(54) /'Xe'/
  DATA Elements(55) /'Cs'/
  DATA Elements(56) /'Ba'/
  DATA Elements(57) /'La'/
  DATA Elements(58) /'Ce'/
  DATA Elements(59) /'Pr'/
  DATA Elements(60) /'Nd'/
  DATA Elements(61) /'Pm'/
  DATA Elements(62) /'Sm'/
  DATA Elements(63) /'Eu'/
  DATA Elements(64) /'Gd'/
  DATA Elements(65) /'Tb'/
  DATA Elements(66) /'Dy'/
  DATA Elements(67) /'Ho'/
  DATA Elements(68) /'Er'/
  DATA Elements(69) /'Tm'/
  DATA Elements(70) /'Yb'/
  DATA Elements(71) /'Lu'/
  DATA Elements(72) /'Hf'/
  DATA Elements(73) /'Ta'/
  DATA Elements(74) /'W '/
  DATA Elements(75) /'Re'/
  DATA Elements(76) /'Os'/
  DATA Elements(77) /'Ir'/
  DATA Elements(78) /'Pt'/
  DATA Elements(79) /'Au'/
  DATA Elements(80) /'Hg'/
  DATA Elements(81) /'Tl'/
  DATA Elements(82) /'Pb'/
  DATA Elements(83) /'Bi'/
  DATA Elements(84) /'Po'/
  DATA Elements(85) /'At'/
  DATA Elements(86) /'Rn'/
  DATA Elements(87) /'Fr'/
  DATA Elements(88) /'Ra'/
  DATA Elements(89) /'Ac'/
  DATA Elements(90) /'Th'/
  DATA Elements(91) /'Pa'/
  DATA Elements(92) /'U '/
  DATA Elements(93) /'Np'/
  DATA Elements(94) /'Pu'/
  DATA Elements(95) /'Am'/
  DATA Elements(96) /'Cm'/
  DATA Elements(97) /'Bk'/
  DATA Elements(98) /'Cf'/
  DATA Elements(99) /'Es'/
  DATA Elements(100) /'Fm'/

  GetElement = Elements(iz)
  RETURN
END FUNCTION GetElement
