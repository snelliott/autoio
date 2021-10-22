Bimolecular ${bimolec_label} 
!---------------------------------------------------
% if not calc_spc1_density:
  Fragment ${spc1_label}
% else:
  Fragment ${spc1_label} density
% endif
${spc1_data}\
% if not isatom1:
      ZeroEnergy[kcal/mol]    0.0
% endif
  End  ! Frag1
!---------------------------------------------------
% if not calc_spc2_density:
  Fragment ${spc2_label}
% else:
  Fragment ${spc2_label} density
% endif
${spc2_data}\
% if not isatom2:
      ZeroEnergy[kcal/mol]    0.0
% endif
  End  ! Frag2
!---------------------------------------------------
  GroundEnergy[kcal/mol]    ${ground_ene}
End  ! Bimol\
