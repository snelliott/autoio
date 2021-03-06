% if rotor_id != '':
Rotor  Hindered   # ${rotor_id}
% else:
Rotor  Hindered
% endif
% if geom:
  Geometry[angstrom]   ${natom}
${geom}
% endif
  Group                ${group}
  Axis                 ${axis}
  Symmetry             ${symmetry}
  Potential[kcal/mol]  ${npotential}
${potential} 
% if use_quantum_weight:
  UseQuantumWeight
% endif
End
