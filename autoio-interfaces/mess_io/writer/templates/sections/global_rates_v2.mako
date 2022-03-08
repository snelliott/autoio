!===================================================
!  GLOBAL KEYWORDS
!===================================================
TemperatureList[K]                     ${temperatures}
PressureList[atm]                      ${pressures}
! 
ReferenceTemperature[K]                ${ref_temperature}
ReferencePressure[atm]                 ${ref_pressure}
!
EnergyStepOverTemperature              .2
EnergyCutoffOverTemperature            ${ene_cutoff_temp}
ExcessEnergyOverTemperature            ${excess_ene_temp}
ModelEnergyLimit[kcal/mol]             800
!
ChemicalTolerance                      ${chem_tol}
WellProjectionThreshold                ${well_projection_thresh}
WellReductionThreshold                 ${well_reduction_thresh}
!
TimePropagationLimit                   ${time_propagation_limit}
TimePropagationStep                    ${time_propagation_step}
!
WellExtension                          ${well_extension}
% if ped_spc_str is not None:
!
PEDOutput                              ped.out
PEDSpecies                             ${ped_spc_str}
% endif
% if hot_ene_str is not None:
!
HotEnergies[kcal/mol]                  ${nhot}
${hot_ene_str}
% endif
% if micro_out_params is not None:
!
MicroRateOutput                        ke.out
MicroEnerMin[kcal/mol]                 ${micro_out_params[0]}
MicroEnerMax[kcal/mol]                 ${micro_out_params[1]}
MicroEnerStep[kcal/mol]                ${micro_out_params[2]}
% endif
!
AtomDistanceMin[bohr]                  1.3
RateOutput                             rate.out
% if float_type == 'quadruple':
!
FloatType                              dd
% endif
