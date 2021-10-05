!===================================================
!  GLOBAL KEYWORDS
!===================================================
TemperatureList[K]                     ${temperatures}
PressureList[atm]                      ${pressures}
!
EnergyStepOverTemperature              .2
% if excess_ene_temp is not None:
ExcessEnergyOverTemperature            ${excess_ene_temp}
% endif
ModelEnergyLimit[kcal/mol]             800
!
CalculationMethod                      direct
!
WellCutoff                             10
% if well_extension is not None:
WellExtension                          ${well_extension}
% endif
ChemicalEigenvalueMax                  0.2
!
ReductionMethod                        ${reduction_method} 
% if reduction_method == 'well_reduction':
WellReductionThreshold                 ${well_reduction_thresh} 
% endif
!
AtomDistanceMin[bohr]                  1.3
RateOutput                             rate.out
% if float_type == 'quadruple':
!
FloatType                              dd
% endif
