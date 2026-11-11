"""Separate near-bound fixtures; do not change the completed baseline recipes."""
from sc1_two_field_cases import case_files as baseline_files, foam
from sc1_two_field_report import AuditContract

VERSION = 1
LEVELS = ('coarse', 'fine')
CASES = tuple(f'bound_{level}_{reference}' for level in LEVELS
              for reference in ('no_reference','gas_reference'))
ALPHA = {'gas': '0.00001', 'liquid': '0.99999'}
FLOOR = 1e-4


def contract(name):
    if name not in CASES:
        raise ValueError('Unknown near-bound case: '+name)
    fine=name.startswith('bound_fine_')
    return AuditContract(
        dt=.0005 if fine else .001, steps=200 if fine else 100,
        volumes=(1e-8,.00099999),
        volume_initial_tol=(1e-16,1e-13), mass_initial_tol=(1e-16,1e-10),
        volume_identity_tol=(1e-16,1e-16), mass_identity_tol=(1e-18,1e-12),
        frozen_mass_tol=(1e-16,1e-11), equation_l1_tol=(1e-16,1e-14),
        budget_arithmetic_tol=(1e-18,1e-12),
        budget_absolute=(1e-15,1e-12), budget_relative=(1e-3,1e-6))


def case_files(name):
    spec=contract(name)
    files=baseline_files('source_gas_reference' if name.endswith('gas_reference')
                         else 'source_no_reference')
    for phase,alpha in ALPHA.items():
        path='0/alpha.'+phase
        if files[path].count('uniform 0.5;')!=3:
            raise ValueError('Unexpected baseline alpha field structure')
        # Initial value and both outlet/reverse-flow values stay complementary.
        files[path]=files[path].replace('uniform 0.5;',f'uniform {alpha};')
    files['system/setFieldsDict']=foam('setFieldsDict',
        'defaultFieldValues (volScalarFieldValue alpha.liquid '+ALPHA['liquid']+
        ' volScalarFieldValue alpha.gas '+ALPHA['gas']+');\nregions ();')
    if spec.dt!=.001:
        text=files['system/controlDict']
        for before,after in (('deltaT 0.001;',f'deltaT {spec.dt};'),
                             ('maxDeltaT 0.001;',f'maxDeltaT {spec.dt};'),
                             ('writeInterval 100;',f'writeInterval {spec.steps};')):
            if text.count(before)!=1:
                raise ValueError('Unexpected baseline time controls: '+before)
            text=text.replace(before,after)
        files['system/controlDict']=text
    return files
