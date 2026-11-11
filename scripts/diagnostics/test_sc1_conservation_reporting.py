#!/usr/bin/env python3
"""Reporting regressions: native endpoints must never become synchronized data."""
import contextlib
import copy
import io
from pathlib import Path
import tempfile
import unittest

from summarize_sc1_conservation import (
    PHASES, compare, load_log, main, native_endpoint_difference,
)


def audit_text(times, start=31.2):
    """Conservative constant-rate fixture with an independently known integral."""
    lines=['SC1_BALANCE_METHOD version=1 exactDiscreteBudget=0']
    for phase in PHASES:
        lines.append(f'SC1_BALANCE_START phase={phase} t_s={start:.17g}'
                     ' mass_kg=1 oxygen_kg=0.23 volume_m3=1')
    for time in times:
        duration=time-start
        for phase in PHASES:
            data=dict(t_s=time,start_s=start,steps=round(duration*100000),
                      mass_kg=1+0.007*duration,oxygen_kg=0.23+2e-6*duration,
                      massResidual_kg=0,oxygenResidual_kg=0,
                      massExternalIntegral_kg=0.01*duration,
                      massTransferIntegral_kg=-0.004*duration,
                      massInIntegral_kg=0.02*duration,
                      massOutIntegral_kg=0.019*duration,
                      oxygenTransferIntegral_kg=2e-6*duration,
                      normalization_kg=0,equationResidualIntegral_kg=0)
            lines.append('SC1_BALANCE phase='+phase+' '+
                         ' '.join(f'{k}={v:.17g}' for k,v in data.items()))
    return '\n'.join(lines+['End',''])


class ReportingTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root=Path(self.directory.name)

    def write(self, name, text):
        path=self.root/name
        path.write_text(text)
        return path

    def run_main(self, *args):
        output=io.StringIO()
        with contextlib.redirect_stdout(output):
            status=main([str(arg) for arg in args])
        return status,output.getvalue()

    def test_default_preserves_synchronized_comparison(self):
        path=self.write('control',audit_text([31.3,31.4]))
        loaded=load_log(path)
        self.assertEqual(compare(loaded,loaded),[(31.3,0,0),(31.4,0,0)])
        status,output=self.run_main(path,'--normal',path)
        self.assertEqual(status,0)
        self.assertIn('Time(s)  normal-control',output)

    def test_late_microsecond_offsets_still_fail_default_after_both_budgets(self):
        c=self.write('control',audit_text([31.3,31.400018618094]))
        n=self.write('normal',audit_text([31.3,31.3999916830166]))
        status,output=self.run_main(c,'--normal',n)
        self.assertEqual(status,1)
        self.assertEqual(output.count('phase             mass residual'),2)
        self.assertIn('Output times differ; do not silently interpolate',output)
        self.assertNotIn('Time(s)  normal-control',output)
        self.assertLess(output.index(str(n)+': cumulative'),output.index('Audit interpretation stopped'))

    def test_native_endpoints_keep_independent_times_durations_and_rates(self):
        ctime=31.6999975003429; ntime=31.6999846530835
        c=self.write('control',audit_text([31.3,ctime]))
        n=self.write('normal',audit_text([31.3,ntime]))
        status,output=self.run_main(c,'--normal',n,'--native-times')
        self.assertEqual(status,0)
        self.assertIn('NOT a synchronized comparison',output)
        self.assertNotIn('Time(s)  normal-control',output)
        for name,time in (('control',ctime),('normal',ntime)):
            self.assertIn(f'{name}_endpoint_s={time:.17g}',output)
            self.assertIn(f'{name}_duration_s={time-31.2:.17g}',output)
        self.assertIn(f'endpoint_offset_s={ntime-ctime:.17g}',output)
        self.assertEqual(output.count('mean rates use only this run'),2)
        difference=native_endpoint_difference(load_log(c),load_log(n))
        self.assertAlmostEqual(difference['bubble_inventory_difference_mg'],
                               .007*(ntime-ctime)*1e6,places=8)
        # Each run has the same known mean O2 transfer, despite different duration.
        for loaded in (load_log(c),load_log(n)):
            r=loaded[1]['liquid'][-1]
            self.assertAlmostEqual(r['oxygenTransferIntegral_kg']/(r['t_s']-r['start_s']),2e-6)

    def test_native_mode_does_not_pair_unequal_report_counts_or_durations(self):
        c=self.write('control',audit_text([31.3,31.4]))
        n=self.write('normal',audit_text([31.35,31.5,31.8]))
        with self.assertRaisesRegex(ValueError,'Different number'):
            compare(load_log(c),load_log(n))
        status,output=self.run_main(c,'--normal',n,'--native-times')
        self.assertEqual(status,0)
        self.assertIn('reports=2',output)
        self.assertIn('reports=3',output)
        self.assertIn('no saved rows paired, no interpolation',output)

    def test_all_initial_phases_and_inventories_are_checked_in_both_modes(self):
        loaded=load_log(self.write('control',audit_text([31.3])))
        for phase in PHASES:
            for key in ('t_s','mass_kg','oxygen_kg','volume_m3'):
                with self.subTest(phase=phase,key=key):
                    changed=copy.deepcopy(loaded)
                    changed[0][phase][key]+=0.001
                    for comparison in (compare,native_endpoint_difference):
                        with self.assertRaisesRegex(ValueError,'Initial'):
                            comparison(loaded,changed)

    def test_even_subtolerance_start_clock_changes_are_rejected(self):
        loaded=load_log(self.write('control',audit_text([31.3])))
        changed=copy.deepcopy(loaded)
        for row in changed[0].values():
            row['t_s']+=1e-12
        for comparison in (compare,native_endpoint_difference):
            with self.assertRaisesRegex(ValueError,'Initial t_s'):
                comparison(loaded,changed)

    def test_initial_mismatch_does_not_suppress_valid_native_budgets(self):
        c=self.write('control',audit_text([31.3]))
        n=self.write('normal',audit_text([31.3001],start=31.2001))
        status,output=self.run_main(c,'--normal',n,'--native-times')
        self.assertEqual(status,1)
        self.assertEqual(output.count('phase             mass residual'),2)
        self.assertIn('not matched conditions',output)
        self.assertNotIn('Native endpoint difference',output)

    def test_missing_initial_volume_on_only_one_side_is_not_matched(self):
        loaded=load_log(self.write('control',audit_text([31.3])))
        changed=copy.deepcopy(loaded)
        del changed[0]['headspaceGas']['volume_m3']
        with self.assertRaisesRegex(ValueError,'volume_m3 missing'):
            native_endpoint_difference(loaded,changed)

    def test_each_valid_run_is_reported_when_its_peer_is_malformed(self):
        valid=self.write('valid',audit_text([31.3]))
        malformed=self.write('malformed',audit_text([31.3]).replace('\nEnd\n','\n'))
        for control,normal in ((valid,malformed),(malformed,valid)):
            with self.subTest(control=control.name):
                status,output=self.run_main(control,'--normal',normal,'--native-times')
                self.assertEqual(status,1)
                self.assertEqual(output.count('phase             mass residual'),1)
                self.assertIn(str(valid)+': cumulative',output)
                self.assertIn('no solver End',output)
                self.assertNotIn('Native endpoint difference',output)

    def test_malformed_records_are_rejected_instead_of_partially_summarized(self):
        valid=audit_text([31.3])
        variants={
            'nonnumeric':valid.replace('massResidual_kg=0','massResidual_kg=bad',1),
            'nonfinite':valid.replace('oxygenResidual_kg=0','oxygenResidual_kg=inf',1),
            'missingphase':valid.replace('SC1_BALANCE phase=liquid ','SC1_BALANCE ',1),
            'duplicatedkey':valid.replace('massResidual_kg=0','massResidual_kg=0 massResidual_kg=0',1),
            'fractionalsteps':valid.replace('steps=10000','steps=10000.5',1),
            'phasesteps':valid.replace('steps=10000','steps=10001',1),
            'initialclock':valid.replace('phase=liquid t_s=31.199999999999999',
                                         'phase=liquid t_s=31.2001',1),
            'initialnumber':valid.replace('oxygen_kg=0.23','oxygen_kg=bad',1),
            'missinginitial':valid.replace('SC1_BALANCE_START phase=gas ', 'OTHER phase=gas ',1),
            'missingreport':valid.replace('SC1_BALANCE phase=gas ', 'OTHER phase=gas ',1),
            'session':valid+'\n'+valid,
            'arithmetic':valid.replace('massResidual_kg=0','massResidual_kg=0.1',1),
        }
        for name,text in variants.items():
            with self.subTest(name=name):
                self.assertNotEqual(text,valid)
                with self.assertRaises(ValueError):
                    load_log(self.write(name,text))

    def test_single_native_run_always_labels_its_own_interval(self):
        path=self.write('control',audit_text([31.30002]))
        status,output=self.run_main(path,'--native-times')
        self.assertEqual(status,0)
        self.assertIn('Own interval: start_s=',output)
        self.assertIn('endpoint_s=',output)
        self.assertIn('duration_s=',output)
        self.assertNotIn('Native endpoint difference',output)


if __name__=='__main__':
    unittest.main(verbosity=2)
