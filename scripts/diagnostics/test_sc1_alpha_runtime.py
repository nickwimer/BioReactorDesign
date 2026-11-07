#!/usr/bin/env python3
"""Read-only runtime provenance checks; mocked ELF-tool output, not an OF build."""
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import verify_sc1_alpha_runtime as runtime
from prepare_sc1_alpha_audit import MARKER, sha256, source_hashes


class RuntimeTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.base=Path(self.directory.name)
        self.root=self.base/'audit'
        self.exe=self.root/'bin/putidaFoamSC1AlphaAudit'
        self.private=self.root/'lib/libbirdphaseSystem.so'
        self.original=self.base/'installed/libbirdphaseSystem.so'
        self.loaded=self.base/'film.so'; self.extra=self.base/'limiter.so'
        for path in (self.exe,self.private,self.original,self.loaded,self.extra):
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(b'fixture '+str(path).encode())
        self.exe.chmod(0o755)
        self.phase=self.base/'original_phase'
        for root in (self.phase,self.root/'phaseSystems'):
            root.mkdir()
            (root/'source.C').write_text('// original\n')
        (self.root/'phaseSystems/source.C').write_text('// instrumented\n')
        (self.root/'solver').mkdir()
        (self.root/'solver/main.C').write_text('// prepared solver\n')
        (self.root/'solver/sc1ConservationAudit.H').write_text('// prepared diagnostic header\n')
        self.manifest={
            'alpha_audit_version':1,'executable':str(self.exe),'private_library':str(self.private),
            'installed_library':str(self.original),'installed_library_sha256':sha256(self.original),
            'phase_source':str(self.phase),'private_phase_source_sha256_before':source_hashes(self.phase),
            'private_phase_source_sha256_after':source_hashes(self.root/'phaseSystems'),
            'solver_sha256_after_alpha':{'main.C':sha256(self.root/'solver/main.C')},
            'solver_source_sha256_after_alpha':source_hashes(self.root/'solver'),
        }
        self.save_manifest()
        self.calls=[]

    def save_manifest(self):
        (self.root/'alpha_audit_manifest.json').write_text(json.dumps(self.manifest))

    def runner(self, command, environment):
        self.calls.append((command,dict(environment)))
        if command[0]=='nm':
            stdout=('00001234 T ' if command[-1]==str(self.private) else '         U ')+MARKER+'\n'
        elif command[-1]==str(self.extra):
            stdout='libfiniteVolume.so => /existing/libfiniteVolume.so (0x1234)\n'
        else:
            stdout='libbirdphaseSystem.so => '+str(self.private)+' (0x1234)\n'
        return {'command':command,'returncode':0,'stdout':stdout,'stderr':''}

    def test_success_checks_all_artifacts_with_process_local_environment(self):
        original_env=os.environ.get('LD_LIBRARY_PATH')
        before={str(p):sha256(p) for p in self.base.rglob('*') if p.is_file()}
        report=runtime.verify(self.root,[self.loaded],[self.extra],runner=self.runner)
        self.assertTrue(report['ok'],report['errors'])
        checked=[r['path'] for r in report['checks'] if r['artifact']=='runtime dependency']
        self.assertEqual(checked,[str(self.exe),str(self.loaded),str(self.extra)])
        self.assertEqual(os.environ.get('LD_LIBRARY_PATH'),original_env)
        for command,environment in self.calls:
            self.assertTrue(environment['LD_LIBRARY_PATH'].startswith(str(self.root/'lib')))
            self.assertNotIn('-r',command)
        self.assertEqual(before,{str(p):sha256(p) for p in self.base.rglob('*') if p.is_file()})
        self.assertEqual(len(report['commands']),5)

    def test_original_and_private_source_or_library_changes_fail(self):
        for path in (self.original,self.phase/'source.C',self.root/'phaseSystems/source.C',
                     self.root/'solver/main.C',self.root/'solver/sc1ConservationAudit.H'):
            with self.subTest(path=path):
                content=path.read_bytes()
                path.write_bytes(content+b'changed')
                self.assertFalse(runtime.verify(self.root,runner=self.runner)['ok'])
                path.write_bytes(content)

    def test_private_library_symlink_fallback_is_rejected(self):
        self.private.unlink()
        self.private.symlink_to(self.original)
        report=runtime.verify(self.root,runner=self.runner)
        self.assertFalse(report['ok'])
        self.assertIn('redirected private library',report['errors'][0])
        self.assertFalse(self.calls)

    def test_missing_or_incorrect_marker_is_rejected(self):
        for symbol in ('','         U '+MARKER+'\n','00001234 W '+MARKER+'\n'):
            def wrong(command,environment):
                result=self.runner(command,environment)
                if command[0]=='nm' and command[-1]==str(self.private):
                    result['stdout']=symbol
                return result
            with self.subTest(symbol=symbol):
                report=runtime.verify(self.root,runner=wrong)
                self.assertFalse(report['ok'])
                self.assertIn('marker',report['errors'][0])

    def test_loaded_library_fallback_duplicate_and_missing_dependencies_fail(self):
        cases=(
            f'libbirdphaseSystem.so => {self.original} (0x1234)\n',
            f'libbirdphaseSystem.so => {self.private} (0x1234)\nlibother.so => not found\n',
            f'libbirdphaseSystem.so => {self.private} (0x1234)\nlibbirdphaseSystem.so => {self.original} (0x2345)\n',
            f'{self.original} (0x1234)\n',
            'libfiniteVolume.so => /existing/libfiniteVolume.so (0x1234)\n',
        )
        for output in cases:
            def wrong(command,environment):
                result=self.runner(command,environment)
                if command==['ldd',str(self.loaded)]:
                    result['stdout']=output
                return result
            with self.subTest(output=output):
                report=runtime.verify(self.root,[self.loaded],runner=wrong)
                self.assertFalse(report['ok'])
                self.assertIn(str(self.loaded),report['errors'][0])
                self.assertEqual(report['commands'][-1]['stdout'],output)

    def test_extra_library_cannot_hide_phase_fallback(self):
        def wrong(command,environment):
            result=self.runner(command,environment)
            if command==['ldd',str(self.extra)]:
                result['stdout']=f'libbirdphaseSystem.so => {self.original} (0x1234)\n'
            return result
        report=runtime.verify(self.root,[],[self.extra],runner=wrong)
        self.assertFalse(report['ok'])
        self.assertIn('fallback',report['errors'][0])

    def test_failed_tool_status_and_malformed_manifest_fail_closed(self):
        def failed(command,environment):
            result=self.runner(command,environment)
            result['returncode']=1
            return result
        self.assertFalse(runtime.verify(self.root,runner=failed)['ok'])
        (self.root/'alpha_audit_manifest.json').write_text('{bad json')
        self.assertFalse(runtime.verify(self.root,runner=self.runner)['ok'])

    def test_cli_default_is_concise_and_verbose_json_preserves_raw_outputs(self):
        report=runtime.verify(self.root,[self.loaded],[self.extra],runner=self.runner)
        args=['--root',str(self.root),'--loaded-library',str(self.loaded),'--extra-library',str(self.extra)]
        for verbose in (False,True):
            output=io.StringIO()
            with patch.object(runtime,'verify',return_value=report) as mocked,contextlib.redirect_stdout(output):
                status=runtime.main(args+(['--verbose-json'] if verbose else []))
            self.assertEqual(status,0)
            mocked.assert_called_once_with(self.root,[self.loaded],[self.extra])
            if verbose:
                self.assertEqual(json.loads(output.getvalue()),report)
            else:
                self.assertIn('PASS (no solver launched)',output.getvalue())
                self.assertNotIn('00001234 T',output.getvalue())


if __name__=='__main__':
    unittest.main(verbosity=2)
