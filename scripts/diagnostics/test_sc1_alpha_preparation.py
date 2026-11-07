#!/usr/bin/env python3
"""Private-library preparation and fallback guards; not an OpenFOAM build."""
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import prepare_sc1_conservation_solver as base
from prepare_sc1_alpha_audit import (
    DEPENDENCIES, HEADERS, MARKER, bind_private_library, main, prepare,
    sha256, source_hashes,
)


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root=Path(self.directory.name)
        self.source=self.root/'live/solver'
        (self.source/'Make').mkdir(parents=True)
        for name in ('putidaFoam.C','YEqns.H'):
            text=subprocess.check_output([
                'git','show','8e10fc70f0048b1e414d7d92701781c5e6122ef0:'
                'OFsolvers/birdmultiphaseEulerFoam/putidaFoam/'+name],cwd=base.ROOT,text=True)
            (self.source/name).write_text(text+'\n// LOCAL_BIOLOGY_SENTINEL\n')
        (self.source/'Make/files').write_text('putidaFoam.C\nEXE = $(FOAM_USER_APPBIN)/localSolver\n')
        (self.source/'Make/options').write_text(
            'EXE_INC = -DLOCAL_BIOLOGY -I../phaseSystems/lnInclude\n'
            'EXE_LIBS = \\' '\n'
            '    -L$(FOAM_USER_LIBBIN) \\' '\n'
            '    -lbirdphaseSystem -llocalBiology\n'
            'EXE_LIBS += -L/local/torch -ltorch\n')
        self.phase=self.root/'installed_sources/phaseSystems'
        shutil.copytree(base.ROOT/'OFsolvers/birdmultiphaseEulerFoam/phaseSystems',self.phase,
                        ignore=shutil.ignore_patterns('lnInclude','linux*'))
        (self.source.parent/'phaseSystems').symlink_to(self.phase,target_is_directory=True)
        for name in DEPENDENCIES:
            if name=='phaseSystems':
                continue
            target=self.root/'installed_sources'/name
            (target/'lnInclude').mkdir(parents=True)
            (target/'lnInclude/existing.H').write_text('// immutable dependency\n')
            (self.source.parent/name).symlink_to(target,target_is_directory=True)
        (self.phase/'lnInclude').mkdir()
        (self.phase/'lnInclude/phaseSystem.H').symlink_to(self.phase/'phaseSystem/phaseSystem.H')
        (self.phase/'Make/linux64Test').mkdir()
        (self.phase/'Make/linux64Test/original.o').write_bytes(b'original object')
        (self.phase/'source_alias.H').symlink_to(self.phase/'phaseSystem/phaseSystem.H')
        self.installed=self.root/'installed_lib/libbirdphaseSystem.so'
        self.installed.parent.mkdir()
        self.installed.write_bytes(b'original installed binary; never a build target')
        self.out=self.root/'private_audit'

    def base_prepare(self, source, out):
        # Exercise the real existing preparer while restricting its scratch guard
        # to this test's isolated temporary fixture.
        with patch.object(base,'scratch',lambda p:Path(p).resolve()), \
             patch.object(sys,'argv',['prepare','--source',str(source),'--out',str(out)]), \
             contextlib.redirect_stdout(io.StringIO()):
            if base.main()!=0:
                raise ValueError('Fixture conservation preparation failed')

    def prepare(self, **kwargs):
        return prepare(self.source,self.out,self.installed,
                       base_prepare=kwargs.get('base_prepare',self.base_prepare))

    def test_private_copy_targets_and_complete_provenance_preserve_originals(self):
        before=source_hashes(self.phase)
        source_before=source_hashes(self.source)
        binary_before=sha256(self.installed)
        manifest=self.prepare()
        private=self.out/'phaseSystems'
        self.assertFalse(private.is_symlink())
        self.assertFalse(any(p.is_symlink() for p in private.rglob('*')))
        self.assertFalse((private/'lnInclude').exists())
        self.assertFalse((private/'Make/linux64Test').exists())
        self.assertEqual((self.phase/'Make/linux64Test/original.o').read_bytes(),b'original object')
        self.assertEqual(source_hashes(self.phase),before)
        self.assertEqual(source_hashes(self.source),source_before)
        self.assertEqual(sha256(self.installed),binary_before)
        self.assertEqual(manifest['private_phase_source_sha256_before'],before)
        self.assertEqual(manifest['private_phase_source_sha256_after'],source_hashes(private))
        self.assertEqual(manifest['installed_library_sha256'],binary_before)
        self.assertEqual(manifest['solver_source_sha256_after_alpha'],source_hashes(self.out/'solver'))
        self.assertEqual(set(manifest['source_dependencies']),set(DEPENDENCIES))
        changed={name for name in before if before[name]!=manifest['private_phase_source_sha256_after'][name]}
        self.assertEqual(changed,{'Make/files','phaseSystem/phaseSystemSolve.C'})
        for name in HEADERS:
            self.assertTrue((private/'phaseSystem'/name).is_file())
        self.assertIn('LIB = '+str(self.out/'lib/libbirdphaseSystem'),(private/'Make/files').read_text())
        self.assertEqual((private/'Make/options').read_bytes(),(self.phase/'Make/options').read_bytes())
        self.assertIn('EXE = '+str(self.out/'bin/putidaFoamSC1AlphaAudit'),
                      (self.out/'solver/Make/files').read_text())
        updated=json.loads((self.out/'source_manifest.json').read_text())
        self.assertEqual(updated['dependencies']['phaseSystems'],str(private))
        self.assertTrue((self.out/'alpha_audit_manifest.json').is_file())
        for name in DEPENDENCIES:
            if name!='phaseSystems':
                self.assertTrue((self.out/name).is_symlink())

    def test_solver_keeps_local_biology_v2_hooks_and_mandatory_marker(self):
        self.prepare()
        main_text=(self.out/'solver/putidaFoam.C').read_text()
        original=(self.source/'putidaFoam.C').read_text()
        for marker in ('LOCAL_BIOLOGY_SENTINEL','#include "bioReact.H"',
                       '#include "phReact.H"','reaction_time += reaction_update_time.value();'):
            self.assertEqual(main_text.count(marker),original.count(marker))
        self.assertIn('sc1MassStages->beforeBlock("alphaPB")',main_text)
        self.assertIn('sc1Conservation->sample',main_text)
        self.assertIn('extern "C" int '+MARKER+'();',main_text)
        self.assertLess(main_text.index('SC1_ALPHA_LIBRARY version='),main_text.index('while (pimple.run'))
        self.assertIn('    }\n\n    const int sc1AlphaLibraryVersion',main_text)
        self.assertIn('LOCAL_BIOLOGY_SENTINEL',(self.out/'solver/YEqns.H').read_text())
        for name in ('sc1ConservationAudit.H','sc1MassStageAudit.H','sc1BudgetLedger.H','sc1MassStageLedger.H'):
            self.assertTrue((self.out/'solver'/name).is_file())

    def test_effective_link_search_order_preserves_existing_flags(self):
        self.prepare()
        make=shutil.which('make')
        if not make:
            self.skipTest('make unavailable')
        options=(self.out/'solver/Make/options').read_text()
        file=self.root/'print.Makefile'
        file.write_text('FOAM_USER_LIBBIN=/original/lib\n'+options+
                        '\n.PHONY: show\nshow:\n\t@printf "%s\\n" "$(EXE_LIBS)"\n')
        flags=subprocess.check_output([make,'--no-print-directory','-s','-f',str(file),'show'],text=True).split()
        self.assertEqual(flags[:3],['-L'+str(self.out/'lib'),'-L/original/lib','-lbirdphaseSystem'])
        self.assertIn('-llocalBiology',flags)
        self.assertIn('-L/local/torch',flags)
        self.assertIn('-lthermophysicalTransportModels',flags)

    def test_existing_and_dangling_destinations_are_never_reused(self):
        self.out.mkdir()
        (self.out/'sentinel').write_text('keep')
        with self.assertRaisesRegex(ValueError,'Destination exists'):
            self.prepare()
        self.assertEqual((self.out/'sentinel').read_text(),'keep')
        dangling=self.root/'dangling'
        dangling.symlink_to(self.root/'absent')
        with self.assertRaisesRegex(ValueError,'Destination exists'):
            prepare(self.source,dangling,self.installed,base_prepare=self.base_prepare)
        self.assertTrue(dangling.is_symlink())

    def test_partial_failure_preserves_tree_and_blocks_retry(self):
        def broken_base(source,out):
            self.base_prepare(source,out)
            (out/'solver/Make/options').write_text('EXE_INC = -DUNKNOWN\n')
        with self.assertRaisesRegex(ValueError,'EXE_LIBS'):
            self.prepare(base_prepare=broken_base)
        self.assertTrue(self.out.exists())
        self.assertFalse((self.out/'alpha_audit_manifest.json').exists())
        with self.assertRaisesRegex(ValueError,'Destination exists'):
            self.prepare()

    def test_missing_dependency_and_library_fail_before_creating_output(self):
        headers=self.source.parent/'multiphaseCompressibleMomentumTransportModels/lnInclude'
        headers.rename(headers.with_name('not_lnInclude'))
        with self.assertRaisesRegex(ValueError,'prebuilt dependency headers'):
            self.prepare()
        self.assertFalse(self.out.exists())
        with self.assertRaises(OSError):
            prepare(self.source,self.out,self.root/'missing/libbirdphaseSystem.so',base_prepare=self.base_prepare)
        self.assertFalse(self.out.exists())

    def test_already_instrumented_phase_source_is_rejected(self):
        from instrument_sc1_alpha_subcycles import instrument_phase_system_solve
        path=self.phase/'phaseSystem/phaseSystemSolve.C'
        path.write_text(instrument_phase_system_solve(path.read_text()))
        with self.assertRaises(ValueError):
            self.prepare()
        self.assertFalse(self.out.exists())

    def test_source_cycle_and_unknown_library_target_fail_closed(self):
        (self.phase/'cycle').symlink_to(self.phase,target_is_directory=True)
        with self.assertRaisesRegex(ValueError,'cycle'):
            self.prepare()
        (self.phase/'cycle').unlink()
        (self.phase/'Make/files').write_text('phaseSystem/phaseSystemSolve.C\n')
        with self.assertRaisesRegex(ValueError,'LIB target'):
            self.prepare()
        self.assertFalse(self.out.exists())

    def test_unsafe_make_paths_are_rejected_before_preparation(self):
        for name in ('unsafe space','unsafe;command','unsafe`command`','unsafe$variable'):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError,'Make path'):
                prepare(self.source,self.root/name,self.installed,base_prepare=self.base_prepare)
            self.assertFalse((self.root/name).exists())

    def test_destination_inside_any_dependency_is_rejected_before_copy(self):
        for name in DEPENDENCIES:
            dependency=(self.source.parent/name).resolve()
            out=dependency/'new_audit'
            with self.subTest(name=name),self.assertRaisesRegex(ValueError,'overlaps source dependency'):
                prepare(self.source,out,self.installed,base_prepare=self.base_prepare)
            self.assertFalse(out.exists())

    def test_destination_cannot_descend_through_installed_library_file(self):
        out=self.installed/'new_audit'
        before=sha256(self.installed)
        with self.assertRaisesRegex(ValueError,'overlaps installed library'):
            prepare(self.source,out,self.installed,base_prepare=self.base_prepare)
        self.assertEqual(sha256(self.installed),before)

    def test_cli_rejects_non_scratch_destination_without_creating_it(self):
        output=io.StringIO()
        with contextlib.redirect_stdout(output):
            status=main(['--source',str(self.source),'--out',str(self.out),
                         '--installed-library',str(self.installed)])
        self.assertEqual(status,1)
        self.assertIn('Must be below nwimer scratch',output.getvalue())
        self.assertFalse(self.out.exists())

    def test_mandatory_marker_rejects_original_library_at_link_and_runtime(self):
        compiler=shutil.which('g++')
        if not compiler:
            self.skipTest('g++ unavailable; runtime fallback guard not verified')
        self.prepare()
        solver=(self.out/'solver/putidaFoam.C').read_text()
        self.assertIn('const int sc1AlphaLibraryVersion = '+MARKER+'();',solver)
        original=self.root/'link_original'; private=self.root/'link_private'
        original.mkdir(); private.mkdir()
        old=self.root/'old.C'; new=self.root/'new.C'; main_source=self.root/'main.C'
        old.write_text('extern "C" int originalLibrary(){return 1;}\n')
        new.write_text('extern "C" int '+MARKER+'(){return 1;}\n')
        main_source.write_text('extern "C" int '+MARKER+'();\nint main(){return '+MARKER+'()==1?0:1;}\n')
        for source,directory in ((old,original),(new,private)):
            subprocess.run([compiler,'-shared','-fPIC',str(source),'-o',str(directory/'libbirdphaseSystem.so')],check=True)
        exe=self.root/'marker_test'
        wrong=subprocess.run([compiler,str(main_source),'-L'+str(original),'-lbirdphaseSystem','-o',str(exe)],
                             text=True,capture_output=True)
        self.assertNotEqual(wrong.returncode,0)
        self.assertIn(MARKER,wrong.stderr)
        subprocess.run([compiler,str(main_source),'-L'+str(private),'-lbirdphaseSystem','-o',str(exe)],check=True)
        for directory,success in ((private,True),(original,False)):
            environment=dict(os.environ,LD_LIBRARY_PATH=str(directory))
            result=subprocess.run([str(exe)],env=environment,text=True,capture_output=True)
            self.assertEqual(result.returncode==0,success)
            if not success:
                self.assertIn(MARKER,result.stderr)


if __name__=='__main__':
    unittest.main(verbosity=2)
