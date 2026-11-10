#!/usr/bin/env python3
"""Prepare/verify a NEW private source-feasibility diagnostic; never launch CFD.

Existing phase equations and biology are preserved. The observer is opt-in via
alpha/sc1AlphaSourceProbe and requires the existing sc1AlphaAudit switch.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess

import prepare_sc1_alpha_audit as alpha
from instrument_sc1_alpha_subcycles import _once
from prepare_sc1_conservation_solver import scratch
from verify_sc1_alpha_runtime import command_result, verify

HERE=Path(__file__).resolve().parent
MARKER='sc1AlphaSourceProbeLibraryVersion'
HEADER='sc1AlphaSourceProbe.H'


def instrument_header(text):
    if 'sc1SourceProbe' in text:
        raise ValueError('Source probe already/partially present')
    text=_once(text,'#include "sc1AlphaSubcycleLedger.H"',
               '#include "sc1AlphaSubcycleLedger.H"\n#include "'+HEADER+'"')
    text=_once(text,'    std::map<label,Table> byCallInStep;',
               '    std::map<label,Table> byCallInStep;\n'
               '    std::map<label,sc1SourceProbe::Table> sourceProbeByCall;')
    text=_once(text,'    label callInStep_;',
               '    label callInStep_;\n    bool sourceProbe_ = false;\n'
               '    sc1SourceProbe::Table sourceProbeTerms_{};')
    text=_once(text,'        callInStep_=++totals_.callInStep;',
               '        callInStep_=++totals_.callInStep;\n'
               '        sourceProbe_=mesh_.solverDict("alpha").lookupOrDefault<Switch>\n'
               '            ("sc1AlphaSourceProbe",false);')
    anchor='            for (label j=0;j<nTerms;++j) add(i,j,c,cell[j]);'
    text=_once(text,anchor,anchor+'''
            if (sourceProbe_)
            {
                const auto probe=sc1SourceProbe::cell
                (
                    old_[i][c],phase[c],subDt,div[c],su[c],sp[c],
                    sourceAlpha_[i][c],originalSu_[i][c],dgdt_[i][c],
                    requestedDilatation_[i][c],mesh_.V()[c],entryRho_[i][c],
                    min(scalar(1),phase.alphaMax())
                );
                sc1SourceProbe::add
                (
                    sourceProbeTerms_[i][sc1SourceProbe::bin(sourceAlpha_[i][c])],
                    probe
                );
            }''')
    anchor='        if (!mesh_.time().writeTime() && mesh_.time().running()) return;'
    text=_once(text,anchor,'''        if (sourceProbe_)
        {
            auto& cumulative=totals_.sourceProbeByCall[callInStep_];
            for (unsigned i=0;i<2;++i)
                for (unsigned b=0;b<4;++b)
                    sc1SourceProbe::add(cumulative[i][b],sourceProbeTerms_[i][b]);
        }
'''+anchor)
    anchor='        Info().precision(oldPrecision);'
    text=_once(text,anchor,'''        if (sourceProbe_)
        {
            // Persistent cumulative tables remain rank-local. Reduce a COPY
            // only at reporting times, never re-add an already global sum.
            const auto& cumulative=totals_.sourceProbeByCall[callInStep_];
            const unsigned width=sc1SourceProbe::maxDtSp;
            scalarField sums(8*width,0), maxima(8,0);
            for (unsigned i=0;i<2;++i)
                for (unsigned b=0;b<4;++b)
                {
                    for (unsigned j=0;j<width;++j)
                        sums[(4*i+b)*width+j]=cumulative[i][b][j];
                    maxima[4*i+b]=cumulative[i][b][sc1SourceProbe::maxDtSp];
                }
            Pstream::listCombineGather(sums,plusEqOp<scalar>());
            Pstream::listCombineScatter(sums);
            Pstream::listCombineGather(maxima,maxEqOp<scalar>());
            Pstream::listCombineScatter(maxima);
            Info<< "SC1_SOURCE_PROBE_METHOD version=1 observerOnly=1"
                << " rawTarget=oldPlusDtTimesBaseSuPlusDgdtMinusLimitedDiv"
                << " pressureTarget=oldPlusDtTimesStoredDivUMinusLimitedDiv"
                << " boundsClassificationTolerance=1e-12 fixedFluxCounterfactual=1"
                << " coupledPhaseFeasibility=0 uniqueCellCounts=0"
                << " scope=cumulativeSameOuterIndex rhoWeight=alphaSolveEntry" << endl;
            const char* bins[4]={"lower","upper","interior","outside"};
            for (unsigned i=0;i<2;++i)
                for (unsigned b=0;b<4;++b)
                {
                    Info<< "SC1_SOURCE_PROBE phase=" << names_[i]
                        << " t_s=" << time_ << " callInStep=" << callInStep_
                        << " bin=" << bins[b];
                    for (unsigned j=0;j<width;++j)
                        Info<< ' ' << sc1SourceProbe::label(j) << '='
                            << sums[(4*i+b)*width+j];
                    Info<< " maxDtAbsSp=" << maxima[4*i+b] << endl;
                }
        }
'''+anchor)
    return text


def prepare(source,out,installed,base_prepare=alpha.run_base_preparer):
    # Validate transformation against repository v1 before creating anything.
    instrument_header((alpha.PHASE_HEADERS/'sc1AlphaSubcycleAudit.H').read_text())
    helper=(HERE/'include'/HEADER).read_bytes()
    manifest=alpha.prepare(source,out,installed,base_prepare=base_prepare)
    out=out.resolve()
    private=out/'phaseSystems/phaseSystem'
    audit=private/'sc1AlphaSubcycleAudit.H'
    audit.write_text(instrument_header(audit.read_text()))
    (private/HEADER).write_bytes(helper)
    solve=private/'phaseSystemSolve.C'
    solve.write_text(_once(solve.read_text(),
        'extern "C" int sc1AlphaAuditLibraryVersion() { return 1; }',
        'extern "C" int sc1AlphaAuditLibraryVersion() { return 1; }\n'
        'extern "C" int '+MARKER+'() { return 1; }'))
    main=out/'solver/putidaFoam.C'
    text=_once(main.read_text(),'extern "C" int sc1AlphaAuditLibraryVersion();',
               'extern "C" int sc1AlphaAuditLibraryVersion();\n'
               'extern "C" int '+MARKER+'();')
    anchor='    Info<< "SC1_ALPHA_LIBRARY version=" << sc1AlphaLibraryVersion << endl;'
    text=_once(text,anchor,anchor+'''
    if (sc1AlphaSourceProbeLibraryVersion()!=1)
        FatalErrorInFunction << "Wrong SC1 source-probe library" << exit(FatalError);
    Info<< "SC1_SOURCE_PROBE_LIBRARY version=1" << endl;''')
    main.write_text(text)
    manifest['private_phase_source_sha256_after']=alpha.source_hashes(out/'phaseSystems')
    manifest['solver_source_sha256_after_alpha']=alpha.source_hashes(out/'solver')
    manifest['solver_sha256_after_alpha']['putidaFoam.C']=alpha.sha256(main)
    manifest['source_probe_version']=1
    (out/'alpha_audit_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    completion={'version':1,'observer_only':True,'alpha_manifest_sha256':
                alpha.sha256(out/'alpha_audit_manifest.json'),
                'preparer_sha256':alpha.sha256(Path(__file__)),
                'header_sha256':alpha.sha256(private/HEADER)}
    (out/'source_probe_manifest.json').write_text(json.dumps(completion,indent=2)+'\n')
    return completion


def verify_probe(root,loaded,extra,runner=command_result):
    root=root.resolve(strict=True)
    manifest=json.loads((root/'source_probe_manifest.json').read_text())
    if (manifest['version']!=1 or manifest['observer_only'] is not True
            or manifest['alpha_manifest_sha256']!=alpha.sha256(root/'alpha_audit_manifest.json')
            or manifest['header_sha256']!=alpha.sha256(root/'phaseSystems/phaseSystem'/HEADER)):
        raise ValueError('Incomplete or changed source-probe preparation')
    report=verify(root,loaded,extra,runner=runner)
    if not report['ok']:
        raise ValueError('; '.join(report['errors']))
    env=os.environ.copy()
    env['LD_LIBRARY_PATH']=report['LD_LIBRARY_PATH']
    for binary,expected in ((root/'lib/libbirdphaseSystem.so','T'),
                            (root/'bin/putidaFoamSC1AlphaAudit','U')):
        result=runner(['nm','-D',str(binary)],env)
        types=re.findall(r'^\s*(?:[0-9a-fA-F]+\s+)?([A-Za-z?])\s+'
                         +MARKER+r'(?:@\S+)?\s*$',result['stdout'],re.M)
        if result['returncode'] or types!=[expected]:
            raise ValueError(f'Missing mandatory source-probe marker in {binary}')
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    prep=sub.add_parser('prepare')
    prep.add_argument('--source',type=Path,default=Path('/scratch/nwimer/sc1_bubble_audit_20261029/solver'))
    prep.add_argument('--out',type=Path,required=True)
    check=sub.add_parser('verify')
    check.add_argument('--root',type=Path,required=True)
    args=parser.parse_args(argv)
    try:
        libbin=os.environ.get('FOAM_USER_LIBBIN')
        if not libbin:
            raise ValueError('Source OpenFOAM and load libfabric first')
        if args.command=='prepare':
            if args.out.is_symlink():
                raise ValueError(f'Destination is a symlink: {args.out}')
            prepare(scratch(args.source),scratch(args.out),Path(libbin)/'libbirdphaseSystem.so')
            print(f'Prepared NEW observer build: {args.out}; no CFD or compilation launched.')
        else:
            report=verify_probe(scratch(args.root),[
                Path('/scratch/nwimer/sc1_headspace_o2_trial_20261103/lib/libsc1HeadspaceFilmModels.so'),
                Path('/scratch/nwimer/sc1_headspace_entrainment_20261104/lib/libsc1HeadspaceEntrainmentModels.so')
            ],[Path(libbin)/'libalphaLimitVelocity.so'])
            for check in report['checks']:
                print(check['artifact']+': '+check['path']+' '+
                      (check.get('sha256') or 'unchanged files='+str(check.get('files'))))
                if check.get('phase_library'):
                    print('  libbirdphaseSystem.so => '+check['phase_library'])
            print('Source-probe runtime verification PASS; no CFD launched.')
        return 0
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print(f'SC1 source probe stopped: {error}\nPreserve all existing/partial build trees.')
        return 1


if __name__=='__main__':
    raise SystemExit(main())
