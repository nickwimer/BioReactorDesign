"""Compile actual OF9 blending methods; test the isolated saved-control correction."""
from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sc1_reactor_startup as run
from sc1_reactor_startup_case import case_files, initialized_fields
from sc1_reactor_startup_flow import statistics, flow
import test_sc1_reactor_startup as startup_tests
from test_sc1_reactor_startup import array_field


class Blending(unittest.TestCase):
    def test_actual_linear_and_none_source_functions(self):
        compiler=shutil.which('g++')
        if not compiler:self.skipTest('C++ compiler unavailable')
        folder=run.base.SOURCE/'phaseSystems/BlendedInterfacialModel/blendingMethods'
        functions=[]
        for model in ('linear','noBlending'):
            text=(folder/model/(model+'.C')).read_text()
            anchor='Foam::tmp<Foam::volScalarField> Foam::blendingMethods::'+model+'::f1'
            self.assertEqual(text.count(anchor),1)
            functions.append(text[text.index(anchor):])
        code=r'''
#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <string>
namespace Foam {
using scalar=double;
constexpr double small=1e-15;
constexpr int dimless=0;
struct fvMesh {};
struct dimensionedScalar {
    double value;
    dimensionedScalar(double x):value(x){}
    dimensionedScalar(int, bool x):value(x){}
    operator double() const {return value;}
};
struct phaseModel {
    std::string n; double a; fvMesh m;
    std::string name() const {return n;}
    const fvMesh& mesh() const {return m;}
    operator double() const {return a;}
};
struct volScalarField {
    static double New(const char*,const fvMesh&,dimensionedScalar a) {return a;}
};
template<class T> using tmp=double;
inline double min(double a,double b){return std::min(a,b);}
inline double max(double a,double b){return std::max(a,b);}
struct values {
    std::map<std::string,double> v;
    double operator[](const std::string& s)const{return v.at(s);}
};
namespace blendingMethods {
struct linear {
    values minFullyContinuousAlpha_{{{"gas",.7},{"liquid",.7}}};
    values minPartlyContinuousAlpha_{{{"gas",.3},{"liquid",.3}}};
    tmp<volScalarField> f1(const phaseModel&,const phaseModel&)const;
    tmp<volScalarField> f2(const phaseModel&,const phaseModel&)const;
};
struct noBlending {
    std::string continuousPhase_="none";
    tmp<volScalarField> f1(const phaseModel&,const phaseModel&)const;
    tmp<volScalarField> f2(const phaseModel&,const phaseModel&)const;
};
}}
'''+ '\n'.join(functions)+r'''
int main() {
    using namespace Foam;
    blendingMethods::linear legacy; blendingMethods::noBlending fixed;
    for(int k=0;k<=10000;++k) {
        phaseModel a{"gas",k/10000.,{}}, b{"liquid",1-k/10000.,{}};
        const double before=1-legacy.f1(a,b)-legacy.f2(a,b);
        const double after=1-fixed.f1(a,b)-fixed.f2(a,b);
        assert(std::abs(before)<1e-12);
        assert(after==1);
    }
    // Choosing an actual phase as 'continuous' would disable unordered drag again.
    fixed.continuousPhase_="liquid";
    phaseModel a{"gas",.4,{}},b{"liquid",.6,{}};
    assert(1-fixed.f1(a,b)-fixed.f2(a,b)==0);
}
'''
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'test.C').write_text(code)
            result=subprocess.run([compiler,'-std=c++14','-Wall','-Wextra','-Werror','-Wno-unused-parameter',str(p/'test.C'),'-o',str(p/'test')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            subprocess.run([str(p/'test')],check=True)

    def test_original_recipe_unchanged_and_only_drag_override_added(self):
        old=case_files(8);new=case_files(8,segregated_drag=True)
        # Canonical recipe captured and compared to the committed e586c47 source.
        self.assertEqual(hashlib.sha256(json.dumps(old,sort_keys=True).encode()).hexdigest(),
            '7cb903ffa76fbca12fa108a7333705fe8c00faef304b32dff28a99865b030403')
        self.assertEqual([k for k in old if old[k]!=new[k]],['constant/phaseProperties'])
        self.assertEqual(new['constant/phaseProperties'].replace(
            ' }\n    drag { type none; continuousPhase none; }\n}', ' } }'),old['constant/phaseProperties'])
        templates=(run.base.SOURCE/'phaseSystems/phaseSystem/phaseSystemTemplates.C').read_text()
        self.assertIn('blendingMethods_.found(modelName)',templates)
        blended=(run.base.SOURCE/'phaseSystems/BlendedInterfacialModel/BlendedInterfacialModel.C').read_text()
        self.assertIn('(scalar(1) - f1() - f2())*(model_().*method)',blended)


class SavedControl(unittest.TestCase):
    def fixture(self,root):
        case=root/'quiescent'
        for rel,text in case_files(8).items():
            p=case/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
        for name in ('checkMesh','writeCellCentres','writeCellVolumes'):(case/('log.'+name)).write_text('saved evidence')
        for name in ('startup_manifest.json','quiescent/result.json','quiescent/log.solver'):(root/name).write_text('saved')
        return dict(variant=run.TAG,root=str(root),ranks=8,inputs=run.base.case_input_hashes(case),
            build=str(root.parent/'build'),archive=str(root.parent/'archive'),regression_root=str(root.parent/'rates'),
            assets={},initial={},mesh_hashes={},regression_sha256='proof')

    def test_copy_only_initial_inputs_and_preserve_failed_control(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old=p/'old';old.mkdir();m=self.fixture(old)
            (old/'quiescent/0.01').mkdir();(old/'quiescent/0.01/U.gas').write_text('completed field')
            before=run.tree_hashes(old)
            report=dict(log_sha256=run.base.sha(old/'quiescent/log.solver'),screen_pass=True,startup_pass=False)
            with patch.object(run,'saved_startup',return_value=(m,report)),patch.object(run,'require_allocation'):
                out=run.prepare_drag(old,p/'new')
                new=json.loads((out/'startup_manifest.json').read_text())
                self.assertEqual(new['variant'],run.DRAG_TAG)
                self.assertEqual(new['baseline'],run.baseline_identity(old,report))
                changed=[k for k in m['inputs'] if m['inputs'][k]!=new['inputs'][k]]
                self.assertEqual(changed,['constant/phaseProperties'])
                self.assertFalse((out/'quiescent/0.01').exists())
                self.assertFalse((out/'quiescent/log.solver').exists())
                self.assertEqual(run.tree_hashes(old),before)
                with self.assertRaises(ValueError):run.prepare_drag(old,p/'new')

    def test_saved_failed_quiescence_is_valid_control_but_tampering_is_not(self):
        with tempfile.TemporaryDirectory() as d,patch.object(run,'CELLS',4):
            root=Path(d);case=root/'quiescent';case.mkdir()
            m=startup_tests.Reports().setup_case(case);m.update(variant=run.TAG,build='/private',assets={},ranks=8)
            u=np.zeros((4,3));u[0,0]=.002
            (case/'0.01/U.liquid').write_text(array_field('U.liquid',u))
            report=run.analyze(case,m);self.assertFalse(report['startup_pass'])
            (case/'result.json').write_text(json.dumps(report))
            (case/'attempt.json').write_text(json.dumps(dict(build='/private',assets={},ranks=8)))
            (root/'startup_manifest.json').write_text(json.dumps(dict(variant=run.TAG)))
            m['inputs']=run.base.case_input_hashes(case)
            with patch.object(run,'verify',return_value=m):
                self.assertEqual(run.saved_startup(root)[1],report)
                (case/'result.json').write_text('{}')
                with self.assertRaisesRegex(ValueError,'disagree'):run.saved_startup(root)
            (root/'startup_manifest.json').write_text(json.dumps(dict(variant=run.DRAG_TAG)))
            with self.assertRaisesRegex(ValueError,'original v1'):run.saved_startup(root)

    def test_new_report_requires_runtime_drag_selection(self):
        with tempfile.TemporaryDirectory() as d,patch.object(run,'CELLS',4):
            p=Path(d);m=startup_tests.Reports().setup_case(p);m['variant']=run.DRAG_TAG
            with self.assertRaisesRegex(ValueError,'runtime selection'):run.analyze(p,m)
            log=p/'log.solver';log.write_text('Selecting drag blending method: none\n'+log.read_text())
            report=run.analyze(p,m)
            self.assertTrue(report['startup_pass']);self.assertEqual(report['suite'],run.DRAG_TAG)

    def test_verifier_rejects_changes_beyond_blending(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old=p/'old';old.mkdir();m=self.fixture(old);case=old/'quiescent'
            c=np.array([[0,y,0] for y in (.999,1.01,1.022,1.033)]);v=np.full(4,.5)
            fields,m['initial']=initialized_fields(c,v)
            for rel,text in fields.items():(case/rel).write_text(text)
            (case/'0/C').write_text(array_field('C',c));(case/'0/V').write_text(array_field('V',v))
            (case/'log.checkMesh').write_text('Mesh OK.')
            m['inputs']=run.base.case_input_hashes(case)
            (p/'archive/constant/polyMesh').mkdir(parents=True)
            report={'log_sha256':run.base.sha(case/'log.solver')}
            with (patch.object(run,'saved_startup',return_value=(m,report)),patch.object(run,'require_allocation'),
                  patch.object(run,'CELLS',4),patch.object(run,'regression',return_value=({},'proof')),
                  patch.object(run.base,'verify_build',return_value={'assets':{}})):
                out=run.prepare_drag(old,p/'new')
                self.assertEqual(run.verify(out)['variant'],run.DRAG_TAG)
                control=out/'quiescent/system/controlDict'
                control.write_text(control.read_text().replace('deltaT 5e-05;','deltaT 1e-05;'))
                manifest=json.loads((out/'startup_manifest.json').read_text())
                manifest['inputs']['system/controlDict']=run.base.sha(control)
                (out/'startup_manifest.json').write_text(json.dumps(manifest))
                with self.assertRaisesRegex(ValueError,'differs from its recipe'):run.verify(out)


class Flow(unittest.TestCase):
    def test_counterflow_can_have_small_mass_velocity(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'0').mkdir();(p/'0.01').mkdir()
            data={'0/C':[[0,1,0],[0,2,0]],'0/V':[1,1],
                '0.01/alpha.liquid':[.5,1],'0.01/alpha.gas':[.5,0],
                '0.01/U.gas':[[0,1,0],[0,0,0]],'0.01/U.liquid':[[0,-.001,0],[0,0,0]]}
            for rel,v in data.items():(p/rel).write_text(array_field(Path(rel).name,v))
            r=flow(p,2)
            self.assertEqual(r['massVelocityRms_m_s'],0)
            self.assertAlmostEqual(r['interfaceSlip']['rms_m_s'],1.001)
            self.assertEqual(r['phases']['liquid']['bulk']['max_m_s'],0)
            self.assertEqual(r['phases']['gas']['present']['peak']['alphaLiquid'],.5)
            self.assertEqual(r['phases']['gas']['bulk'],{'cells':0})


if __name__=='__main__':unittest.main()
