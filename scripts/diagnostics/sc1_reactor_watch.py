"""Stop only this new startup step on numerical failure; keep its log and fields."""
import json
import math
import os
import re
import signal
import subprocess
import time


class Monitor:
    def __init__(self):
        self.line=0;self.time=None;self.reason=None

    def feed(self,text):
        self.line+=1;text=text.strip()
        if self.reason:return self.reason
        match=re.fullmatch(r'Time = (\S+)',text)
        if match:self.time=match[1]
        reason=None
        match=re.fullmatch(r'Courant Number mean: (\S+) max: (\S+)',text)
        if match:
            try:values=[float(x) for x in match.groups()]
            except ValueError:values=[math.nan]
            if not all(math.isfinite(x) for x in values):reason='nonfinite Courant number'
            elif values[-1]>.2:reason='mixture Courant number exceeded 0.2'
        match=re.fullmatch(r'(gas|liquid) fraction, min, max = (\S+) (\S+) (\S+)',text)
        if match:
            try:values=[float(x) for x in match.groups()[1:]]
            except ValueError:values=[math.nan]
            if not all(math.isfinite(x) for x in values):reason='nonfinite alpha subcycle'
            elif values[1]<-1e-8 or values[2]>1+1e-8:reason='alpha subcycle exceeded [0,1] by more than 1e-8'
        if text.startswith('limitPressure:'):reason='pressure limiter activated in quiescent vessel'
        lower=text.lower()
        banner=re.match(r'^sigfpe\s*:\s*enabling floating point exception trapping\b',lower)
        if not banner and any(s in lower for s in ('floating point exception','sigfpe::sighandler','foam fatal','segmentation fault')):
            reason='solver failure'
        if (text.startswith('SC1_') or 'Solving for p_rgh,' in text) and re.search(r'(?<![a-z])(?:nan|[+-]?inf)(?![a-z])',lower):
            reason='nonfinite audit or pressure solve'
        if reason:self.reason=dict(reason=reason,line=self.line,lastPrintedTime=self.time,text=text)
        return self.reason


def stop(process):
    """Signal this child session only; srun forwards termination to its own step."""
    try:os.killpg(process.pid,signal.SIGTERM)
    except ProcessLookupError:pass
    try:process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        try:os.killpg(process.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        process.wait(timeout=10)


def execute(command,cwd,log,env,*,wall_seconds=1800,poll_seconds=.2):
    if wall_seconds<=0 or poll_seconds<=0:raise ValueError('Invalid startup monitor timing')
    print('Executing with startup monitor: '+' '.join(map(str,command)),flush=True)
    monitor=Monitor();start=time.monotonic();reason=None;pending=''
    with log.open('x') as output:
        process=subprocess.Popen(list(map(str,command)),cwd=cwd,env=env,stdout=output,
                                 stderr=subprocess.STDOUT,start_new_session=True)
        try:
            with log.open(errors='replace') as reader:
                while True:
                    # Poll before draining: final writes are then visible even on a fast exit.
                    code=process.poll();pending+=reader.read()
                    lines=pending.split('\n');pending=lines.pop()
                    if code is not None and pending:lines.append(pending);pending=''
                    for line in lines:
                        if monitor.feed(line):reason=monitor.reason;break
                    if reason:break
                    if code is not None:break
                    if time.monotonic()-start>=wall_seconds:
                        reason=dict(reason='startup exceeded wall-clock limit',wall_seconds=wall_seconds);break
                    time.sleep(poll_seconds)
            if reason:
                stop(process)
                reason.update(returncode=process.returncode,meaning='Stopped failed/incomplete diagnostic; no valid restart implied.')
                with (cwd/'guard_stop.json').open('x') as f:json.dump(reason,f,indent=2);f.write('\n')
                print(json.dumps(reason,indent=2),flush=True)
                raise ValueError('Startup stopped: '+reason['reason']+'; preserve '+str(log))
            if process.returncode:
                raise ValueError(f'Command failed ({process.returncode}); preserve {log}')
        except BaseException:
            if process.poll() is None:stop(process)
            raise
