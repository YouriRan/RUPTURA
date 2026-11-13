"""Export the canonical Lab/CLI configuration and portable cluster launchers."""
import io
import json
from copy import deepcopy
from pathlib import Path
import zipfile

def cluster_bundle(config, directory=None):
    if Path(config.get('Output',{}).get('Directory','results')).is_absolute():
        raise ValueError('Use a relative Output.Directory, such as results, before exporting a portable cluster bundle.')
    from ruptura.process import canonical_chemisorption
    portable = canonical_chemisorption(config)
    dependencies = {}
    def check(value, parent=None):
        if isinstance(value,dict):
            for key,item in value.items():
                if key.lower() == 'filename' and parent == 'mpdsettings' and item:
                    source = Path(item).expanduser()
                    if not source.is_absolute():
                        source = Path(directory or '.') / source
                    source = source.resolve()
                    if not source.is_file():
                        raise ValueError(f'MPD distribution file does not exist: {source}')
                    if source not in dependencies:
                        dependencies[source] = f'mpd/distribution_{len(dependencies)+1}.data'
                    value[key] = dependencies[source]
                elif key in ('ReadColumnFile','FileName') and item:
                    raise ValueError('Cluster export requires inline material/grid data except MPD distribution files, which are included automatically.')
                else:
                    check(item, key.lower())
        elif isinstance(value,list):
            for item in value:check(item,parent)
    check(portable)
    run='''#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
RUPTURA_EXECUTABLE="${RUPTURA_EXECUTABLE:-ruptura}"
"$RUPTURA_EXECUTABLE" simulation.json --validate
exec "$RUPTURA_EXECUTABLE" simulation.json "$@"
'''
    submit='''#!/usr/bin/env bash
#SBATCH --job-name=ruptura
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --output=ruptura-%j.out
#SBATCH --error=ruptura-%j.err
set -euo pipefail
cd -- "${SLURM_SUBMIT_DIR:-$(dirname -- "${BASH_SOURCE[0]}")}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
srun bash run.sh "$@"
'''
    files={'simulation.json':json.dumps(portable,indent=2,allow_nan=False)+'\n','run.sh':run,'submit_hal.sh':submit,'README.txt':'This simulation.json uses the same canonical configuration edited in Lab.\nMPD distribution files, when configured, are included under mpd/; only their FileName paths change for portability.\nBuild RUPTURA on the cluster with native CPU optimization disabled.\nSet RUPTURA_EXECUTABLE to its absolute path, then run bash run.sh locally\nor sbatch --export=ALL,RUPTURA_EXECUTABLE=/path/to/ruptura submit_hal.sh.\nAdjust memory, wall time and site module directives to the cluster.\nThe exporter does not submit jobs. For restart: bash run.sh --restart.\n'}
    stream=io.BytesIO()
    with zipfile.ZipFile(stream,'w',zipfile.ZIP_DEFLATED) as archive:
        for name,text in files.items():
            info=zipfile.ZipInfo(name);info.external_attr=(0o100755 if name.endswith('.sh') else 0o100644)<<16;info.compress_type=zipfile.ZIP_DEFLATED
            archive.writestr(info,text)
        for source,name in dependencies.items():
            archive.write(source,name)
    return stream.getvalue()
