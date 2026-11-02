"""Installed worker for production whole-study CoSplice; no experimental engines."""
from pathlib import Path
import json,sys
import pandas as pd
from product import run,write_report
from enrichment42 import resources,calculate,ROOT
target=Path(sys.argv[1]);cfg=json.loads((target/'config.json').read_text())
run(cfg,target)
summary=json.loads((target/'summary.json').read_text())
# Gene-set databases are explicit user resources, never bundled study caches.
(target/'report_revision.json').write_text(json.dumps({'schema':'study-report-2.9','presentation_revision':39}))
write_report(target,summary)
if not (target/'cosplice_chapter6.html').exists():
    from chapter6_42 import build
    build(target,target)
print(json.dumps({'status':'completed','samples':summary['samples'] if 'samples' in summary else summary['study']['samples'],'network_events':summary['network_events']}))
