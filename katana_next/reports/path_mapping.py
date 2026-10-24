"""Explicit consumer-side path bindings; producer files and identities stay intact."""
from pathlib import Path
import copy
import re


class PathMappings:
    def __init__(self, mappings=(), *, base=None):
        self.rows = []
        for row in mappings or ():
            source = str(row['source_prefix']).replace('\\', '/').rstrip('/')
            target = Path(row['target_prefix']).expanduser()
            if not target.is_absolute(): target = Path(base or Path.cwd()) / target
            if not source: raise ValueError('A path mapping needs a nonempty source prefix')
            self.rows.append((source, target.resolve(), bool(re.match(r'^[A-Za-z]:/', source))))
        if len({s.casefold() if w else s for s, _, w in self.rows}) != len(self.rows):
            raise ValueError('Path mapping source prefixes must be unique')
        self.rows.sort(key=lambda r: len(r[0]), reverse=True)

    def value(self, value):
        if not isinstance(value, str): return value
        normalized = value.replace('\\', '/')
        for prefix, target, windows in self.rows:
            actual, wanted = (normalized.casefold(), prefix.casefold()) if windows else (normalized, prefix)
            if actual != wanted and not actual.startswith(wanted + '/'): continue
            tail = normalized[len(prefix):].lstrip('/')
            path = (target / tail).resolve()
            if not path.is_relative_to(target): raise ValueError('Mapped source path escapes its declared target')
            return str(path)
        return value

    def project(self, value):
        """Return a runtime projection only; never rewrite a source JSON file."""
        if isinstance(value, dict): return {k:self.project(v) for k,v in value.items()}
        if isinstance(value, list): return [self.project(v) for v in value]
        return self.value(value)

    def resolve(self, value, base):
        path = Path(self.value(str(value))).expanduser()
        return (path if path.is_absolute() else Path(base) / path).resolve()


def identity_projection(project):
    """Biological design identity excludes machine-specific resource locations."""
    return {
        'samples': [{k:s.get(k) for k in ('sample_id','biological_sample_id','condition','batch','attributes')}
                    for s in project.get('samples', [])],
        'contrasts': copy.deepcopy(project.get('contrasts', [])),
        'genome_build': project.get('annotation', {}).get('genome_build'),
    }
