"""Explicit consumer-side path bindings; producer files and identities stay intact."""
from pathlib import Path
import copy
import os
import re


class PathMappings:
    def __init__(self, mappings=(), *, base=None):
        self.rows = []
        self._values = {}
        self._resolved = {}
        self._value_targets = {}
        self._projected_targets = {}
        for row in mappings or ():
            source = str(row['source_prefix']).replace('\\', '/').rstrip('/')
            target = Path(row['target_prefix']).expanduser()
            if not target.is_absolute(): target = Path(base or Path.cwd()) / target
            if not source: raise ValueError('A path mapping needs a nonempty source prefix')
            self.rows.append((source, target.resolve(), bool(re.match(r'^[A-Za-z]:/', source))))
        if len({s.casefold() if w else s for s, _, w in self.rows}) != len(self.rows):
            raise ValueError('Path mapping source prefixes must be unique')
        self.rows.sort(key=lambda r: len(r[0]), reverse=True)
        self._target_boundaries = sorted({target for _, target, _ in self.rows},
                                         key=lambda target: len(str(target)), reverse=True)

    def value(self, value):
        if not isinstance(value, str): return value
        if value in self._values: return self._values[value]
        normalized = value.replace('\\', '/')
        for prefix, target, windows in self.rows:
            actual, wanted = (normalized.casefold(), prefix.casefold()) if windows else (normalized, prefix)
            if actual != wanted and not actual.startswith(wanted + '/'): continue
            tail = normalized[len(prefix):].lstrip('/')
            # Metadata can contain millions of producer paths that are not
            # opened. Normalize their consumer spelling without stat/lstat or
            # a symlink walk (especially expensive on mounted filesystems).
            path = Path(os.path.normpath(str(target / tail)))
            if not path.is_relative_to(target): raise ValueError('Mapped source path escapes its declared target')
            self._values[value] = str(path)
            self._value_targets[value] = target
            prior = self._projected_targets.get(str(path))
            if prior is None or target.is_relative_to(prior):
                self._projected_targets[str(path)] = target
            return self._values[value]
        self._values[value] = value
        return value

    def project(self, value):
        """Return a runtime projection only; never rewrite a source JSON file."""
        if isinstance(value, dict): return {k:self.project(v) for k,v in value.items()}
        if isinstance(value, list): return [self.project(v) for v in value]
        return self.value(value)

    def resolve(self, value, base):
        key = (str(value), str(base))
        if key in self._resolved: return self._resolved[key]
        path = Path(self.value(str(value))).expanduser()
        path = Path(os.path.normpath(str(path if path.is_absolute() else Path(base) / path)))
        # project() may already have replaced the producer prefix. Keep that
        # declared boundary when the projected value is later opened, and also
        # recognize direct consumer paths into an explicitly mapped target.
        boundary = self._value_targets.get(str(value), self._projected_targets.get(str(path)))
        if boundary is None:
            boundary = next((target for target in self._target_boundaries
                             if path.is_relative_to(target)), None)
        resolved = path.resolve()
        if boundary is not None and not resolved.is_relative_to(boundary):
            raise ValueError('Mapped source path escapes its declared target through a symbolic link')
        self._resolved[key] = resolved
        return self._resolved[key]


def identity_projection(project):
    """Biological design identity excludes machine-specific resource locations."""
    return {
        'samples': [{k:s.get(k) for k in ('sample_id','biological_sample_id','condition','batch','attributes')}
                    for s in project.get('samples', [])],
        'contrasts': copy.deepcopy(project.get('contrasts', [])),
        'genome_build': project.get('annotation', {}).get('genome_build'),
    }
