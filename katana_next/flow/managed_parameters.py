"""Use the native workflow's public parameters in the storage-aware launcher."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys


def schema():
    root = Path(__file__).parent / '_pipeline'
    if not (root / 'public_parameters.json').is_file():
        root = Path(__file__).resolve().parents[2] / 'pipeline'
    return json.loads((root / 'public_parameters.json').read_text(encoding='utf-8'))


def convert(value, kind):
    if value is None:
        return None
    if kind == 'bool':
        if isinstance(value, bool):
            return value
        if isinstance(value, str) and value.lower() in ('true', 'false'):
            return value.lower() == 'true'
        raise ValueError('expected true or false')
    if kind in ('int', 'number'):
        if isinstance(value, bool):
            raise ValueError('expected a number, not a boolean')
        parsed = float(value)
        if not math.isfinite(parsed) or (kind == 'int' and not parsed.is_integer()):
            raise ValueError('expected a finite ' + kind)
        return int(parsed) if kind == 'int' else parsed
    if kind == 'path':
        return Path(value).expanduser().resolve()
    if kind == 'map':
        parsed = json.loads(value) if isinstance(value, str) else value
        if not isinstance(parsed, dict):
            raise ValueError('expected a mapping')
        return parsed
    if kind == 'selection' and isinstance(value, list):
        if any(not isinstance(v, str) for v in value):
            raise ValueError('expected a list of contrast IDs')
        return value
    if not isinstance(value, str):
        raise ValueError('expected text')
    return value


class PublicParser(argparse.ArgumentParser):
    """Parameter file < explicit CLI, with identical public option names."""

    def parse_args(self, args=None, namespace=None):
        argv = list(sys.argv[1:] if args is None else args)
        probe = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
        probe.add_argument('--params-file', '-params-file', type=Path)
        supplied, _ = probe.parse_known_args(argv)
        declared = schema()
        values = {}
        if supplied.params_file:
            content = supplied.params_file.read_text(encoding='utf-8-sig')
            try:
                document = json.loads(content)
            except json.JSONDecodeError:
                import yaml
                document = yaml.safe_load(content)
            document = {} if document is None else document
            if not isinstance(document, dict):
                self.error('The parameter file must contain a mapping')
            unknown = set(document) - set(declared)
            if unknown:
                self.error('Unknown public parameters: ' + ', '.join(sorted(map(str, unknown))))
            for key, value in document.items():
                if value is not None:
                    try:
                        values[key] = convert(value, declared[key]['type'])
                    except (TypeError, ValueError) as exc:
                        self.error(f'Invalid {key}: {exc}')
        target = namespace or argparse.Namespace()
        for key, value in values.items():
            setattr(target, key, value)
        result = super().parse_args(argv, target)
        if result.mode not in ('complete', 'downstream-only'):
            self.error('Managed cleanup supports complete or downstream-only; use nextflow run for upstream-only')
        if not result.outdir:
            self.error('Supply --outdir or outdir in the parameter file')
        # Routing fields are controlled by the two native stage invocations.
        # retain/lean in a shared native template never changes cleanup policy.
        result.native_parameters = {}
        for key in self.native_keys:
            value = vars(result).pop(key, None)
            if value is not None:
                if declared[key]['type'] == 'path':
                    value = str(Path(value).resolve(strict=True))
                result.native_parameters[key] = value
        return result


def add_native_parameters(parser):
    declared = schema()
    present = {action.dest for action in parser._actions}
    parser.native_keys = set(declared) - present - {'storage_mode', 'publish_dir_mode'}
    for key in sorted(parser.native_keys):
        kind = declared[key]['type']
        def parse(value, kind=kind):
            try:
                return convert(value, kind)
            except (TypeError, ValueError) as exc:
                raise argparse.ArgumentTypeError(str(exc)) from exc
        options = {'type': parse, 'default': None}
        if kind == 'bool':
            options.update(nargs='?', const=True)
        parser.add_argument('--' + key, **options)


def native_paths(parameters):
    declared = schema()
    return [Path(value) for key, value in parameters.items()
            if declared[key]['type'] == 'path' and value is not None]
