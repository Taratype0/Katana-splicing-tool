"""Standalone whole-study CoSplice. Uses the public Katana workflow service."""
import argparse
import json
from pathlib import Path
import sys

from .study_api import CoSpliceStudy


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', required=True, type=Path, help='Katana project manifest')
    commands = parser.add_subparsers(dest='action', required=True)
    for name in ('plan', 'run'):
        sub = commands.add_parser(name)
        sub.add_argument('--parameters', type=Path, help='Optional parameter JSON')
        sub.add_argument('--paired', action='store_true', default=None,
                         help='Explicit repeated-individual design; otherwise use declared comparison design')
        if name == 'run':
            sub.add_argument('--force', action='store_true')
    commands.add_parser('artifacts')
    args = parser.parse_args(argv)
    try:
        parameters = json.loads(args.parameters.read_text(encoding='utf8')) if getattr(args, 'parameters', None) else {}
        if not isinstance(parameters, dict):
            raise ValueError('Parameters must be a JSON object')
        if getattr(args, 'paired', None) is not None:
            if parameters.get('paired') is False:
                raise ValueError('--paired conflicts with paired=false in parameter file')
            parameters['paired'] = True
        with CoSpliceStudy(args.project) as client:
            if args.action == 'plan':
                result = client.plan(parameters)
            elif args.action == 'run':
                result = client.run(parameters, force=args.force).model_dump(mode='json')
            else:
                result = [a.model_dump(mode='json') for a in client.artifacts()]
        print(json.dumps(result, indent=2, ensure_ascii=True))
        if args.action == 'plan' and not result['ready']:
            return 2
        if args.action == 'run' and result['status'] not in {'complete', 'cache_hit'}:
            return 1
        return 0
    except (OSError, ValueError, TypeError, KeyError, RuntimeError, ImportError) as exc:
        print(json.dumps({'error': str(exc)}), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
