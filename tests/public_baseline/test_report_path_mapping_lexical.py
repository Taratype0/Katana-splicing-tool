"""Metadata projection is lexical; actual file reads retain symlink bounds."""
from pathlib import Path
from unittest.mock import patch

import pytest

from katana_next.reports.path_mapping import PathMappings


def test_metadata_projection_never_resolves_or_stats_paths(tmp_path):
    target = tmp_path / 'consumer'
    mapper = PathMappings([{'source_prefix': 'E:/study', 'target_prefix': str(target)}])
    original = {'annotation': {'inputs': ['E:/study/genes/../gene-sets/GO.tsv',
                                          r'e:\STUDY\tables\KEGG.tsv']},
                'source_sha256': 'literal-hash', 'sample': 'sample1',
                'nested': [{'path': 'E:/study/tables/a.tsv'}]}
    with patch.object(Path, 'resolve', side_effect=AssertionError('metadata resolve I/O')), \
         patch.object(Path, 'stat', side_effect=AssertionError('metadata stat I/O')), \
         patch.object(Path, 'lstat', side_effect=AssertionError('metadata lstat I/O')):
        projected = mapper.project(original)
        repeated = mapper.project(original)
    assert projected == repeated
    assert projected['annotation']['inputs'] == [str(target / 'gene-sets/GO.tsv'),
                                                 str(target / 'tables/KEGG.tsv')]
    assert projected['nested'][0]['path'] == str(target / 'tables/a.tsv')
    assert original['nested'][0]['path'] == 'E:/study/tables/a.tsv'
    assert projected['source_sha256'] == 'literal-hash'


@pytest.mark.parametrize('path', ['E:/study/../../outside', 'E:/study/sub/../../../outside'])
def test_lexical_parent_escape_is_rejected_without_resolve(tmp_path, path):
    mapper = PathMappings([{'source_prefix': 'E:/study', 'target_prefix': str(tmp_path / 'consumer')}])
    with patch.object(Path, 'resolve', side_effect=AssertionError('metadata resolve I/O')):
        with pytest.raises(ValueError, match='escapes'):
            mapper.value(path)


def test_windows_source_case_and_longest_prefix_component_guard(tmp_path):
    main, nested = tmp_path / 'main', tmp_path / 'nested'
    mapper = PathMappings([{'source_prefix': 'E:/Study', 'target_prefix': str(main)},
                           {'source_prefix': 'E:/Study/nested', 'target_prefix': str(nested)}])
    assert mapper.value(r'e:\STUDY\nested\a.tsv') == str(nested / 'a.tsv')
    assert mapper.value('E:/Study/nested-other/a.tsv') == str(main / 'nested-other/a.tsv')
    assert mapper.value('E:/Study-other/a.tsv') == 'E:/Study-other/a.tsv'
    with pytest.raises(ValueError, match='unique'):
        PathMappings([{'source_prefix': 'E:/Study', 'target_prefix': str(main)},
                      {'source_prefix': 'e:/study', 'target_prefix': str(nested)}])


def test_posix_source_prefix_remains_case_sensitive(tmp_path):
    mapper = PathMappings([{'source_prefix': '/producer/Study', 'target_prefix': str(tmp_path)}])
    assert mapper.value('/producer/Study/tables/a.tsv') == str(tmp_path / 'tables/a.tsv')
    assert mapper.value('/producer/study/tables/a.tsv') == '/producer/study/tables/a.tsv'


def test_actual_file_resolution_is_cached_per_mapper_and_base(tmp_path):
    target = tmp_path / 'consumer'
    target.mkdir()
    table = target / 'a.tsv'
    table.write_text('gene\tFDR\ng1\t0.001\n', encoding='utf-8')
    mapper = PathMappings([{'source_prefix': 'E:/study', 'target_prefix': str(target)}])
    native_resolve = Path.resolve
    calls = []
    def observed(path, *args, **kwargs):
        calls.append(str(path))
        return native_resolve(path, *args, **kwargs)
    with patch.object(Path, 'resolve', observed):
        first = mapper.resolve('E:/study/a.tsv', tmp_path)
        assert mapper.resolve('E:/study/a.tsv', tmp_path) == first
        assert mapper.resolve('E:/study/a.tsv', str(tmp_path)) == first
    assert first == table
    assert calls == [str(table)]
    assert table.read_text('utf-8') == 'gene\tFDR\ng1\t0.001\n'


@pytest.mark.parametrize('spelling', ['producer', 'projected', 'direct-consumer'])
def test_actual_symlink_escape_is_rejected_after_metadata_projection(tmp_path, spelling):
    target, outside = tmp_path / 'consumer', tmp_path / 'outside'
    target.mkdir()
    outside.mkdir()
    (outside / 'result.tsv').write_text('private original\n', encoding='utf-8')
    try:
        (target / 'escape').symlink_to(outside, target_is_directory=True)
    except (OSError, NotImplementedError) as error:
        pytest.skip('Native symlink creation unavailable: ' + str(error))
    mapper = PathMappings([{'source_prefix': 'E:/study', 'target_prefix': str(target)}])
    producer = 'E:/study/escape/result.tsv'
    value = producer
    if spelling == 'projected':
        value = mapper.project({'path': producer})['path']
        assert value == str(target / 'escape/result.tsv')
    elif spelling == 'direct-consumer':
        value = str(target / 'escape/result.tsv')
    with pytest.raises(ValueError, match='symbolic link'):
        mapper.resolve(value, tmp_path)
    assert (outside / 'result.tsv').read_text('utf-8') == 'private original\n'


def test_actual_in_target_symlink_resolves_to_readable_file(tmp_path):
    target = tmp_path / 'consumer'
    target.mkdir()
    actual = target / 'actual.tsv'
    actual.write_text('native complete result\n', encoding='utf-8')
    try:
        (target / 'alias.tsv').symlink_to(actual)
    except (OSError, NotImplementedError) as error:
        pytest.skip('Native symlink creation unavailable: ' + str(error))
    mapper = PathMappings([{'source_prefix': 'E:/study', 'target_prefix': str(target)}])
    assert mapper.resolve('E:/study/alias.tsv', tmp_path) == actual
