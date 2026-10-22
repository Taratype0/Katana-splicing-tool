from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from scripts import public_motiftrace5712_20261016 as public
from katana_next.science.motiftrace_full_execution import profiles
from katana_next.science.motiftrace_exact import ExactSignedRank, full_family_adjust
from katana_next.science.motiftrace_library_execution import DatasetPlan, PackedPWM


def snapshot(tmp_path):
    (tmp_path / 'snapshots').mkdir()
    for name in ('motiftrace_custom.py', 'motiftrace_full_execution.py',
                 'motiftrace_library_execution.py', 'motiftrace_exact.py'):
        shutil.copyfile(public.REPO / 'katana_next/science' / name, tmp_path / 'snapshots' / name)
    shutil.copyfile(public.REPO / 'scripts/run_motiftrace_exact5712_20261016.py',
                    tmp_path / 'snapshots/run_motiftrace_exact5712_20261016.py')
    return public.kernels(tmp_path)


def test_original_AST_exact_scores_and_profiles_are_identical(tmp_path):
    functions = snapshot(tmp_path)
    for word in ('ACGTACGT', 'CTCTCTCTNCTCTCT', 'A'*60, 'A'*3):
        expected = profiles(word, '[CT][CT]T')
        actual = functions['profiles'](word, '[CT][CT]T')
        for width in expected:
            for first, second in zip(expected[width], actual[width]):
                np.testing.assert_array_equal(first, second)
    values = np.random.default_rng(17).integers(-3, 4, (43, 12)).astype(float)
    values[2, :3] = np.nan
    for expected, actual in zip(ExactSignedRank().evaluate(values), functions['ExactSignedRank']().evaluate(values)):
        np.testing.assert_array_equal(expected, actual)


def test_original_AST_dataset_and_continuous_PWM_are_identical(tmp_path):
    functions = snapshot(tmp_path)
    records = []
    pairs = []
    for i in range(25):
        pairs.append({'event_type': 'SE', 'direction': 'inclusion_up', 'global_pair_index': i,
                      'foreground_event': f'f{i}', 'background_event': f'b{i}'})
        for role, word in (('f', 'CTCTCTAACG'*9), ('b', 'ACGTACGTAA'*9)):
            for region in ('affected_inclusion', 'inclusion_junction_1_acceptor_intron'):
                records.append({'event_id': f'{role}{i}', 'region': region, 'sequence': word,
                                'event_type': 'SE', 'layer': 'affected_sequence'})
    rna, pair_table = pd.DataFrame(records), pd.DataFrame(pairs)
    original, actual = DatasetPlan(rna, pair_table), functions['DatasetPlan'](rna, pair_table)
    matrix = np.array([[.6, .1, .2, .1]]*4)
    packed1, packed2 = PackedPWM(rna.sequence.tolist()), functions['PackedPWM'](rna.sequence.tolist())
    packed1.activate(matrix)
    packed2.activate(matrix)
    for model in ({'kind': 'consensus', 'pattern': 'CTCT'}, {'kind': 'PWM'}):
        for first, second in zip(original.datasets(model, packed1), actual.datasets(model, packed2)):
            assert first['meta'] == second['meta']
            np.testing.assert_array_equal(first['differences'], second['differences'])


def test_disk_format_retains_string_identity_and_all_numeric_rows(tmp_path):
    frame = pd.DataFrame({'method': ['count', 'PWM'], 'P_exact': [1., 2.**-100], 'n': [20, 100]})
    public.write_arrays(tmp_path / 'test.npz', frame)
    actual = public.load_arrays(tmp_path / 'test.npz')
    pd.testing.assert_frame_equal(actual, frame)


def test_complete_family_adjustments_preserve_ties_and_match_original():
    values = np.r_[np.repeat([2.**-40, .001, 1.], [15, 7, 20]), [.03, .07]]
    for expected, actual in zip(full_family_adjust(values), public.adjust(values)):
        np.testing.assert_allclose(expected, actual, rtol=2e-14, atol=1e-16)


def test_window_coordinates_are_RNA_oriented_not_invented_genomic():
    frame = pd.DataFrame({'window_nt': [50, 50, 20, 0], 'window_start0': [0, 3, 4, -1],
                          'coordinate_system': ['masked_acceptor_3prime_minus20', 'masked_donor_5prime_plus6',
                                                'affected_RNA_5prime', 'affected_RNA_5prime']})
    actual = public.physical_coordinates(frame)
    np.testing.assert_array_equal(actual.RNA_start_nt[:3], [-70, 9, 4])
    np.testing.assert_array_equal(actual.RNA_end_nt[:3], [-20, 59, 24])
    assert pd.isna(actual.RNA_start_nt.iloc[-1])


def test_source_reader_accepts_UTF8_BOM_without_changing_definition(tmp_path):
    source = tmp_path / 'with_bom.py'
    source.write_bytes(b'\xef\xbb\xbfdef original(value):\n    return value + 7\n')
    before = source.read_bytes()
    api = public.named_source(source, {'original'}, {})
    assert api['original'](12) == 19
    assert source.read_bytes() == before


def test_actual_archived_author_API_counts_and_rank_sum_execute(tmp_path):
    root = public.LOCAL / 'bundle'
    api = public.author_kernel(root)
    for group in ('up', 'dn', 'bg'):
        for region in api['region']:
            (tmp_path / f'{group}.{region}.fasta').write_text(
                ''.join(f'>runtime_{i}\n' + 'CTACGT'*50 + '\n' for i in range(24)),
                encoding='utf-8')
    api.update(totalExonCount={'up': 24, 'dn': 24, 'bg': 24}, fastaPath=str(tmp_path))
    counts = api['countMotif']({}, 'runtime', 'CT', ['CT'])
    assert counts['up'] == counts['dn'] == counts['bg']
    assert sum(sum(values) for region in counts['up'].values() for values in region.values()) > 0
    probabilities = api['computePValues'](counts['up'], counts['bg'], {})
    assert all(0 <= p <= 1 for region in probabilities.values() for p in region.values())
