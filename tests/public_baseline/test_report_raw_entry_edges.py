"""Report-only edge cases observed after actual FASTQ/BAM calculation."""
import pandas as pd
import pytest

from katana_next.flow.cli import required_figure_chapters
from katana_next.reports.engine import Builder


def test_complete_table_preserves_signed_infinity_missing_and_raw_file(tmp_path):
    source = tmp_path / 'native.tsv'
    original = b'gene\teffect\nA\tinf\nB\t-inf\nC\t\nD\t0\n'
    source.write_bytes(original)
    builder = Builder({}, tmp_path / 'report', 'stable_significant')
    frame = pd.read_csv(source, sep='\t')
    record = {'name': 'Native effect', 'source': 'outputs/native.tsv'}
    output = builder.output / builder.table_browser(frame, record)
    content = output.read_text('utf-8')
    # Browser JSON retains the sign and distinguishes missing from zero.
    assert '"Infinity"' in content and '"-Infinity"' in content
    assert '["C", null]' in content and '["D", 0.0]' in content
    assert 'original downloaded table is unchanged' in content
    assert source.read_bytes() == original


def single_contrast_plan():
    return {'contrast_ids': ['B_vs_A'],
            'applicability': {'cross_comparison': {'status': 'not_applicable', 'contrast_ids': ['B_vs_A']}},
            'configuration': {'design_factors': ['condition']}}


def test_declared_single_contrast_omits_only_meaningless_overlap_figure():
    assert required_figure_chapters(single_contrast_plan()) == {1, 2, 3, 5, 6, 7}


@pytest.mark.parametrize('case', ['two_contrasts', 'interactions', 'multifactor', 'missing_declaration', 'wrong_scope'])
def test_chapter_four_cannot_be_waived_for_other_designs(case):
    plan = single_contrast_plan()
    if case == 'two_contrasts':
        plan['contrast_ids'].append('C_vs_A')
    elif case == 'interactions':
        plan['configuration']['interactions'] = [{'id': 'treatment_by_genotype'}]
    elif case == 'multifactor':
        plan['configuration']['design_factors'].append('genotype')
    elif case == 'missing_declaration':
        plan.pop('applicability')
    else:
        plan['applicability']['cross_comparison']['contrast_ids'] = ['X_vs_Y']
    assert required_figure_chapters(plan) == set(range(1, 8))
