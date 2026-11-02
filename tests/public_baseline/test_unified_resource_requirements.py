"""MotifTrace's own scanner does not require the optional CLIP/FIMO path."""
from katana_next.flow.selection import DEFAULTS, resource_kinds


def config(**enabled):
    return {'analyses': {**dict.fromkeys(DEFAULTS, False), **enabled}}


def test_custom_eventgraph_library_needs_no_second_database_or_fimo():
    settings = config(motiftrace=True)
    settings['parameters'] = {'prioritize.motiftrace_consequences': {'custom_motif_path': '/data/models.meme'}}
    assert resource_kinds(settings) == set()


def test_registered_motif_library_and_optional_binding_have_distinct_resources():
    assert resource_kinds(config(motiftrace=True)) == {'rbp_motif_database'}
    assert resource_kinds(config(motif_binding=True)) == {'rbp_motif_database', 'fimo_runtime', 'rbp_binding_peaks'}
