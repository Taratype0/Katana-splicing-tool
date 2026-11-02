"""Explicit engine selection; no fallback, model filtering or data migration."""
ENGINES = ('eventgraph_20261030', 'legacy')


def analyze_consequences(*args, analysis_engine='eventgraph_20261030', **kwargs):
    if analysis_engine == 'legacy':
        from katana_next.science.motiftrace_execution import analyze_consequences as run
    elif analysis_engine == 'eventgraph_20261030':
        from katana_next.science.motiftrace_eventgraph_execution import analyze_consequences as run
    else:
        raise ValueError('Unknown MotifTrace analysis_engine: ' + str(analysis_engine))
    return run(*args, **kwargs)
