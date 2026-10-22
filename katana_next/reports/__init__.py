"""Read-only assembly of the approved seven-section Katana report."""


def build_report(project, output_dir, **options):
    """Build from completed project files, never starting upstream analyses."""
    from .engine import build_report as implementation
    return implementation(project, output_dir, **options)


__all__ = ['build_report']
