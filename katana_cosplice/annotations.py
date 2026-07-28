from __future__ import annotations

from dataclasses import dataclass, field

from katana_cosplice.architecture import AnnotationProvider, AnnotationProviderSpec


ANNOTATION_SPECS: tuple[AnnotationProviderSpec, ...] = (
    AnnotationProviderSpec(
        key="deg",
        title="DEG support",
        scope="gene",
        required_inputs=("standardized DEG table",),
        description="Annotate modules with expression evidence; never redefine AS-event membership.",
    ),
    AnnotationProviderSpec(
        key="dtu",
        title="DTU / isoform support",
        scope="gene and transcript",
        required_inputs=("quant.sf or DTU results", "GTF"),
        description="Annotate module genes with isoform-usage support.",
    ),
    AnnotationProviderSpec(
        key="rbp",
        title="RBP support",
        scope="gene and module",
        required_inputs=("RBP reference list",),
        description="Prioritize plausible RBP-linked modules without claiming causality.",
    ),
    AnnotationProviderSpec(
        key="motif",
        title="Splicing motif support",
        scope="event and module",
        required_inputs=("genome FASTA", "GTF"),
        description="Test sequence-motif enrichment around event regions.",
    ),
    AnnotationProviderSpec(
        key="go_kegg",
        title="GO / KEGG enrichment",
        scope="module gene set",
        required_inputs=("gene identifiers", "organism annotation source"),
        description="Describe biological functions represented by module genes.",
    ),
)


@dataclass
class AnnotationRegistry:
    """Provider registry kept separate from core network discovery."""

    providers: dict[str, AnnotationProvider] = field(default_factory=dict)

    def register(self, provider: AnnotationProvider) -> None:
        key = provider.spec.key
        if key in self.providers:
            raise ValueError(f"Annotation provider already registered: {key}")
        self.providers[key] = provider

    def get(self, key: str) -> AnnotationProvider:
        try:
            return self.providers[key]
        except KeyError as exc:
            raise KeyError(f"Annotation provider is not configured: {key}") from exc

    def status(self) -> list[dict[str, object]]:
        return [
            {
                "key": spec.key,
                "title": spec.title,
                "scope": spec.scope,
                "configured": spec.key in self.providers,
                "required_inputs": list(spec.required_inputs),
                "description": spec.description,
            }
            for spec in ANNOTATION_SPECS
        ]
