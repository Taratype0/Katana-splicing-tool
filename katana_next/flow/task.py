"""Stable task entry point used by Python and Nextflow."""
from .runtime import execute_task, validate_receipt

__all__ = ["execute_task", "validate_receipt"]
