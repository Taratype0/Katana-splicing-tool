"""Release alignment and separately owned Python-backend index descriptors."""
from contextlib import contextmanager


def close_alignment(handle):
    # bamnostic's BGZF close() does not close Bai/Csi._io. On Windows that
    # descriptor keeps the index locked even after its BAM is closed.
    index = getattr(handle, "_index", None) if type(handle).__module__.startswith("bamnostic") else None
    index_stream = getattr(index, "_io", None)
    try:
        handle.close()
    finally:
        if index_stream is not None:
            index_stream.close()


@contextmanager
def closing_alignment(handle):
    try:
        yield handle
    finally:
        close_alignment(handle)
