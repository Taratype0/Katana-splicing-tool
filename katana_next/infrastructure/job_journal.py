"""Project-local durable job history; never automatically re-execute saved jobs."""
import os
from pathlib import Path
import sqlite3
from contextlib import contextmanager

from katana_next.domain.models import JobRecord, JobStatus, JobEvent, utc_now


def process_alive(pid):
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.OpenProcess.restype = wintypes.HANDLE
        kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        kernel.GetExitCodeProcess.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
        kernel.CloseHandle.argtypes = [wintypes.HANDLE]
        handle = kernel.OpenProcess(0x1000, False, pid)
        if not handle:
            return ctypes.get_last_error() == 5
        try:
            code = ctypes.c_ulong()
            if not kernel.GetExitCodeProcess(handle, ctypes.byref(code)):
                return True
            return code.value == 259
        finally:
            kernel.CloseHandle(handle)
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


class JobJournal:
    def __init__(self, output_root):
        self.directory = Path(output_root)/".katana"
        self.directory.mkdir(parents=True, exist_ok=True)
        self.path = self.directory/"jobs.sqlite3"
        self.lease = None
        with self.connect() as db:
            db.execute("CREATE TABLE IF NOT EXISTS jobs (id TEXT PRIMARY KEY, pid INTEGER NOT NULL, record TEXT NOT NULL)")

    @contextmanager
    def connect(self):
        connection = sqlite3.connect(self.path, timeout=30)
        try:
            with connection:
                yield connection
        finally:
            connection.close()

    def acquire_executor(self):
        if self.lease is not None:
            return
        handle = (self.directory/"executor.lock").open("a+b")
        if handle.tell() == 0:
            handle.write(b"0"); handle.flush()
        handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            handle.close()
            raise ValueError("Another Katana executor owns this output directory. Close its runtime before starting jobs here; cached results remain readable.") from exc
        self.lease = handle

    def close(self):
        if self.lease is not None:
            self.lease.close()
            self.lease = None

    def save(self, record):
        with self.connect() as db:
            db.execute("INSERT INTO jobs VALUES (?,?,?) ON CONFLICT(id) DO UPDATE SET pid=excluded.pid, record=excluded.record",
                       (record.job_id, os.getpid(), record.model_dump_json()))

    def read(self, job_id=None):
        result = []
        with self.connect() as db:
            query = "SELECT id,pid,record FROM jobs"
            entries = db.execute(query if job_id is None else query+" WHERE id=?",
                                 () if job_id is None else (job_id,)).fetchall()
            for identity, pid, serialized in entries:
                record = JobRecord.model_validate_json(serialized)
                if record.status in {JobStatus.QUEUED, JobStatus.RUNNING} and not process_alive(pid):
                    message = "Interrupted by process exit; retry explicitly with unchanged inputs"
                    event = JobEvent(sequence=len(record.events)+1, progress=record.progress, message=message)
                    record = record.model_copy(update={"status": JobStatus.FAILED, "message": message,
                        "error": "Executor exited before completion", "completed_at": utc_now(), "events": (*record.events,event)})
                    db.execute("UPDATE jobs SET record=? WHERE id=? AND record=?", (record.model_dump_json(),identity,serialized))
                result.append(record)
        return result
