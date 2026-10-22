from __future__ import annotations

import threading
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Callable

from katana_next.domain.models import JobEvent, JobRecord, JobStatus


class JobCancelled(RuntimeError):
    pass


@dataclass
class JobContext:
    job_id: str
    _manager: "JobManager"
    _cancel_event: threading.Event

    def report(self, progress: float, message: str) -> None:
        self._manager._report(self.job_id, progress, message)

    def check_cancelled(self) -> None:
        if self._cancel_event.is_set():
            raise JobCancelled("Job cancelled by user.")


JobRunner = Callable[[JobContext], list[str]]


class JobManager:
    def __init__(self, max_workers: int = 2) -> None:
        self._executor = ThreadPoolExecutor(
            max_workers=max_workers, thread_name_prefix="katana-next"
        )
        self._jobs: dict[str, JobRecord] = {}
        self._futures: dict[str, Future[list[str]]] = {}
        self._cancel_events: dict[str, threading.Event] = {}
        self._lock = threading.RLock()
        self._journals = {}

    def bind_project(self, project_id, output_root, *, execute=False):
        from katana_next.infrastructure.job_journal import JobJournal
        with self._lock:
            if project_id not in self._journals:
                self._journals[project_id] = JobJournal(output_root)
            journal = self._journals[project_id]
            if journal.directory.parent.resolve() != output_root.resolve():
                raise ValueError("The same project ID cannot use two output directories in one runtime")
            if execute:
                journal.acquire_executor()
            for record in journal.read():
                if record.job_id not in self._futures:
                    self._jobs[record.job_id] = record

    def _persist(self, job_id):
        record = self._jobs[job_id]
        journal = self._journals.get(record.project_id)
        if journal is not None:
            journal.save(record)

    def submit(
        self,
        *,
        project_id: str,
        analysis_id: str,
        run_fingerprint: str,
        runner: JobRunner,
        parameters: dict | None = None,
    ) -> JobRecord:
        job = JobRecord(
            project_id=project_id,
            analysis_id=analysis_id,
            status=JobStatus.QUEUED,
            run_fingerprint=run_fingerprint,
            parameters=parameters or {},
            message="Queued for rebuild",
            cache_status="rebuild",
            events=(JobEvent(sequence=1, progress=0.0, message="Queued for rebuild"),),
        )
        cancel_event = threading.Event()
        with self._lock:
            self._jobs[job.job_id] = job
            self._persist(job.job_id)
            self._cancel_events[job.job_id] = cancel_event
            future = self._executor.submit(
                self._execute, job.job_id, runner, cancel_event
            )
            self._futures[job.job_id] = future
        return job

    def cache_hit(
        self,
        *,
        project_id: str,
        analysis_id: str,
        run_fingerprint: str,
        artifact_ids: list[str],
        parameters: dict | None = None,
    ) -> JobRecord:
        now = datetime.now(timezone.utc)
        event = JobEvent(sequence=1, progress=1.0, message="Compatible cache hit")
        job = JobRecord(
            project_id=project_id,
            analysis_id=analysis_id,
            status=JobStatus.CACHE_HIT,
            run_fingerprint=run_fingerprint,
            parameters=parameters or {},
            progress=1.0,
            message=event.message,
            cache_status="hit",
            artifact_ids=tuple(artifact_ids),
            events=(event,),
            started_at=now,
            completed_at=now,
        )
        with self._lock:
            self._jobs[job.job_id] = job
            self._persist(job.job_id)
        return job

    def _execute(
        self, job_id: str, runner: JobRunner, cancel_event: threading.Event
    ) -> list[str]:
        self._transition(job_id, JobStatus.RUNNING, "Rebuilding artifact", started=True)
        context = JobContext(job_id=job_id, _manager=self, _cancel_event=cancel_event)
        try:
            artifact_ids = runner(context)
            context.check_cancelled()
        except JobCancelled as exc:
            self._transition(job_id, JobStatus.CANCELLED, str(exc), completed=True)
            return []
        except Exception as exc:
            self._transition(
                job_id,
                JobStatus.FAILED,
                "Analysis failed",
                completed=True,
                error=str(exc),
            )
            return []
        self._complete(job_id, artifact_ids)
        return artifact_ids

    def _complete(self, job_id: str, artifact_ids: list[str]) -> None:
        with self._lock:
            current = self._jobs[job_id]
            event = JobEvent(
                sequence=len(current.events) + 1,
                progress=1.0,
                message="Rebuild complete",
            )
            self._jobs[job_id] = current.model_copy(
                update={
                    "status": JobStatus.COMPLETE,
                    "progress": 1.0,
                    "message": event.message,
                    "artifact_ids": tuple(artifact_ids),
                    "events": (*current.events, event),
                    "completed_at": datetime.now(timezone.utc),
                }
            )
            self._persist(job_id)

    def _report(self, job_id: str, progress: float, message: str) -> None:
        with self._lock:
            current = self._jobs[job_id]
            event = JobEvent(
                sequence=len(current.events) + 1,
                progress=max(0.0, min(1.0, progress)),
                message=message,
            )
            self._jobs[job_id] = current.model_copy(
                update={
                    "progress": event.progress,
                    "message": message,
                    "events": (*current.events, event),
                }
            )
            self._persist(job_id)

    def _transition(
        self,
        job_id: str,
        status: JobStatus,
        message: str,
        *,
        started: bool = False,
        completed: bool = False,
        error: str | None = None,
    ) -> None:
        with self._lock:
            current = self._jobs[job_id]
            event = JobEvent(
                sequence=len(current.events) + 1,
                progress=current.progress,
                message=message,
            )
            updates: dict[str, object] = {
                "status": status,
                "message": message,
                "events": (*current.events, event),
                "error": error,
            }
            now = datetime.now(timezone.utc)
            if started:
                updates["started_at"] = now
            if completed:
                updates["completed_at"] = now
            self._jobs[job_id] = current.model_copy(update=updates)
            self._persist(job_id)

    def get(self, job_id: str) -> JobRecord | None:
        with self._lock:
            if job_id not in self._futures:
                record = self._jobs.get(job_id)
                journal = self._journals.get(record.project_id) if record is not None else None
                if journal is not None:
                    for saved in journal.read(job_id):
                        self._jobs[job_id] = saved
            return self._jobs.get(job_id)

    def list(self, project_id: str | None = None) -> list[JobRecord]:
        with self._lock:
            for identity, journal in self._journals.items():
                if project_id is None or identity == project_id:
                    for record in journal.read():
                        if record.job_id not in self._futures:
                            self._jobs[record.job_id] = record
            jobs = list(self._jobs.values())
        if project_id:
            jobs = [job for job in jobs if job.project_id == project_id]
        return sorted(jobs, key=lambda item: item.created_at, reverse=True)

    def cancel(self, job_id: str) -> JobRecord:
        with self._lock:
            current = self.get(job_id)
            if current is None:
                raise KeyError(job_id)
            if current.status not in {JobStatus.QUEUED, JobStatus.RUNNING}:
                return current
            event = self._cancel_events.get(job_id)
            if event is None:
                raise ValueError("This job belongs to another executor. Cancel it from the runtime that started it.")
            event.set()
            return self._jobs[job_id]

    def wait(self, job_id: str, timeout: float = 30.0) -> JobRecord:
        with self._lock:
            future = self._futures.get(job_id)
        if future is not None:
            future.result(timeout=timeout)
        job = self.get(job_id)
        if job is None:
            raise KeyError(job_id)
        return job

    def shutdown(self) -> None:
        with self._lock:
            for event in self._cancel_events.values():
                event.set()
        self._executor.shutdown(wait=True, cancel_futures=True)
        with self._lock:
            for identity, future in self._futures.items():
                if future.cancelled():
                    self._transition(identity, JobStatus.CANCELLED, "Cancelled before execution during shutdown", completed=True)
            for journal in self._journals.values():
                journal.close()
