import { useState } from "react";
import { Check, ChevronDown, ChevronUp, CircleX, Clock3, DatabaseZap, LoaderCircle } from "lucide-react";
import type { Job } from "../types";
import { useMutation, useQueryClient } from "@tanstack/react-query";
import { api } from "../api";
import { displayModuleName } from "../branding";

type Props = { jobs: Job[] };

function JobIcon({ job }: { job: Job }) {
  if (job.status === "running" || job.status === "queued") return <LoaderCircle className="spin" size={16} />;
  if (job.status === "failed" || job.status === "cancelled") return <CircleX size={16} />;
  if (job.status === "cache_hit") return <DatabaseZap size={16} />;
  return <Check size={16} />;
}

export function JobDrawer({ jobs }: Props) {
  const client = useQueryClient();
  const retry = useMutation({mutationFn: api.retryJob, onSuccess: () => {void client.invalidateQueries({queryKey: ["jobs"]});}});
  const [expanded, setExpanded] = useState(false);
  const active = jobs.find((job) => job.status === "queued" || job.status === "running");
  const latest = active ?? jobs[0];
  return (
    <section className={expanded ? "job-drawer expanded" : "job-drawer"}>
      <button type="button" className="job-summary" onClick={() => setExpanded(!expanded)}>
        <span className="job-summary-title"><Clock3 size={16} />Activity</span>
        {latest ? (
          <span className={`job-state ${latest.status}`}>
            <JobIcon job={latest} />
            {latest.cache_status === "hit" ? "Cache hit" : displayModuleName(latest.message)}
          </span>
        ) : <span className="job-state idle">No analysis has run</span>}
        {expanded ? <ChevronDown size={16} /> : <ChevronUp size={16} />}
      </button>
      {expanded && (
        <div className="job-list">
          {retry.error && <p role="alert">{retry.error.message}</p>}
          {jobs.length === 0 ? <p>No recorded jobs. Opening a project does not start an analysis.</p> : jobs.map((job) => (
            <article key={job.job_id}>
              <JobIcon job={job} />
              <div>
                <strong>{job.analysis_id}</strong>
                <span>{job.cache_status === "hit" ? "Compatible cache reused" : displayModuleName(job.message)}</span>
                {job.error && <small>{job.error}</small>}
                {["failed", "cancelled"].includes(job.status) && <button className="secondary-action" disabled={retry.isPending}
                  title="Reuse existing stage checkpoints where supported; otherwise restart this workflow. Inputs and method version must be unchanged."
                  onClick={() => retry.mutate(job.job_id)}>Retry / resume saved run</button>}
              </div>
              <div className="job-progress" aria-label={`${Math.round(job.progress * 100)} percent`}>
                <span style={{ width: `${job.progress * 100}%` }} />
              </div>
            </article>
          ))}
        </div>
      )}
    </section>
  );
}
