import { useEffect, useRef, useState } from "react";
import { ArrowUp, File, Folder, FolderOpen, X } from "lucide-react";
import { api } from "../api";

export type LocalDirectory = {path:string; parent:string|null; entries:{name:string;path:string;directory:boolean}[]; locations:{name:string;path:string}[]};

export function PathPicker({label, value, onChange, kind = "either", extensions = [], placeholder = "Choose a file or folder"}:
  {label:string;value:string;onChange:(path:string)=>void;kind?:"file"|"directory"|"either";extensions?:string[];placeholder?:string}) {
  const [open, setOpen] = useState(false);
  const [directory, setDirectory] = useState<LocalDirectory>();
  const [selected, setSelected] = useState("");
  const [search, setSearch] = useState("");
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);
  const dialog = useRef<HTMLDialogElement>(null);
  const serial = useRef(0);
  const trigger = useRef<HTMLButtonElement>(null);

  async function browse(path:string) {
    const request = ++serial.current;
    setBusy(true); setError(""); setSelected(""); setSearch("");
    try {const result = await api.browseLocalFiles(path); if(request === serial.current) {
      setDirectory(result);
      try {localStorage.setItem("katana:last-browse-folder", result.path);} catch { /* Browsing works without persistent storage. */ }
    }}
    catch(error) {if(request === serial.current) setError(error instanceof Error ? error.message : String(error));}
    finally {if(request === serial.current) setBusy(false);}
  }
  function close() {++serial.current; setOpen(false); trigger.current?.focus();}
  useEffect(()=>{if(open) dialog.current?.showModal();}, [open]);
  const choice = selected || (kind !== "file" ? directory?.path : "");
  const entries = directory?.entries.filter(entry =>
    (entry.directory || (kind !== "directory" && (!extensions.length || extensions.some(ext=>entry.name.toLowerCase().endsWith(ext))))) &&
    entry.name.toLowerCase().includes(search.toLowerCase())) || [];

  return <div className="path-picker">
    <div className="path-picker-field"><input aria-label={label} value={value} onChange={event=>onChange(event.target.value)} placeholder={placeholder}/>
      <button ref={trigger} type="button" aria-label={`Browse ${label}`} onClick={()=>{setOpen(true); let initial = value;
        if(!initial) {try {initial = localStorage.getItem("katana:last-browse-folder") || "";} catch { /* Use Home. */ }}
        void browse(initial);}}><FolderOpen size={16}/>Browse</button></div>
    {open && <dialog ref={dialog} className="file-dialog" aria-label={`Choose ${label}`} onCancel={event=>{event.preventDefault();close();}}>
      <header><h2>{kind === "directory" ? "Choose folder" : kind === "file" ? "Choose file" : "Choose file or folder"}</h2><button type="button" aria-label="Close file browser" onClick={close}><X size={20}/></button></header>
      <div className="file-browser"><aside>{directory?.locations.map(location=><button type="button" key={location.path} onClick={()=>void browse(location.path)}><Folder size={16}/>{location.name}</button>)}</aside>
        <div className="file-browser-main"><div className="file-browser-location"><button type="button" aria-label="Parent folder" disabled={busy || !directory?.parent} onClick={()=>void browse(directory!.parent!)}><ArrowUp size={16}/></button><span title={directory?.path}>{directory?.path || "Loading…"}</span></div>
          <input aria-label="Filter files and folders" placeholder="Find in this folder" value={search} onChange={event=>setSearch(event.target.value)}/>
          {error && <p className="form-error" role="alert">{error} <button type="button" onClick={()=>void browse("")}>Go to Home</button></p>}
          <div className="file-browser-list" aria-busy={busy}>{busy ? <p>Loading…</p> : entries.map(entry=><button type="button" className={selected === entry.path ? "selected" : ""} key={entry.path}
            onClick={()=>entry.directory ? void browse(entry.path) : setSelected(entry.path)}>
            {entry.directory ? <Folder size={18}/> : <File size={18}/>}<span>{entry.name}</span></button>)}
            {!busy && !error && !entries.length && <p>No matching files or folders</p>}</div>
        </div></div>
      <footer><span title={choice}>{selected.split(/[\\/]/).pop() || (kind !== "file" ? "Current folder" : "Select a file")}</span>
        <button type="button" onClick={close}>Cancel</button><button type="button" className="primary-action" disabled={busy || !!error || !choice} onClick={()=>{onChange(choice!);close();}}>Choose{selected || kind === "file" ? " file" : " folder"}</button></footer>
    </dialog>}
  </div>;
}
