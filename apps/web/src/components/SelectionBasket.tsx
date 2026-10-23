import { BookmarkPlus, X } from "lucide-react";

type Props = {
  count: number;
  saving: boolean;
  onSave: () => void;
  onClear: () => void;
  onAddGenes?: () => void;
  onAddEvents?: () => void;
};

export function SelectionBasket({ count, saving, onSave, onClear, onAddGenes, onAddEvents }: Props) {
  if (!count) return null;
  return (
    <div className="selection-basket" role="status">
      <div>
        <BookmarkPlus size={18} />
        <span><strong>{count}</strong> source row{count === 1 ? "" : "s"} in selection</span>
      </div>
      <div>
        <button type="button" className="ghost-action" onClick={onClear}><X size={15} />Clear</button>
        {onAddGenes && <button type="button" className="ghost-action" onClick={onAddGenes}>Add genes to targets</button>}
        {onAddEvents && <button type="button" className="ghost-action" onClick={onAddEvents}>Add events to targets</button>}
        <button type="button" className="primary-action compact" disabled={saving} onClick={onSave}>
          {saving ? "Saving..." : "Save selection"}
        </button>
      </div>
    </div>
  );
}
