import { BookmarkPlus, X } from "lucide-react";

type Props = {
  count: number;
  saving: boolean;
  onSave: () => void;
  onClear: () => void;
};

export function SelectionBasket({ count, saving, onSave, onClear }: Props) {
  if (!count) return null;
  return (
    <div className="selection-basket" role="status">
      <div>
        <BookmarkPlus size={18} />
        <span><strong>{count}</strong> event{count === 1 ? "" : "s"} in selection</span>
      </div>
      <div>
        <button type="button" className="ghost-action" onClick={onClear}><X size={15} />Clear</button>
        <button type="button" className="primary-action compact" disabled={saving} onClick={onSave}>
          {saving ? "Saving..." : "Save selection"}
        </button>
      </div>
    </div>
  );
}
