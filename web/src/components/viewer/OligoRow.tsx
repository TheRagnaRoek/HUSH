import { CopyButton } from "../common/CopyButton";

export function OligoRow({ label, seq }: { label: string; seq: string }) {
  return (
    <div className="oligo">
      <span className="ol">{label}</span>
      <span className="os mono">{seq}</span>
      <CopyButton text={seq} />
    </div>
  );
}
