import { LOOPS, loopInfo, MIN_TERMINATOR_LENGTH } from "../../lib/loops";
import { useSession } from "../../session/useSession";

/** Binds the shRNA scaffold options (loop, poly-U terminator) to the session's
 * ConvertOptions, which runConvert sends to the API. */
export function ScaffoldOptions() {
  const { convertOptions, updateConvertOptions } = useSession();
  const loop = convertOptions.loop ?? "jensen";
  const include = convertOptions.include_terminator ?? true;
  const termLen = convertOptions.terminator_length ?? 6;
  const info = loopInfo(loop);

  return (
    <div className="scaffold-opts">
      <label className="so-field">
        <span className="so-label">Loop</span>
        <select value={loop} onChange={(e) => updateConvertOptions({ loop: e.target.value })}>
          {LOOPS.map((l) => (
            <option key={l.name} value={l.name}>
              {l.label}
            </option>
          ))}
        </select>
      </label>

      <label className="so-field so-check">
        <input
          type="checkbox"
          checked={include}
          onChange={(e) => updateConvertOptions({ include_terminator: e.target.checked })}
        />
        <span className="so-label">Pol III poly-U terminator</span>
      </label>

      <label className="so-field">
        <span className="so-label">length</span>
        <input
          type="number"
          min={MIN_TERMINATOR_LENGTH}
          value={termLen}
          disabled={!include}
          onChange={(e) => {
            const n = e.target.valueAsNumber;
            if (!Number.isNaN(n)) updateConvertOptions({ terminator_length: Math.round(n) });
          }}
        />
      </label>

      {info && <span className="so-cite">{info.citation}</span>}
    </div>
  );
}
