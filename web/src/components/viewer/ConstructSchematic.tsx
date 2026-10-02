export function ConstructSchematic({
  sense,
  loop,
  antisense,
  terminator,
}: {
  sense: string;
  loop: string;
  antisense: string;
  terminator: string;
}) {
  return (
    <div className="schematic">
      <div className="segbar">
        <div className="seg" style={{ background: "var(--seg-sense)" }}>
          <span className="sl">sense</span>
          {sense}
        </div>
        <div className="seg" style={{ background: "var(--seg-loop)" }}>
          <span className="sl">loop</span>
          {loop}
        </div>
        <div className="seg" style={{ background: "var(--seg-anti)" }}>
          <span className="sl">antisense (guide)</span>
          {antisense}
        </div>
        <div className="seg" style={{ background: "var(--seg-term)" }}>
          <span className="sl">Pol III term</span>
          {terminator}
        </div>
      </div>
    </div>
  );
}
