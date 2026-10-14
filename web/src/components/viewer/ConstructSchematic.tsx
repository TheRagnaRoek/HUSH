export function ConstructSchematic({
  sense,
  senseOverhang = "",
  loop,
  antisense,
  guideOverhang = "",
  terminator,
}: {
  sense: string;
  senseOverhang?: string;
  loop: string;
  antisense: string;
  guideOverhang?: string;
  terminator: string;
}) {
  return (
    <div className="schematic">
      <div className="segbar">
        <div className="seg" style={{ background: "var(--seg-sense)" }}>
          <span className="sl">sense</span>
          {sense}
        </div>
        {senseOverhang && (
          <div className="seg" style={{ background: "var(--seg-oh)" }}>
            <span className="sl">3&prime; overhang</span>
            {senseOverhang}
          </div>
        )}
        <div className="seg" style={{ background: "var(--seg-loop)" }}>
          <span className="sl">loop</span>
          {loop}
        </div>
        <div className="seg" style={{ background: "var(--seg-anti)" }}>
          <span className="sl">antisense (guide)</span>
          {antisense}
        </div>
        {guideOverhang && (
          <div className="seg" style={{ background: "var(--seg-oh)" }}>
            <span className="sl">3&prime; overhang</span>
            {guideOverhang}
          </div>
        )}
        {terminator && (
          <div className="seg" style={{ background: "var(--seg-term)" }}>
            <span className="sl">Pol III term</span>
            {terminator}
          </div>
        )}
      </div>
    </div>
  );
}
