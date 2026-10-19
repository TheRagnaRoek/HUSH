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
        <div className="seg sense">
          <span className="sl">Sense (Passenger)</span>
          {sense}
        </div>
        {senseOverhang && (
          <div className="seg oh">
            <span className="sl">3&prime; Overhang</span>
            {senseOverhang}
          </div>
        )}
        <div className="seg loop">
          <span className="sl">Loop</span>
          {loop}
        </div>
        <div className="seg anti">
          <span className="sl">Antisense (Guide)</span>
          {antisense}
        </div>
        {guideOverhang && (
          <div className="seg oh">
            <span className="sl">3&prime; Overhang</span>
            {guideOverhang}
          </div>
        )}
        {terminator && (
          <div className="seg term">
            <span className="sl">Pol III Terminator</span>
            {terminator}
          </div>
        )}
      </div>
    </div>
  );
}
