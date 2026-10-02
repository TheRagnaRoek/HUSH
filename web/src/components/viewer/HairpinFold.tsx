import type { ConstructResult } from "../../api/types";

/**
 * Schematic stem-loop. TODO(viz): replace with a real fold rendered by fornac
 * (ViennaRNA) once ConvertResponse carries the dot-bracket structure; the engine
 * already computes it, it just isn't serialized yet.
 */
export function HairpinFold({ construct }: { construct: ConstructResult }) {
  const rungs = Array.from({ length: 9 }, (_, i) => 50 + i * 10);
  return (
    <svg width="140" height="150" viewBox="0 0 140 150" aria-label="stem-loop fold">
      <path
        d="M50 142 L50 42 Q50 20 70 20 Q90 20 90 42 L90 142"
        fill="none"
        stroke="var(--seg-anti)"
        strokeWidth="7"
        strokeLinecap="round"
      />
      <path
        d="M50 142 L50 42 Q50 20 70 20 Q90 20 90 42 L90 142"
        fill="none"
        stroke="var(--seg-sense)"
        strokeWidth="7"
        strokeLinecap="round"
        strokeDasharray="2 12"
      />
      {rungs.map((y, i) => (
        <line key={i} x1="53" y1={y} x2="87" y2={y} stroke="var(--grey)" strokeWidth="1.3" />
      ))}
      <text x="70" y="15" textAnchor="middle" fontSize="8" fill="var(--seg-loop)" fontFamily="DM Mono, monospace">
        {construct.loop_seq}
      </text>
      <line x1="90" y1="142" x2="108" y2="142" stroke="var(--seg-term)" strokeWidth="7" strokeLinecap="round" />
      <text x="116" y="145" fontSize="8" fill="var(--seg-term)" fontFamily="DM Mono, monospace">
        UUUU
      </text>
    </svg>
  );
}
