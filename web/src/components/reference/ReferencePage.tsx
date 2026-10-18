import { useEffect, useRef, useState } from "react";
import { Link, useLocation } from "react-router-dom";
import { GROUPS, SECTIONS, type RefSection } from "../../reference/content";
import { SOURCES } from "../../reference/sources";

function Section({ s, active }: { s: RefSection; active: boolean }) {
  return (
    <section id={s.id} className={active ? "ref-sec active" : "ref-sec"}>
      <h3>{s.title}</h3>
      <p className="ref-sum">{s.summary}</p>

      {s.table && (
        <div className="ref-tablewrap">
        <table className="ref-table">
          <thead>
            <tr>
              {s.table.head.map((h) => (
                <th key={h}>{h}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {s.table.rows.map((row) => (
              <tr key={row.join("|")}>
                {row.map((cell, i) => (
                  <td key={i}>{cell}</td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
        </div>
      )}

      <ul className="ref-points">
        {s.points.map((p) => (
          <li key={p}>{p}</li>
        ))}
      </ul>

      {s.legend && (
        <dl className="ref-legend">
          {s.legend.map(([term, text]) => (
            <div key={term}>
              <dt>{term}</dt>
              <dd>{text}</dd>
            </div>
          ))}
        </dl>
      )}

      {s.defaults && (
        <p className="ref-def">
          <b>Default.</b> {s.defaults}
        </p>
      )}
      {s.cites.length > 0 && (
        <p className="ref-src">
          <b>Sources.</b> {s.cites.map((k) => SOURCES[k].full).join("; ")}
        </p>
      )}
    </section>
  );
}

// A section counts as being read once its top has passed this line (below the sticky bar).
const READ_LINE = 96;
const SOURCES_ID = "sources";

/** The in-app reference: design rules, why they exist, and what they rest on. */
export function ReferencePage() {
  const { hash } = useLocation();
  const [current, setCurrent] = useState<string>(SECTIONS[0].id);
  const tocRef = useRef<HTMLElement>(null);

  // The router does not scroll to a hash by itself.
  useEffect(() => {
    const id = hash.replace(/^#/, "");
    if (!id) {
      window.scrollTo({ top: 0, behavior: "instant" });
      return;
    }
    document.getElementById(id)?.scrollIntoView({ block: "start", behavior: "instant" });
  }, [hash]);

  // Follow the reading position. The last sections cannot scroll up to the line, so at
  // the bottom of the page a linked section that is on screen counts as the one read.
  useEffect(() => {
    const ids = [...SECTIONS.map((x) => x.id), SOURCES_ID];
    const linked = hash.replace(/^#/, "");

    const update = () => {
      let index = 0;
      ids.forEach((id, i) => {
        const el = document.getElementById(id);
        if (el && el.getBoundingClientRect().top <= READ_LINE) index = i;
      });
      const linkedIndex = ids.indexOf(linked);
      const atBottom =
        window.innerHeight + window.scrollY >= document.documentElement.scrollHeight - 4;
      if (atBottom && linkedIndex > index) {
        const el = document.getElementById(linked);
        if (el && el.getBoundingClientRect().top < window.innerHeight) index = linkedIndex;
      }
      setCurrent(ids[index]);
    };
    // Cheap enough to run on every scroll event; React skips an unchanged state.
    update();
    window.addEventListener("scroll", update, { passive: true });
    window.addEventListener("resize", update);
    return () => {
      window.removeEventListener("scroll", update);
      window.removeEventListener("resize", update);
    };
  }, [hash]);

  // Keep the current entry inside the contents list, which scrolls on its own.
  useEffect(() => {
    const toc = tocRef.current;
    const link = toc?.querySelector<HTMLElement>('a[aria-current="true"]');
    if (!toc || !link) return;
    if (link.offsetTop < toc.scrollTop) {
      toc.scrollTop = link.offsetTop - 8;
    } else if (link.offsetTop + link.offsetHeight > toc.scrollTop + toc.clientHeight) {
      toc.scrollTop = link.offsetTop + link.offsetHeight - toc.clientHeight + 8;
    }
  }, [current]);

  const sources = Object.values(SOURCES)
    .map((x) => x.full)
    .sort((a, b) => a.localeCompare(b));
  const currentTitle =
    current === SOURCES_ID ? "All sources" : SECTIONS.find((x) => x.id === current)?.title;

  return (
    <div className="workspace ref">
      <div className="ref-bar">
        <Link to="/" className="hbtn ghost">
          &larr; Back to the designer
        </Link>
        {currentTitle && <span className="ref-here">{currentTitle}</span>}
      </div>

      <div className="whead">
        <h2>Design rules and sources</h2>
        <p>
          A brief overview of what each rule tests, why it exists and what it rests on. It is
          condensed from the HUSH documentation, which is the full reference, and is not part of
          the design process.
        </p>
      </div>

      <div className="ref-layout">
        <nav className="ref-toc" aria-label="Sections" ref={tocRef}>
          {GROUPS.map((g) => (
            <div key={g}>
              <div className="funnel-h">{g}</div>
              <ul>
                {SECTIONS.filter((x) => x.group === g).map((x) => (
                  <li key={x.id}>
                    <Link
                      to={{ pathname: "/reference", hash: `#${x.id}` }}
                      aria-current={current === x.id ? "true" : undefined}
                    >
                      {x.title}
                    </Link>
                  </li>
                ))}
              </ul>
            </div>
          ))}
          <ul>
            <li>
              <Link
                to={{ pathname: "/reference", hash: `#${SOURCES_ID}` }}
                aria-current={current === SOURCES_ID ? "true" : undefined}
              >
                All sources
              </Link>
            </li>
          </ul>
        </nav>

        <div className="ref-body">
          {SECTIONS.map((x) => (
            <Section key={x.id} s={x} active={x.id === current} />
          ))}

          <section id={SOURCES_ID} className={current === SOURCES_ID ? "ref-sec active" : "ref-sec"}>
            <h3>All sources</h3>
            <ul className="ref-points">
              {sources.map((x) => (
                <li key={x}>{x}</li>
              ))}
            </ul>
          </section>
        </div>
      </div>
    </div>
  );
}
