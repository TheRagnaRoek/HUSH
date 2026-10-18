import { Link } from "react-router-dom";
import type { RefId } from "../../reference/content";
import { shortCites, type SourceKey } from "../../reference/sources";

/** "Sources: ..." in brief, then a link to the matching reference section. */
export function RefNote({ cites, to }: { cites: SourceKey[]; to: RefId }) {
  return (
    <>
      {cites.length > 0 && <span className="dtcite">Sources: {shortCites(cites)}</span>}
      <Link className="reflink" to={{ pathname: "/reference", hash: `#${to}` }}>
        Read more in References &rarr;
      </Link>
    </>
  );
}
