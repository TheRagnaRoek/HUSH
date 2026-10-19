import type { ReactNode } from "react";
import { useLocation } from "react-router-dom";
import { useSession } from "../../session/useSession";
import { GuideBar } from "./GuideBar";
import { Header } from "./Header";
import { WorkflowRail } from "./WorkflowRail";

export function AppShell({ children }: { children: ReactNode }) {
  const { expert } = useSession();
  // The reference page sits outside the design process, so it has no steps or guide.
  const onReference = useLocation().pathname.startsWith("/reference");
  return (
    <>
      <div className="brandbar" />
      <div className={expert && !onReference ? "app wide" : "app"}>
        <Header />
        {!expert && !onReference && <WorkflowRail />}
        {!expert && !onReference && <GuideBar />}
        {children}
        <div className="foot">
          HUSH Designer &middot; Team HIDE &middot; allele-aware siRNA/shRNA design
        </div>
      </div>
    </>
  );
}
