import type { ReactNode } from "react";
import { GuideBar } from "./GuideBar";
import { Header } from "./Header";
import { WorkflowRail } from "./WorkflowRail";

export function AppShell({ children }: { children: ReactNode }) {
  return (
    <>
      <div className="brandbar" />
      <div className="app">
        <Header />
        <WorkflowRail />
        <GuideBar />
        {children}
        <div className="foot">
          HUSH Designer &middot; Project HIDE &middot; allele-aware siRNA/shRNA design
        </div>
      </div>
    </>
  );
}
