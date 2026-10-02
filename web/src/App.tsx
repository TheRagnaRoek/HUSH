// Bootstrap first so our tokens/theme/app styles override it.
import "bootstrap/dist/css/bootstrap.min.css";
import "./styles/fonts.css";
import "./styles/tokens.css";
import "./styles/theme.css";
import "./styles/app.css";

import { Route, Routes } from "react-router-dom";
import { AppShell } from "./components/layout/AppShell";
import { DesignerPage } from "./DesignerPage";
import { SessionProvider } from "./session/SessionProvider";

export default function App() {
  return (
    <SessionProvider>
      <AppShell>
        <Routes>
          <Route path="/" element={<DesignerPage />} />
          <Route path="*" element={<DesignerPage />} />
        </Routes>
      </AppShell>
    </SessionProvider>
  );
}
