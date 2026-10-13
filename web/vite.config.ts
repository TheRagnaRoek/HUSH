import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// Standalone SPA served at the site root; it calls the HUSH API at VITE_API_URL
// (default http://localhost:8000), which already sends permissive CORS, so no
// dev proxy is needed. Unlike the wiki build, there is no team-slug base path.
export default defineConfig({
  base: "/",
  plugins: [react()],
});
