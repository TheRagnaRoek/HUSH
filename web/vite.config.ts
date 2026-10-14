import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// Standalone SPA served at the site root. The client calls its own origin, so
// dev proxies the API paths to the backend and a container serves both from one
// port - the same relative URLs either way, with no API host baked into the
// bundle. Unlike the wiki build, there is no team-slug base path.
const API_TARGET = process.env.HUSH_API_URL ?? "http://localhost:8000";
const API_PATHS = ["/health", "/inspect", "/design", "/sessions", "/openapi.json", "/docs"];

export default defineConfig({
  base: "/",
  plugins: [react()],
  server: {
    proxy: Object.fromEntries(
      API_PATHS.map((path) => [path, { target: API_TARGET, changeOrigin: true }]),
    ),
  },
});
