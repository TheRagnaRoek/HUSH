/**
 * Minimal typed fetch client for the HUSH API. Thin on purpose: the backend
 * already handles CORS, errors come back as { detail }, and the app's wrappers
 * (hush.ts) own the method/path/type mapping. Swap in openapi-fetch later if a
 * generated client is preferred - `npm run gen:api` emits the schema for it.
 */

export const API_URL: string =
  (import.meta.env.VITE_API_URL as string | undefined)?.replace(/\/$/, "") ??
  "http://localhost:8000";

export class ApiError extends Error {
  status: number;
  constructor(status: number, message: string) {
    super(message);
    this.name = "ApiError";
    this.status = status;
  }
}

async function parseError(res: Response): Promise<never> {
  let detail = `${res.status} ${res.statusText}`;
  try {
    const body = (await res.json()) as { detail?: string };
    if (body?.detail) detail = body.detail;
  } catch {
    // non-JSON error body; keep the status line
  }
  throw new ApiError(res.status, detail);
}

export async function apiGet<T>(path: string): Promise<T> {
  const res = await fetch(`${API_URL}${path}`, {
    headers: { accept: "application/json" },
  });
  if (!res.ok) return parseError(res);
  return (await res.json()) as T;
}

export async function apiPost<T>(path: string, body?: unknown): Promise<T> {
  const res = await fetch(`${API_URL}${path}`, {
    method: "POST",
    headers: { "content-type": "application/json", accept: "application/json" },
    body: body === undefined ? undefined : JSON.stringify(body),
  });
  if (!res.ok) return parseError(res);
  return (await res.json()) as T;
}
