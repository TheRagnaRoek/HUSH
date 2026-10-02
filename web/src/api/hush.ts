/**
 * Typed wrappers over the HUSH API endpoints. The session flow mirrors
 * service.HushSession: POST the inputs once to /sessions, then design/convert by
 * id so the server reuses its cached panel alignment across re-designs.
 */

import { apiGet, apiPost } from "./client";
import type {
  ConvertOptions,
  ConvertResponse,
  CreateSessionRequest,
  DesignParams,
  DesignResponse,
  HealthResponse,
  InspectRequest,
  InspectResponse,
  SessionResponse,
} from "./types";

export function health(): Promise<HealthResponse> {
  return apiGet<HealthResponse>("/health");
}

export function inspect(req: InspectRequest): Promise<InspectResponse> {
  return apiPost<InspectResponse>("/inspect", req);
}

export function createSession(req: CreateSessionRequest): Promise<SessionResponse> {
  return apiPost<SessionResponse>("/sessions", req);
}

export function getSession(id: string): Promise<SessionResponse> {
  return apiGet<SessionResponse>(`/sessions/${encodeURIComponent(id)}`);
}

export function sessionDesign(id: string, params?: DesignParams): Promise<DesignResponse> {
  return apiPost<DesignResponse>(`/sessions/${encodeURIComponent(id)}/design`, params ?? {});
}

export function sessionConvert(id: string, options?: ConvertOptions): Promise<ConvertResponse> {
  return apiPost<ConvertResponse>(`/sessions/${encodeURIComponent(id)}/convert`, options ?? {});
}

/** Stateless one-shot design (throwaway session server-side). */
export function design(
  target_fasta: string,
  panel_fasta: string | null,
  params?: DesignParams,
): Promise<DesignResponse> {
  return apiPost<DesignResponse>("/design", { target_fasta, panel_fasta, params });
}
