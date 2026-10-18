import { useCallback, useEffect, useMemo, useRef, useState, type ReactNode } from "react";
import {
  createSession,
  inspect as apiInspect,
  sessionConvert,
  sessionDesign,
} from "../api/hush";
import type {
  ConvertOptions,
  ConvertResponse,
  DesignParams,
  DesignResponse,
  InspectResponse,
} from "../api/types";
import { DEFAULT_LOOP, DEFAULT_TERMINATOR_LENGTH } from "../lib/loops";
import { DEFAULT_PARAMS } from "../lib/params";
import { constructForRegion, sampleConvert, sampleDesign, sampleEnrichment, sampleInspect } from "../lib/sample";
import { SessionContext, type SessionValue, type Step } from "./context";

/**
 * The working session. Opens on sample data; the run actions call the API when
 * FASTA input is present and fall back to the sample, flagged offline, on failure.
 */
export function SessionProvider({ children }: { children: ReactNode }) {
  const [step, setStep] = useState<Step>("candidates");
  const [expert, setExpert] = useState(false);

  const [targetFasta, setTargetFasta] = useState("");
  const [panelFasta, setPanelFasta] = useState("");
  const [params, setParams] = useState<DesignParams>(DEFAULT_PARAMS);
  const updateParams = useCallback(
    (patch: Partial<DesignParams>) => setParams((p) => ({ ...p, ...patch })),
    [],
  );
  const [convertOptions, setConvertOptions] = useState<ConvertOptions>({
    loop: DEFAULT_LOOP,
    stem_overhang: "keep",
    include_terminator: true,
    terminator_length: DEFAULT_TERMINATOR_LENGTH,
  });
  const updateConvertOptions = useCallback(
    (patch: Partial<ConvertOptions>) => setConvertOptions((o) => ({ ...o, ...patch })),
    [],
  );

  const [inspect, setInspect] = useState<InspectResponse | null>(sampleInspect);
  const [design, setDesign] = useState<DesignResponse | null>(sampleDesign);
  const [convert, setConvert] = useState<ConvertResponse | null>(sampleConvert);

  const [selectedRegion, setSelectedRegion] = useState<string | null>(
    sampleDesign.selected[0]?.region ?? null,
  );
  const [offline, setOffline] = useState(true);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const sessionId = useRef<string | null>(null);

  // Refresh the inspect summary (the targeting step reads panel members from
  // it) when the inputs change; debounced because textareas fire per keystroke.
  useEffect(() => {
    if (!targetFasta.trim()) {
      setInspect(sampleInspect);
      return;
    }
    let cancelled = false;
    const timer = setTimeout(() => {
      apiInspect({
        target_fasta: targetFasta,
        panel_fasta: panelFasta.trim() ? panelFasta : null,
      })
        .then((res) => {
          if (!cancelled) setInspect(res);
        })
        .catch(() => {
          // Show no panel rather than a stale one.
          if (!cancelled) setInspect(null);
        });
    }, 500);
    return () => {
      cancelled = true;
      clearTimeout(timer);
    };
  }, [targetFasta, panelFasta]);

  const ensureSession = useCallback(async (): Promise<string> => {
    if (sessionId.current) return sessionId.current;
    const res = await createSession({
      target_fasta: targetFasta,
      panel_fasta: panelFasta.trim() ? panelFasta : null,
    });
    sessionId.current = res.id;
    setInspect(
      await apiInspect({
        target_fasta: targetFasta,
        panel_fasta: panelFasta.trim() ? panelFasta : null,
      }),
    );
    return res.id;
  }, [targetFasta, panelFasta]);

  const runDesign = useCallback(
    async (override?: Partial<DesignParams>) => {
      setError(null);
      // Applied to state as well, so the criteria step shows what ran.
      const effective = override ? { ...params, ...override } : params;
      if (override) setParams(effective);
      if (!targetFasta.trim()) {
        // No input provided: stay on the representative sample.
        setDesign(sampleDesign);
        setSelectedRegion(sampleDesign.selected[0]?.region ?? null);
        setOffline(true);
        setStep("candidates");
        return;
      }
      setBusy(true);
      try {
        const id = await ensureSession();
        const res = await sessionDesign(id, effective);
        setDesign(res);
        setSelectedRegion(res.selected[0]?.region ?? null);
        setConvert(null);
        setOffline(false);
        setStep("candidates");
      } catch (e) {
        setError(e instanceof Error ? e.message : String(e));
        setDesign(sampleDesign);
        setSelectedRegion(sampleDesign.selected[0]?.region ?? null);
        setOffline(true);
        setStep("candidates");
      } finally {
        setBusy(false);
      }
    },
    [targetFasta, ensureSession, params],
  );

  const runConvert = useCallback(async () => {
    setError(null);
    if (offline || !sessionId.current) {
      setConvert(sampleConvert);
      setStep("construct");
      return;
    }
    setBusy(true);
    try {
      const res = await sessionConvert(sessionId.current, convertOptions);
      setConvert(res);
      setStep("construct");
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
      setConvert(sampleConvert);
      setStep("construct");
    } finally {
      setBusy(false);
    }
  }, [offline, convertOptions]);

  const selectedCandidate = useMemo(
    () => design?.selected.find((c) => c.region === selectedRegion),
    [design, selectedRegion],
  );
  const selectedConstruct = useMemo(() => {
    if (!selectedRegion) return undefined;
    return (
      convert?.constructs.find((c) => c.region === selectedRegion) ??
      (offline ? constructForRegion(selectedRegion) : undefined)
    );
  }, [convert, selectedRegion, offline]);
  const enrichment = useMemo(
    () => (offline && selectedRegion ? sampleEnrichment[selectedRegion] : undefined),
    [offline, selectedRegion],
  );

  const value: SessionValue = {
    step,
    setStep,
    expert,
    setExpert,
    targetFasta,
    setTargetFasta,
    panelFasta,
    setPanelFasta,
    params,
    updateParams,
    convertOptions,
    updateConvertOptions,
    inspect,
    design,
    convert,
    selectedRegion,
    setSelectedRegion,
    selectedCandidate,
    selectedConstruct,
    enrichment,
    offline,
    busy,
    error,
    runDesign,
    runConvert,
  };

  return <SessionContext.Provider value={value}>{children}</SessionContext.Provider>;
}
