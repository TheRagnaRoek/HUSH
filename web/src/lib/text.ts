/** Upper-case the first letter, leaving the rest as given. */
export function cap(s: string): string {
  return s.charAt(0).toUpperCase() + s.slice(1);
}

/** Keep GC as an acronym wherever an engine name is shown. */
export function acronyms(s: string): string {
  return s.replace(/\bgc\b/gi, "GC");
}

/** An engine name ("gc", "allele hit") as display text at the start of a line. */
export function gateLabel(s: string): string {
  return cap(acronyms(s));
}
