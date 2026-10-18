/**
 * The publications behind the design rules, as cited in docs/design-rules.md and
 * docs/shrna-conversion.md. "short" is for tooltips, "full" for the reference page.
 */

export interface Source {
  short: string;
  full: string;
}

export const SOURCES = {
  elbashir2001: {
    short: "Elbashir 2001",
    full: "Elbashir et al. 2001, Nature 411:494; EMBO J 20:6877",
  },
  uitei2004: { short: "Ui-Tei 2004", full: "Ui-Tei et al. 2004, NAR 32:936" },
  reynolds2004: { short: "Reynolds 2004", full: "Reynolds et al. 2004, Nat Biotechnol 22:326" },
  amarzguioui2004: {
    short: "Amarzguioui & Prydz 2004",
    full: "Amarzguioui & Prydz 2004, BBRC 316:1050",
  },
  khvorova2003: { short: "Khvorova 2003", full: "Khvorova et al. 2003, Cell 115:209" },
  schwarz2003: { short: "Schwarz 2003", full: "Schwarz et al. 2003, Cell 115:199" },
  naito2009: { short: "Naito 2009", full: "Naito et al. 2009, BMC Bioinformatics 10:392" },
  uitei2008: { short: "Ui-Tei 2008", full: "Ui-Tei et al. 2008, NAR 36:7100" },
  salomon2015: { short: "Salomon 2015", full: "Salomon et al. 2015, Cell 162:84" },
  wee2012: { short: "Wee 2012", full: "Wee et al. 2012, Cell 151:1055" },
  freier1986: { short: "Freier 1986", full: "Freier et al. 1986, PNAS 83:9373" },
  xia1998: { short: "Xia 1998", full: "Xia et al. 1998, Biochemistry 37:14719" },
  groebe1988: { short: "Groebe & Uhlenbeck 1988", full: "Groebe & Uhlenbeck 1988, NAR 16:11725" },
  huang2009: { short: "Huang 2009", full: "Huang et al. 2009, NAR 37:7560" },
  schwarz2006: { short: "Schwarz 2006", full: "Schwarz et al. 2006, PLoS Genet 2:e140" },
  du2005: { short: "Du 2005", full: "Du et al. 2005, NAR 33:1671" },
  brummelkamp2002: {
    short: "Brummelkamp 2002",
    full: "Brummelkamp et al. 2002, Science 296:550",
  },
  jensen2012: { short: "Jensen 2012", full: "Jensen et al. 2012, PLoS ONE 7:e43095" },
  moffat2006: { short: "Moffat 2006", full: "Moffat et al. 2006, Cell 124:1283" },
  wang2013: { short: "Wang 2013", full: "Wang et al. 2013, PLoS ONE 8:e56110" },
  gu2012: { short: "Gu 2012", full: "Gu et al. 2012, Cell 151:900" },
  ge2010: { short: "Ge 2010", full: "Ge et al. 2010, RNA 16:106" },
  bogenhagen1981: {
    short: "Bogenhagen & Brown 1981",
    full: "Bogenhagen & Brown 1981, Cell 24:261",
  },
  gao2018: { short: "Gao 2018", full: "Gao et al. 2018, Mol Ther Nucleic Acids 10:36" },
  gao2017: { short: "Gao 2017", full: "Gao et al. 2017, Transcription 8:275" },
} satisfies Record<string, Source>;

export type SourceKey = keyof typeof SOURCES;

/** "Ui-Tei 2004; Reynolds 2004", for a tooltip. */
export function shortCites(keys: SourceKey[]): string {
  return keys.map((k) => SOURCES[k].short).join("; ");
}
