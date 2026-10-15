"""Paper-wide palette — single source of truth (chosen by Alejandra 2026-08-28, option C).

Three semantic roles, used identically in every figure:
  FOCAL   #1F8A8C  teal   — this study's dark/novel/focal set (dark families, candidate models,
                            the catalogue's own genes)
  CONTROL #EE6C4D  coral  — matched/blinded control populations (annotated twins, known folds,
                            control absences); also the warm pole where a panel needs a pair
                            with no control present (e.g. AMR vs metal, sun-exposed surface)
  GREY    #999999         — known/reference objects and de-emphasis (superposition targets,
                            alternative denominators)

pLDDT ramp: single-hue teal, pale (low confidence) -> deep (high); monotonic lightness, never
collides with CONTROL. Compartment trio harmonised to the palette (no focal/control in those
panels): Surface warm coral, Deep teal, PlantRoot soft violet.
CVD: teal/coral separation passes with wide margin (approx. Lab dE 83; see palette sheet
2026-08-28); re-validate if any value changes.
"""
FOCAL = "#1F8A8C"
FOCAL_LIGHT = "#63AFAD"   # variant of the same measure (e.g. space-partial beside marginal)
CONTROL = "#EE6C4D"
REFKNOWN = "#8A6FB8"   # known reference OBJECT (superposition target, DB reference tips)
GREY = "#999999"
LIGHTGREY = "#bbbbbb"
FILLGREY = "#e2e2e2"
INK = "#222222"
MUTED = "#666666"
GRID = "#dddddd"

COMPARTMENT = {"Surface": "#EE6C4D", "Deep": "#1F8A8C", "PlantRoot": "#8A6FB8"}

PLDDT_STOPS = ["#E8F2F1", "#63AFAD", "#0F5D5F"]   # low -> high

def plddt_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("plddt_teal", PLDDT_STOPS)

# aliases matching the generators' local names
C_DARK = C_FOCAL = C_HI = C = FOCAL
C_CTRL = C_COMP = CONTROL
