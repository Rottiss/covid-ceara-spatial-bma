"""Display labels for the main figures; data identifiers remain unchanged."""
import os
import re

LANGUAGE = os.environ.get("FIGURE_LANGUAGE", "en").lower()
if LANGUAGE not in {"en", "pt"}:
    raise ValueError("FIGURE_LANGUAGE must be 'en' or 'pt'.")


def text(english, portuguese):
    return english if LANGUAGE == "en" else portuguese


OUTCOMES = {
    "casos": text("Cases", "Casos"),
    "obitos": text("Deaths", "Óbitos"),
    "letalidade": text("Case fatality", "Letalidade"),
}

COVARIATES = {
    "z_pop": "Population",
    "z_area": "Land area",
    "z_densidade": "Population density",
    "z_pop.rural": "Rural population",
    "z_idosos": "Population aged 60+",
    "metrop": "Fortaleza metro region",
    "semiarido": "Semi-arid region",
    "z_sus1k": "SUS health facilities",
    "z_leitos1k": "SUS hospital beds",
    "z_prof1k": "SUS health professionals",
    "z_imuni": "Immunization coverage",
    "z_int.circ": "Circulatory admissions",
    "z_int.resp": "Respiratory admissions",
    "z_int.diab": "Diabetes admissions",
    "z_int.asma": "Asthma admissions",
    "z_tax.agua": "Water-supply coverage",
    "z_desp.saude": "Health and sanitation spending",
    "z_pibpc": "GDP per capita",
    "z_emprego": "Formal employment",
    "z_ex.pobr": "Extreme poverty",
    "z_bolsaf": "Bolsa Família beneficiaries",
    "z_ivs.infra": "IVS urban infrastructure",
    "z_ideb.fund": "IDEB (grades 5/9)",
    "z_ideb3": "IDEB, upper secondary",
    "z_desp.educ": "Education and culture spending",
    "z_energia100k": "Electricity consumption",
    "z_tax.hom": "Homicide rate",
    "aliadogov": "Federal party alignment",
    "z_votos": "Mayor's vote share",
    "z_eleitas.fem": "Elected women",
}


def covariate_label(identifier, portuguese):
    return COVARIATES[identifier] if LANGUAGE == "en" else portuguese


def scenario_label(label):
    """Translate exported scenario labels without changing grouping keys."""
    if LANGUAGE == "pt":
        return label
    for source, target in (
        ("Realista, quatro verdadeiras", "Realistic, 4 active"),
        ("Realista, oito verdadeiras", "Realistic, 8 active"),
        ("Baixa autocorrelação", "Low autocorrelation"),
        ("Alta autocorrelação", "High autocorrelation"),
        ("Controle", "Control"),
        ("Realista", "Realistic"),
        ("verdadeiras", "true covariates"),
    ):
        label = label.replace(source, target)
    return re.sub(r"(?<=\d),(?=\d)", ".", label)
