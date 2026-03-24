import re
from io import StringIO

import pandas as pd
import matplotlib.pyplot as plt
import requests


# ---------- FONTES OFICIAIS HARMONIZADAS ----------
PRIMARY_URL = (
    "https://climate.copernicus.eu/sites/default/files/2026-01/"
    "GCH2025_PR_Fig1_timeseries_annual_global_temperature_anomalies_preindustrial_data_0.csv"
)

# fallback: o Met Office mantém um painel oficial com as mesmas 6 séries
# já em diferença anual em relação ao pré-industrial
FALLBACK_URL = "https://www.metoffice.gov.uk/hadobs/monitoring/global-temperature.html"


EXPECTED = [
    "Year",
    "ERA5",
    "JRA-3Q",
    "GISTEMPv4",
    "NOAAGlobalTempv6",
    "Berkeley Earth",
    "HadCRUT5",
]


def download_text(url: str) -> str:
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.text


def find_header_and_read_csv(txt: str) -> pd.DataFrame:
    """
    Procura a linha real de cabeçalho dentro de um texto CSV que possa conter
    linhas descritivas antes da tabela.
    """
    lines = txt.splitlines()
    header_idx = None

    for i, line in enumerate(lines):
        line_clean = line.strip().replace('"', "")
        if "Year" in line_clean and "ERA5" in line_clean:
            header_idx = i
            break

    if header_idx is None:
        raise ValueError("Cabeçalho com 'Year' e 'ERA5' não encontrado.")

    csv_text = "\n".join(lines[header_idx:])

    # tenta vírgula
    df = pd.read_csv(StringIO(csv_text))
    if "Year" not in [str(c).strip() for c in df.columns]:
        # tenta ponto e vírgula
        df = pd.read_csv(StringIO(csv_text), sep=";")

    df.columns = [str(c).strip() for c in df.columns]
    return df


def try_copernicus() -> pd.DataFrame:
    txt = download_text(PRIMARY_URL)
    df = find_header_and_read_csv(txt)
    return df


def try_metoffice_fallback() -> pd.DataFrame:
    """
    Fallback conservador:
    baixa a página do Met Office e tenta localizar um link .csv associado
    ao painel de 'Annual global mean near-surface temperature difference
    from pre-industrial conditions'.
    """
    html = download_text(FALLBACK_URL)

    # procura links CSV absolutos ou relativos
    matches = re.findall(r'https?://[^"\']+\.csv', html)
    if not matches:
        rel = re.findall(r'["\']([^"\']+\.csv)["\']', html)
        matches = [
            "https://www.metoffice.gov.uk" + m if m.startswith("/") else m
            for m in rel
        ]

    # tenta encontrar algum csv que contenha as séries esperadas
    for url in matches:
        try:
            txt = download_text(url)
            df = pd.read_csv(StringIO(txt))
            df.columns = [str(c).strip() for c in df.columns]
            if "Year" in df.columns and ("ERA5" in df.columns or "HadCRUT5" in df.columns):
                return df
        except Exception:
            continue

    raise ValueError("Fallback do Met Office não encontrou um CSV utilizável.")


def load_harmonized_dataframe() -> pd.DataFrame:
    try:
        df = try_copernicus()
        source = "Copernicus"
    except Exception as e1:
        print(f"Falha na fonte principal ({e1}). Tentando fallback...")
        df = try_metoffice_fallback()
        source = "Met Office fallback"

    print(f"Fonte usada: {source}")
    print("Colunas originais:", df.columns.tolist())

    # mantém só colunas esperadas, se existirem
    keep = [c for c in EXPECTED if c in df.columns]
    if "Year" not in keep:
        raise ValueError("A coluna 'Year' não foi encontrada após a leitura.")

    df = df[keep].copy()

    # conversões numéricas
    df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
    for c in df.columns:
        if c != "Year":
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["Year"]).copy()
    df["Year"] = df["Year"].astype(int)

    return df.sort_values("Year").reset_index(drop=True)


def save_outputs(df: pd.DataFrame) -> pd.DataFrame:
    df.to_excel("ecmwf_like_temperature_wide.xlsx", index=False)

    long_df = df.melt(
        id_vars="Year",
        var_name="Dataset",
        value_name="Temp_C_above_preindustrial",
    ).dropna()

    long_df.to_excel("ecmwf_like_temperature_long.xlsx", index=False)
    return long_df


def plot_ecmwf_style(long_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(11.5, 7.2))

    # outras fontes em cinza
    other = long_df[long_df["Dataset"] != "ERA5"]
    for _, g in other.groupby("Dataset"):
        ax.scatter(
            g["Year"],
            g["Temp_C_above_preindustrial"],
            s=10,
            facecolors="white",
            edgecolors="gray",
            linewidths=0.9,
            zorder=2,
        )

    # ERA5 em barras coloridas + linha
    era5 = long_df[long_df["Dataset"] == "ERA5"].sort_values("Year")
    if not era5.empty:
        vals = era5["Temp_C_above_preindustrial"].clip(lower=0)
        vmin, vmax = vals.min(), vals.max()
        norm = (vals - vmin) / (vmax - vmin) if vmax != vmin else vals * 0
        colors = plt.cm.YlOrRd(norm)

        ax.bar(
            era5["Year"],
            era5["Temp_C_above_preindustrial"],
            width=0.9,
            color=colors,
            edgecolor="none",
            zorder=1,
        )

        ax.plot(
            era5["Year"],
            era5["Temp_C_above_preindustrial"],
            linewidth=1.2,
            zorder=3,
        )

    # linha de 1,5 °C
    ax.axhline(1.5, color="gray", linestyle="--", linewidth=1)
    ax.text(
        long_df["Year"].max() + 1,
        1.5,
        "+1.5°C",
        va="center",
        ha="left",
        color="gray",
    )

    ax.set_xlim(long_df["Year"].min() - 3, long_df["Year"].max() + 5)
    ax.set_ylim(-0.25, 1.75)
    ax.set_xlabel("")
    ax.set_ylabel("Global annual surface air temperature increase above pre-industrial (°C)")
    ax.set_title(
        "2025 was the third-warmest year on record according to ERA5",
        loc="left",
        weight="bold",
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="0.85", linewidth=0.8)

    ax.text(1885, 0.22, "Other sources*", color="gray")
    ax.text(2004, 0.55, "ERA5 data", color="black", weight="bold")

    foot = (
        "*Other sources comprise JRA-3Q, GISTEMPv4, NOAAGlobalTempv6, "
        "Berkeley Earth, HadCRUT5. Data for 2025 are only available for ERA5 and JRA-3Q.\n"
        "Reference period: pre-industrial (1850-1900) • Credit: C3S/ECMWF"
    )
    fig.text(0.06, 0.02, foot, fontsize=9)

    plt.tight_layout(rect=[0, 0.06, 1, 1])
    plt.savefig("ecmwf_like_temperature_figure.png", dpi=300, bbox_inches="tight")
    plt.savefig("ecmwf_like_temperature_figure.svg", bbox_inches="tight")
    plt.show()


def main():
    df = load_harmonized_dataframe()
    print(df.head())
    print(df.tail())
    print(df.columns.tolist())

    long_df = save_outputs(df)
    plot_ecmwf_style(long_df)


if __name__ == "__main__":
    main()
