# Numer albumu: 30969
# Data: 22.06.2026
# Opis programu:
# Generator losowych sekwencji nukleotydowych (DNA) z zapisem w formacie FASTA.
# Program pobiera od użytkownika długość sekwencji, procentowy rozkład
# nukleotydów (A/C/G/T), ID, opis, imię oraz motyw do wyszukania, a następnie:
#  - generuje pseudolosową sekwencję DNA o zadanym rozkładzie,
#  - wstawia imię użytkownika (małymi literami) w losowej pozycji,
#  - zapisuje sekwencję wraz z nicią komplementarną i odwrotnie
#    komplementarną do pliku multi-FASTA ({ID}.fasta),
#  - oblicza i wyświetla statystyki składu nukleotydowego oraz GC-content,
#  - generuje wykres GC-content (analiza okna przesuwnego) i zapisuje go
#    do pliku PNG,
#  - wyszukuje wskazany motyw i wypisuje pozycje jego wystąpień (od 1).


import random

# Import matplotlib zabezpieczony - jeśli biblioteki nie ma, program nadal działa,
# pomijając jedynie generowanie wykresu (reszta funkcji pozostaje dostępna)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_DOSTEPNY = True
except ImportError:
    MATPLOTLIB_DOSTEPNY = False


def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """Pobiera od użytkownika liczbę całkowitą z zakresu [min_val, max_val].
    W przypadku błędnego wejścia (nie-liczba lub wartość spoza zakresu)
    ponawia pytanie zamiast kończyć program wyjątkiem.
    """
    while True:
        dane = input(prompt).strip()
        try:
            value = int(dane)
        except ValueError:
            print(f"Błąd: wartość musi być liczbą całkowitą z zakresu [{min_val}, {max_val}].")
            continue
        if value < min_val or value > max_val:
            print(f"Błąd: wartość musi być liczbą całkowitą z zakresu [{min_val}, {max_val}].")
            continue
        return value


def get_nucleotide_distribution() -> dict:
    """Pobiera procentowy udział każdego nukleotydu (A, C, G, T).
    Waliduje, czy suma udziałów wynosi 100%. Zwraca słownik {nukleotyd: %}.
    W razie błędu (nie-liczba, wartość spoza [0,100], suma != 100) ponawia pytanie.
    """
    while True:
        dist = {}
        valid = True
        for nuc in "ACGT":
            # replace(",", ".") - akceptuje też polski przecinek dziesiętny (np. 30,5)
            dane = input(f"Podaj udział procentowy {nuc} [%]: ").strip().replace(",", ".")
            try:
                value = float(dane)
            except ValueError:
                print(f"Błąd: udział {nuc} musi być liczbą.")
                valid = False
                break
            if value < 0 or value > 100:
                print(f"Błąd: udział {nuc} musi być z zakresu [0, 100].")
                valid = False
                break
            dist[nuc] = value
        if not valid:
            continue

        total = sum(dist.values())
        # tolerancja 1e-9 chroni przed błędami zaokrąglenia liczb zmiennoprzecinkowych
        if abs(total - 100) > 1e-9:
            print(f"Błąd: suma udziałów musi wynosić 100% (podano {total:.2f}%).")
            continue
        return dist


def generate_sequence(length: int) -> str:
    """Zwraca losową sekwencję DNA o zadanej długości (rozkład jednostajny A/C/G/T)."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def generate_sequence_weighted(length: int, dist: dict) -> str:
    """Generuje sekwencję DNA o zadanym rozkładzie procentowym nukleotydów.
    Wykorzystuje random.choices z wagami - losowanie jest probabilistyczne,
    więc faktyczny skład sekwencji jest tylko zbliżony do zadanego rozkładu.
    """
    nucleotides = list(dist.keys())
    weights = [dist[n] for n in nucleotides]
    return ''.join(random.choices(nucleotides, weights=weights, k=length))


def calculate_stats(sequence: str) -> dict:
    """Zwraca słownik ze statystykami sekwencji.
    Klucze: "A", "C", "G", "T" (udziały procentowe, float) oraz
    "GC" (zawartość GC, czyli suma udziałów G i C, float).
    """
    length = len(sequence)
    stats = {nuc: sequence.count(nuc) / length * 100 for nuc in "ACGT"}
    stats["GC"] = stats["G"] + stats["C"]
    return stats


def complement(sequence: str) -> str:
    """Zwraca nić komplementarną (zamiana par A<->T, C<->G), w tym samym kierunku."""
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return ''.join(pairs[nuc] for nuc in sequence.upper())


def reverse_complement(sequence: str) -> str:
    """Zwraca nić odwrotnie komplementarną.
    To nić komplementarna odczytana od końca (odwrócona przez [::-1]),
    co odpowiada antyrównoległej drugiej nici DNA w kierunku 5'->3'.
    """
    return complement(sequence)[::-1]


def insert_name(sequence: str, name: str) -> str:
    """Wstawia imię w losową pozycję sekwencji (zapisane małymi literami).
    Dzięki małym literom imię jest wizualnie odróżnialne od nukleotydów
    (wielkie litery) i nie jest traktowane jako część sekwencji biologicznej.
    """
    pos = random.randint(0, len(sequence))
    return sequence[:pos] + name.lower() + sequence[pos:]


def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """Zwraca sformatowany rekord FASTA jako string.
    Nagłówek zaczyna się od '>' i zawiera ID oraz opis rozdzielone spacją.
    Sekwencja jest łamana na linie o szerokości line_width (domyślnie 80 znaków).
    """
    header = f">{seq_id} {description}".rstrip()
    lines = [header]
    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i:i + line_width])
    return "\n".join(lines) + "\n"


def find_motif(sequence: str, motif: str) -> list:
    """Wyszukuje wszystkie wystąpienia motywu i zwraca ich pozycje (indeksowanie od 1).
    Przesunięcie start = idx + 1 (a nie idx + len(motif)) sprawia, że
    wykrywane są również wystąpienia nakładające się (np. "AA" w "AAAA").
    """
    positions = []
    motif = motif.upper()
    start = 0
    while True:
        idx = sequence.find(motif, start)
        if idx == -1:
            break
        positions.append(idx + 1)   # +1 -> konwencja biologiczna (od 1, nie od 0)
        start = idx + 1
    return positions


def get_valid_id() -> str:
    """Pobiera ID sekwencji. ID nie może być puste ani zawierać białych znaków."""
    while True:
        seq_id = input("Podaj ID sekwencji: ").strip()
        if not seq_id:
            print("Błąd: ID nie może być puste.")
            continue
        if any(c.isspace() for c in seq_id):
            print("Błąd: ID nie może zawierać białych znaków.")
            continue
        return seq_id


def get_valid_name() -> str:
    """Pobiera imię użytkownika. Imię nie może być puste."""
    while True:
        name = input("Podaj imię: ").strip()
        if not name:
            print("Błąd: imię nie może być puste.")
            continue
        return name


def get_valid_motif() -> str:
    """Pobiera motyw do wyszukania (niepusty, złożony tylko ze znaków A, C, G, T)."""
    while True:
        motif = input("Podaj motyw do wyszukania (np. ATG): ").strip().upper()
        if not motif:
            print("Błąd: motyw nie może być pusty.")
            continue
        if any(c not in "ACGT" for c in motif):
            print("Błąd: motyw może zawierać tylko znaki A, C, G, T.")
            continue
        return motif


def sliding_window_gc(sequence: str, window: int):
    """Liczy GC-content w oknie przesuwnym wzdłuż sekwencji.
    Zwraca dwie listy: pozycje (środek każdego okna) oraz wartości GC% dla okien.
    Jeśli okno jest dłuższe niż sekwencja, zostaje skrócone do jej długości.
    """
    n = len(sequence)
    if window > n:
        window = n
    positions = []
    gc_values = []
    for i in range(0, n - window + 1):
        frag = sequence[i:i + window]
        gc = (frag.count("G") + frag.count("C")) / window * 100
        positions.append(i + window // 2)  # środek okna jako pozycja na osi X
        gc_values.append(gc)
    return positions, gc_values


def plot_gc_content(sequence: str, window: int, filename: str, seq_id: str = ""):
    """Generuje wykres liniowy GC-content (analiza okna przesuwnego) i zapisuje do pliku PNG."""
    positions, gc_values = sliding_window_gc(sequence, window)
    overall_gc = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100

    plt.figure(figsize=(10, 5))
    plt.plot(positions, gc_values, color="tab:blue", linewidth=1, label="GC% (okno)")
    plt.axhline(overall_gc, color="tab:red", linestyle="--", linewidth=1,
                label=f"Średnia GC = {overall_gc:.2f}%")
    plt.title(f"GC-content wzdłuż sekwencji {seq_id} (okno = {window} nt)".rstrip())
    plt.xlabel("Pozycja w sekwencji [nt]")
    plt.ylabel("GC-content [%]")
    plt.ylim(0, 100)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()


def main():
    """Główna funkcja programu - obsługuje interakcję z użytkownikiem i scala wszystkie kroki."""
    #pobranie danych wejściowych z walidacją
    length = validate_positive_int("Podaj długość sekwencji: ")
    dist = get_nucleotide_distribution()
    seq_id = get_valid_id()
    description = input("Podaj opis sekwencji: ").strip()
    name = get_valid_name()
    motif = get_valid_motif()

    #generacja sekwencji o zadanym rozkładzie (czysta sekwencja = podstawa analiz)
    sequence = generate_sequence_weighted(length, dist)
    stats = calculate_stats(sequence)

    #nić komplementarna i odwrotnie komplementarna (z czystej sekwencji DNA)
    # liczone z 'sequence', NIE z wersji z imieniem - imię nie jest nukleotydem
    # i nie posiada pary komplementarnej (complement rzuciłby KeyError).
    comp_seq = complement(sequence)
    revcomp_seq = reverse_complement(sequence)

    # wersja z imieniem (tylko do zapisu rekordu głównego, nie wpływa na statystyki)
    output_sequence = insert_name(sequence, name)

    #zapis do pliku multi-FASTA: rekord główny + komplementarny + odwrotnie komplementarny
    filename = f"{seq_id}.fasta"
    with open(filename, "w") as f:
        f.write(format_fasta(seq_id, description, output_sequence))
        f.write(format_fasta(f"{seq_id}_comp", "nic komplementarna", comp_seq))
        f.write(format_fasta(f"{seq_id}_revcomp", "nic odwrotnie komplementarna", revcomp_seq))

    #wykres GC-content (okno ~10% długości, min. 10 nt)
    window = max(1, min(length, max(10, length // 10)))
    plot_filename = f"{seq_id}_gc.png"
    if MATPLOTLIB_DOSTEPNY:
        plot_gc_content(sequence, window, plot_filename, seq_id)
        wykres_info = f"Wykres GC-content zapisany do pliku: {plot_filename}"
    else:
        wykres_info = "Wykres pominięty (brak biblioteki matplotlib - zainstaluj: pip install matplotlib)."

    #wyszukiwanie motywu (w czystej sekwencji, bez imienia)
    motif_positions = find_motif(sequence, motif)

    #wypisanie wyników
    print(f"Sekwencja zapisana do pliku: {filename}")
    print(f"  (rekordy: {seq_id}, {seq_id}_comp, {seq_id}_revcomp)")
    print(wykres_info)
    print(f"Statystyki sekwencji (n={length}):")
    for nuc in "ACGT":
        print(f"{nuc}: {stats[nuc]:.2f}%")
    print(f"GC-content: {stats['GC']:.2f}%")

    #podgląd nici skracany do 60 znaków dla czytelności (pełne wersje są w pliku)
    print(f"\nNić komplementarna:           {comp_seq if length <= 60 else comp_seq[:60] + '...'}")
    print(f"Nić odwrotnie komplementarna: {revcomp_seq if length <= 60 else revcomp_seq[:60] + '...'}")

    print(f"\nWyszukiwanie motywu '{motif}':")
    if motif_positions:
        print(f"Liczba wystąpień: {len(motif_positions)}")
        print(f"Pozycje (od 1): {motif_positions}")
    else:
        print("Brak wystąpień motywu w sekwencji.")


if __name__ == "__main__":
    main()
