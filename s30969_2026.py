# s30969
# 12.05.2026
# generator losowych sekwencji DNA

import random
import csv


def validate_positive_int(prompt: str,
                          min_val: int = 1,
                          max_val: int = 100_000) -> int:
                            
    while True:
        value = input(prompt)
        try:
            number = int(value)

            if min_val <= number <= max_val:
                return number

            print(f"Błąd: wartość musi być liczbą całkowitą "
                  f"z zakresu [{min_val}, {max_val}].")

        except ValueError:
            print(f"Błąd: wartość musi być liczbą całkowitą "
                  f"z zakresu [{min_val}, {max_val}].")

#walidacja id sekwencji fasta
def validate_seq_id() -> str:
    while True:
        seq_id = input("Podaj ID sekwencji: ")

        if " " in seq_id or "\t" in seq_id:
            print("Błąd: ID nie może zawierać białych znaków.")
        elif seq_id == "":
            print("Błąd: ID nie może być puste.")
        else:
            return seq_id


def generate_sequence(length: int) -> str:
    nucleotides = ["A", "C", "G", "T"]

    return "".join(random.choice(nucleotides)
                   for _ in range(length))

#oblicza statystyki sekwencji
def calculate_stats(sequence: str) -> dict:
    length = len(sequence)

    stats = {}

    for nucleotide in ["A", "C", "G", "T"]:
        count = sequence.count(nucleotide)
        stats[nucleotide] = round((count / length) * 100, 2)

    gc = sequence.count("G") + sequence.count("C")

    stats["gc_ratio_A"] = round((gc / length) * 100, 2)

    return stats

#wstawia imie w losowe miejsce sekwencji
def insert_name(sequence: str, name: str) -> str:
    position = random.randint(0, len(sequence))

    return (sequence[:position] +
            name.lower() +
            sequence[position:])

#formatuje rekord fasta
def format_fasta(seq_id: str,
                 description: str,
                 sequence: str,
                 line_width: int = 80) -> str:
    header = f">{seq_id}"

    if description:
        header += f" {description}"

    lines = [header]

    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i:i + line_width])

    lines.append("# EOF_1")

    return "\n".join(lines)

#zapisuje fasta do pliku
def save_fasta(filename: str, content: str):
    with open(filename, "w") as file:
        file.write(content)

#wyszukuje motyw w sekwencji
def find_motif(sequence: str, motif: str) -> list:
    positions = []

    for i in range(len(sequence) - len(motif) + 1):

        if sequence[i:i + len(motif)] == motif:
            positions.append(i + 1)

    return positions

#tworzy sekwencjie komplementarna
def complementary_sequence(sequence: str) -> str:
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    result = ""

    for nucleotide in sequence:
        result += complement[nucleotide]

    return result

#tworzy sekwencjie odwrotnie komplementarna
def reverse_complement(sequence: str) -> str:
    comp = complementary_sequence(sequence)

    return comp[::-1]

#Tworzy sekwencję mRNA
def transcribe_mrna(sequence: str) -> str:
    return sequence.replace("T", "U")


def sliding_window_gc(sequence: str,
                      window_size: int,
                      output_csv: str):
                        
    #oblicza GC - content w oknach przesuwnych
    with open(output_csv, "w", newline="") as csvfile:

        writer = csv.writer(csvfile)

        writer.writerow(["pozycja_startu", "gc_content"])

        for i in range(len(sequence) - window_size + 1):

            window = sequence[i:i + window_size]

            gc = ((window.count("G") +
                   window.count("C")) / window_size) * 100

            writer.writerow([i + 1, round(gc, 2)])


def main():
    length = validate_positive_int(
        "Podaj długość sekwencji: "
    )

    seq_id = validate_seq_id()

    description = input("Podaj opis sekwencji: ")

    name = input("Podaj imię: ")

    sequence = generate_sequence(length)

    stats = calculate_stats(sequence)

    sequence_with_name = insert_name(sequence, name)

    fasta_content = format_fasta(
        seq_id,
        description,
        sequence_with_name
    )

    filename = f"{seq_id}.fasta"

    save_fasta(filename, fasta_content)

    print(f"\nSekwencja zapisana do pliku: {filename}")

    print(f"\nStatystyki sekwencji (n={length}):")

    for nucleotide in ["A", "C", "G", "T"]:
        print(f"  {nucleotide}: {stats[nucleotide]:.2f}%")

    print(f"  GC-content: {stats['gc_ratio_A']:.2f}%")

  #dodatkowe funkcjonalności:

    motif = input("\nPodaj motyw do wyszukania: ").upper()

    positions = find_motif(sequence, motif)

    if positions:
        print("Motyw znaleziony na pozycjach:")
        print(positions)
    else:
        print("Nie znaleziono motywu.")

    comp_seq = complementary_sequence(sequence)

    rev_comp_seq = reverse_complement(sequence)

    mrna = transcribe_mrna(sequence)

    extra_records = "\n\n"

    extra_records += format_fasta(
        seq_id + "_complement",
        "Komplementarna",
        comp_seq
    )

    extra_records += "\n\n"

    extra_records += format_fasta(
        seq_id + "_reverse_complement",
        "Odwrotnie komplementarna",
        rev_comp_seq
    )

    extra_records += "\n\n"

    extra_records += format_fasta(
        seq_id + "_mRNA",
        "mRNA",
        mrna
    )

    with open(filename, "a") as file:
        file.write(extra_records)

    print("\nDodano:")
    print("- sekwencję komplementarną")
    print("- sekwencję odwrotnie komplementarną")
    print("- sekwencję mRNA")

    window_size = validate_positive_int(
        "\nPodaj rozmiar okna sliding window: ",
        1,
        length
    )

    csv_name = f"{seq_id}_gc.csv"

    sliding_window_gc(
        sequence,
        window_size,
        csv_name
    )

    print(f"Wyniki sliding window zapisano do: {csv_name}")


if __name__ == "__main__":
    main()
