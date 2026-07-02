"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""


def minimise_allele(alt: str, ref: str) -> str:
    """Converts a VCF allele string into a minimised SPDI format.

    This function converts a VCF allele representation that omits anchoring bases for prediction
    scores in the INFO column. It removes the anchoring base if the first base of the reference
    and alternate alleles is identical. If the resulting allele is empty, a hyphen ('-') is returned.

    Args:
        alt (str): The alternate allele from the VCF data.
        ref (str): The reference allele from the VCF data.

    Returns:
        str: The minimised allele in SPDI format.
    """
    minimised_allele_string = alt
    if ref[0] == alt[0]:
        minimised_allele_string = alt[1:] if len(alt) > 1 else "-"
    return minimised_allele_string


def minimise_protein_sequence(
    ref: str, alt: str, start: int, end: int, length: int
) -> tuple[str, str, int, int, int]:
    """Minimise a protein sequence change

    Args:
        ref (str): The reference amino acid sequence.
        alt (str): The alternate amino acid sequence.
        start (int): The start position of sequence.
        end (int): The end position of sequence.
        length (int): Length of sequence

    Returns:
        tuple[str, str, int, int, int]: The minimised reference and alternate protein sequences, positions and length.
    """

    # handling only protein deletion as this can be ambiguous in UI
    if len(ref) <= len(alt) or ref == "-" or alt == "-":
        return (ref, alt, start, end, length)

    prefix_length = 0
    while prefix_length < len(alt) and ref[prefix_length] == alt[prefix_length]:
        prefix_length += 1

    min_ref = ref[prefix_length:]
    min_alt = alt[prefix_length:]

    suffix_length = 0
    while min_ref and min_alt and min_ref[-1] == min_alt[-1]:
        suffix_length += 1
        min_ref = min_ref[:-1]
        min_alt = min_alt[:-1]

    if prefix_length or suffix_length:
        start = int(start) + prefix_length
        end = int(end) - suffix_length
        length = length - prefix_length - suffix_length

    min_ref = min_ref or "-"
    min_alt = min_alt or "-"
    return (min_ref, min_alt, start, end, length)


def decode_population_name(name: str):
    return name.replace("$2C", ",")
