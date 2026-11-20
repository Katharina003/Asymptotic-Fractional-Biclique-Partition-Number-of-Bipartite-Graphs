"""Verify an exact primal certificate for the Domino graph's Kronecker power.

Usage: python tools/verify_primal.py path/to/resultK1.txt
Checks 
- biclique realizability, 
- matrix columns, 
- nonnegative rational weights,
- zero equality slacks, 
- Mx=edgeorbitsizes, and
- the exact objective.
This proves primal feasibility (an upper bound), not global optimality.
"""

import argparse
from fractions import Fraction
from itertools import combinations_with_replacement
from math import factorial
from pathlib import Path
import re


def require(condition, message):
    """Stop verification with an explanation if a check fails."""
    if not condition:
        raise ValueError(message)


def read_field(text, label):
    """Extracts all lines that start with the label."""
    matches = []
    for line in text.splitlines():
        if line.startswith(label):
            matches.append(line[len(label):].strip())
    require(len(matches) == 1, f"Expected exactly one field: {label}")
    return matches[0]


def read_section(text, start, end=None):
    """Extracts the text between two section headings (start and end contain starting index and ending index)."""
    require(text.count(start) == 1,
            f"Missing or repeated section: {start.strip()}")
    remainder = text.split(start, 1)[1]

    if end is None:
        return remainder

    require(end in remainder, f"Missing section end: {end.strip()}")
    return remainder.split(end, 1)[0]


def read_nodes(text, side, k):
    """Return a dictionary: bit string -> number of selected nodes."""
    # side is either "Top" or "Bot" such that string matches
    line = read_field(text, side + " Nodes: ")
    nodes = {}
    # Example node: (N: 010, #2, K: 3)
    # find all nodes
    pattern = r"\(N\s*:\s*([01]+)\s*,\s*#(\d+)\s*,\s*K:\s*(\d+)\s*\)"
    matches = list(re.finditer(pattern, line))
    require(matches, f"Missing {side} nodes")
    position = 0

    for match in matches:
        # no broken strings inbetween
        require(not line[position:match.start()].strip(), f"Malformed {side} nodes")
        name, count, power = match.groups()
        #check whether node format is okay
        require(int(power) == k and len(name) == k, f"Wrong dimension in {side} node {name}")
        require(int(count) > 0, f"Nonpositive count in {side} node {name}")
        # nodes should be completely compressed
        require(name not in nodes, f"Duplicate {side} node type {name}")
        # save how many nodes of orbit exist in biclique
        nodes[name] = int(count)
        #go to next position in string
        position = match.end()
    require(not line[position:].strip(), f"Malformed {side} nodes")
    return nodes


def read_result(text):
    """Read the current Solver.cpp output."""
    k = int(read_field(text, "power: "))
    require(k >= 1, "Kronecker power must be positive")
    numEdgeOrbits = int(read_field(text, "MBasis: ")) # in the result file behind MBasis stands its dimension. As it is quadratic only one number is printed.
    require(int(read_field(text, "numEdgeOrbits: ")) == numEdgeOrbits,# part of the domino graph fields
            "Graph and matrix dimensions differ")
    
    # extract matrix
    matrix_text = read_section(text, f"MBasis: {numEdgeOrbits}\n", "Dominograph:")
    matrix = [[int(number) for number in line.split()]
              for line in matrix_text.strip().splitlines()]

    # parse all bicliques that are part of the final solution
    heading = re.search(r"(?m)^Solution Bicliques:[ \t]*\n", text)
    require(heading is not None, "Missing Solution Bicliques section")
    biclique_text = read_section(text, heading.group(), "MBasis:")
    bicliques = []
    position = 0
    for match in re.finditer(r"biclique number (\d+):\{(.*?)\n\}", biclique_text, re.S):
        # check whether there is undefined text between two bicliques
        require(not biclique_text[position:match.start()].strip(), "Malformed biclique block")
        index, block = match.groups()
        bicliques.append({
            "index": int(index),
            "top": read_nodes(block, "Top", k),
            "bottom": read_nodes(block, "Bot", k),
        })
        position = match.end()
    require(not biclique_text[position:].strip(), "Malformed biclique block")

    # check whether the number of bicliques is the expected one
    numBasicVars = len(bicliques)
    require(numBasicVars <= numEdgeOrbits, "Invalid number of basic variables")

    # read exact fractions for each variable
    headings = list(re.finditer(r"(?m)^Solution values:[ \t]*\n", text))
    require(len(headings) == 1, "Expected exactly one exact solution section")
    value_text = text[headings[0].end():].split("Dual solution:", 1)[0]
    entries = []
    for line in value_text.strip().splitlines():
        match = re.fullmatch(r"([xs])_(\d+) = (-?\d+) / (-?\d+)", line.strip())
        require(match is not None, f"Malformed solution value: {line}")
        #kind: x (normal variable) or s (slackvariable)
        kind, index, numerator, denominator = match.groups()
        require(int(denominator) != 0, "Invalid denominator")
        index = int(index)
        if kind == "s":
            column = len(entries)
            # A slack s_i belongs to edge orbit i, so its column must be e_i.
            require(0 <= index < numEdgeOrbits, "Invalid slack orbit index")
            require(len(matrix) == numEdgeOrbits and
                    all(len(row) == numEdgeOrbits for row in matrix), "Invalid matrix dimensions") #is checked later again (at a point that is not in the parser, and thus prettier) but needed here to prevent index errors
            # slack should have exactly one 1 at i
            require(column < numEdgeOrbits and
                    all(row[column] == int(i == index) for i, row in enumerate(matrix)),
                    "Slack column does not match its orbit index")
            # choose slack indices 0 to num slack variables instead of i (which is not consecutive between variables)
            index = column - numBasicVars
        entries.append((kind, index, Fraction(int(numerator), int(denominator))))
    return k, numEdgeOrbits, numBasicVars, matrix, bicliques, entries


def read_certificate(path):
    """Read the certificate."""
    text = Path(path).read_text(encoding="utf-8")
    return read_result(text)





# verification methods

def check_dimensions(numEdgeOrbits, numBasicVars, matrix, bicliques, entries):
    """Check that the certificate contains the expected rows, columns and values."""
    require(len(matrix) == numEdgeOrbits, "base should have |[E]| rows")
    for row in matrix:
        require(len(row) == numEdgeOrbits, "base should have |[E]| columns")
    require(len(bicliques) == numBasicVars, "not every variable (ignoring slack variables) has a corresponding biclique")
    require(len(entries) == numEdgeOrbits, "x vector not long enough")

    # check that every biclique identifier is unique, such that the solution cannot be misread
    indices = set()
    for biclique in bicliques:
        require(biclique["index"] not in indices, "to bicliques have the same index")
        indices.add(biclique["index"])


def check_node_capacities(nodes, opposite_nodes, k):
    """Check how many distinct nodes of each type can fit in this biclique."""

    # check for each node orbit, whether it is possible to select orbitCount many
    # different nodes in that orbit
    for orbit, orbitCount in nodes.items():
        capacity = 1 # numer of distinct nodes in the orbit
        # multiply the number of possibilities for each coordinate
        for coordinate in range(k):

            if orbit[coordinate] == "0":
                # check whether all nodes need to be assigned to 0, because of edges in
                #(. . 1)
                #(. . .)
                #(1 . .)
                # (if we do not fix it to 0 then some of the nodes will not be neighbors)
                oppositeInV_0 = False
                for opposite_bits in opposite_nodes:
                    if opposite_bits[coordinate] == "0":
                        oppositeInV_0 = True
                # when we do not need to fix the node, we can either choose 0 or 2
                if not oppositeInV_0:
                    capacity *= 2
        require(orbitCount <= capacity,
                f"node type {orbit} needs {orbitCount} nodes, but only {capacity} fit")


def check_bicliques(bicliques, k):
    """Check that each descriptor can be realized as a biclique in D^k.

    Nodes in the program are represented as condensed nodes in Section 4.3.3 of the thesis.
    We now explain how we select concrete nodes in each orbit.

    If both left and right nodes (top and bottom) have outer vertices in a coordinate, 
    It means an edge in a corner of the matrix is selected (one of the two shown below.)
    ( . . 1)
    ( . . .)
    ( 1 . .)
    then we can use mirroring symmetrie to project the biclique to a top left biclique.
    So we fix all of them to (left 0, right 5). This is important, else some combinations
    of left and right nodes might not be neighboring.

    If in one coordinate all combinations of left and right nodes
    have at least one entry L_1^1 or R_1^1, then all nodes are adjacent and we can select any node.
    ( . 1 .)
    ( 1 . 1)
    ( . 1 .)
    we count whether it is possible to select enough distinct nodes fulfilling 
    the just described rules.
    """
    for biclique in bicliques:
        top = biclique["top"]
        bottom = biclique["bottom"]
        require(top and bottom, "the zero matrix is not a valid biclique")

        # to reduce duplicated code
        #checks node capacities for first entrie, thus 2 calls.
        check_node_capacities(top, bottom, k)
        check_node_capacities(bottom, top, k)


def edge_orbit(top_node, bottom_node):
    """Find the edge orbit from the two node types."""
    # bit pairs 00=a, 10=b, 01=d, 11=c
    letters = []
    for top, bottom in zip(top_node, bottom_node):
        if top == "0" and bottom == "0":
            letters.append("a")
        elif top == "1" and bottom == "0":
            letters.append("b")
        elif top == "0" and bottom == "1":
            letters.append("d")
        else:
            letters.append("c")
    return "".join(sorted(letters))


def check_columns(matrix, orbitNames, bicliques):
    """Recompute each column by counting its biclique's edges."""
    # iterate over all bicliques and reconstruct edge orbit names from node orbits
    for j, biclique in enumerate(bicliques):
        column = {}
        for top, top_count in biclique["top"].items():
            for bottom, bottom_count in biclique["bottom"].items():
                # construct name and multiplicity
                name = edge_orbit(top, bottom)
                column[name] = column.get(name, 0) + top_count * bottom_count
        # compare new matrix with given matrix
        for i, name in orbitNames.items():
            require(matrix[i][j] == column.get(name, 0),
                    f"Column {j}, row {i} does not match its biclique")


def check_weights(entries, bicliques, numBasicVars):
    """Check variable indices, nonnegative weights, and zero equality slacks."""
    values = []
    for j, (kind, index, value) in enumerate(entries):
        if j < numBasicVars:
            require(kind == "x", "found slack variable instead of normal one")
            require(index == bicliques[j]["index"], "variable claims to be basic but does not belong to basis")
            require(value >= 0, "All weights must be positive")
        else:
            require(kind == "s", "Expected a slack variable, not a normal variable")
            require(value == 0, "slack variables must be zero, else Mx<1 for at least one row in M")
        values.append(value)
    return values


def check_equations(matrix, orbitSizes, values, numBasicVars):
    """Check Mx=b exactly, one row at a time. All slack values are already zero."""
    for i in range(len(matrix)):
        total = Fraction(0)
        for j in range(numBasicVars):
            total += matrix[i][j] * values[j]
        require(total == orbitSizes[i], f"Mx=b fails in row {i}")


def check_objective(values, numBasicVars, objective):
    """The objective is the sum of the biclique weights."""
    total = Fraction(0)
    for j in range(numBasicVars):
        total += values[j]
    require(total == objective, "Wrong exact objective")


def verify(path):
    """Read the data, then perform each check in order."""
    # read the data
    k, numEdgeOrbits, numBasicVars, matrix, bicliques, entries = read_certificate(path)
    # compute additional needed data
    # the orbit names and sizes are not printed in result file, so we need to compute them, as described in the master's thesis
    orbitNames = dict(enumerate("".join(letters)
                                   for letters in combinations_with_replacement("abcd", k)))
    orbitSizes = []
    for name in orbitNames.values():
        denominator = 2 ** name.count("c")
        for letter in "abcd":
            denominator *= factorial(name.count(letter))
        orbitSizes.append(2**k * factorial(k) // denominator)
    require(numEdgeOrbits == len(orbitSizes), "Wrong number of orbitNames")
    require(numEdgeOrbits == len(orbitNames), "Wrong number of orbitNames")

    # The printed minimized value is rounded, the results of the thesis are calculated by using the verifier. Compute exact objective.
    objective = sum((value for kind, _, value in entries if kind == "x"), Fraction(0))

    # now check everything
    # first check that dimensions of M and x fit, and that every entry has a corresponding biclique entry
    check_dimensions(numEdgeOrbits, numBasicVars, matrix, bicliques, entries)
    # then check wether all bicliques are really bicliques
    check_bicliques(bicliques, k)
    # then recompute M from the now verified bicliques, to approve its columns
    check_columns(matrix, orbitNames, bicliques)
    # check x>=0 and slack = 0. 
    values = check_weights(entries, bicliques, numBasicVars)
    # check that Mx = edgeorbit sizes
    check_equations(matrix, orbitSizes, values, numBasicVars)
    # check the objective value.
    check_objective(values, numBasicVars, objective)
    return objective


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)

    # read file names (at least one, but also multiple at the same time possible)
    parser.add_argument("certificates", nargs="+")
    certificates = parser.parse_args().certificates
    with Path("certificates.txt").open("w", encoding="utf-8") as output:
        for certificate in certificates:
            objective = verify(certificate)
            line = f"{certificate}: VERIFIED; exact primal objective = {objective} (approximately {float(objective):.6f})"
            print(line)
            print(line, file=output)
