# Searching for reaction pathways in formose data

## 1. Generate reaction network:
Contains reactions for:
- enolate formation,
- enolate protonation,
- aldol addition of formaldehyde to enolates,
- Cannizzaro reactions and
- enolate reactions with sugars of C<sub><4</sub>.

## 2. Remove nodes
Formaldehyde, hydroxide and water are ignored as nodes in pathways searches (CaCl<sub>2</sub> is not included in the framework).

## 3. For each data set:
  - Make a list of detected compounds, include which compounds cannot be detected.
  - This list is ordered by amplitude (highest amplitude first)
  - Only reactions which include compounds in this list as reactants are included in the pathway search.
  - Remove reactions in which inputs are products

## 4. Connect inputs (reactants) to products
Find the shortest pathways from the carbon input to the compound with the highest amplitude .

## 5. [Version 2] Connect products to each other.

iterate down the amplitude list, find shortest pathways to every compound downstream
of the node.

## 6. Check that all products have reactions leading to them
If they do not, cycle backwards from them in the sorted amplitude list until a connection is found.
