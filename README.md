# Sbox-depth-evaluation
Depth evaluation of 4-bit Sbox

There are four files in this folder:
sbox_depth_test.py
all_expr_hw8_within_depth4.5
readme.txt
depth3.txt

The sbox_depth_test.py program outputs the minimum depth of a given S-box. 
all_expr_hw8_within_depth4.5.txt contains all the 4-bit balanced boolean functions (from depth 0) up to depth 4.5 that is loaded in python data structure called dictionary.
depth3.txt contains all 4-bit optimal Sboxes whose depth are 3. Here "optimal" we mean the highest probability of differential characteristics and linear bias are all 2^{-2}.

-------------------------------------------------------------------
How to test depth of a given S-box:
-------------------------------------------------------------------

1. Write the given S-box in a list as input to the function find_expression().
2. Run sbox_depth_test.py in a python GUI or a terminal.

We give an example for testing depth of Sb0 in Midori in the program. Running the program you will get:

Given S-box:
[12, 10, 13, 3, 14, 11, 15, 7, 8, 9, 1, 5, 0, 2, 4, 6]
a' = 
['( ( a ) AND ( b ) ) NOR ( ( NOT ( c ) ) NOR ( ( a ) NOR ( d ) ) )', 3, 5.0]
b' = 
['( ( ( a ) NOR ( d ) ) NOR ( ( b ) AND ( c ) ) ) NAND ( ( a ) NAND ( ( c ) AND ( d ) ) )', 3.5, 7.0]
c' = 
['( ( b ) NOR ( d ) ) NOR ( NOT ( ( a ) NAND ( ( b ) NAND ( d ) ) ) )', 3.5, 4.5]
d' = 
['( ( ( a ) NAND ( b ) ) NAND ( ( c ) OR ( d ) ) ) NOR ( ( a ) NOR ( ( b ) OR ( c ) ) )', 3.5, 7.0]
Depth is 3.5.

where a|b|c|d is the 4-bit input to the given S-box and a'|b'|c'|d' is the 4-bit output. Each logical expression is followed by its depth and gate size.
