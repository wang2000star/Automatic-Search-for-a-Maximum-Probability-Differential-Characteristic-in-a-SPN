Reference: "Automatic Search for a Maximum Probability Differential Characteristic in a Substitution-Permutation Network"
<https://eprint.iacr.org/2016/652.pdf>
Code description: Based on SPN structure using branch and bound method for optimal differential feature search, I wrote the code mainly for 64 bit PRESENT password search, but not limited to PRESENT. The PRESENT uses N identical S-boxes, so the differential probability transfer table is the same. If the cases of different S-boxes are dealt with, the code can be appropriately adjusted. In addition, considering that the linear layer of the PRESENT is a bit replacement, the decision conditions are increased, so it can be said that the algorithm designed in the paper is fully implemented.



