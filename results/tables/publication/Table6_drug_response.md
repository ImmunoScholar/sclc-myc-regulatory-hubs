# Table 6. BET-inhibitor sensitivity versus MYC-family expression

| Screen | Drug | Gene | n | rho | p | rho_NEadj | p_NEadj |
|---|---|---|---|---|---|---|---|
| GDSC2 | BIRABRESIB | MYCL | 30 |  0.621 | 0.000253 |  0.631 | 0.000187 |
| GDSC2 | MOLIBRESIB | MYCL | 30 |  0.617 | 0.000281 |  0.587 | 0.000648 |
| GDSC1 | PFI-1 | MYCL | 29 |  0.426 | 0.021300 |  0.429 | 0.020400 |
| GDSC2 | BIRABRESIB | MYC | 30 | -0.390 | 0.033000 | -0.314 | 0.091200 |
| GDSC1 | MOLIBRESIB | MYCN | 31 | -0.351 | 0.052900 | -0.307 | 0.092700 |
| GDSC2 | MOLIBRESIB | MYC | 30 | -0.326 | 0.078400 | -0.231 | 0.219000 |
| GDSC1 | PFI-1 | MYC | 29 | -0.279 | 0.142000 | -0.183 | 0.341000 |
| PRISM | MOLIBRESIB | MYCN | 21 |  0.306 | 0.178000 |  0.266 | 0.243000 |

**Note.** AUC is a sensitivity metric on which lower means more sensitive, so positive rho means MYCL-high lines are LESS sensitive. This is the only association in the project that survives lineage adjustment, and is reported as exploratory: it does not replicate in PRISM, and GDSC1 and GDSC2 are the same platform over overlapping lines rather than two independent screens.
