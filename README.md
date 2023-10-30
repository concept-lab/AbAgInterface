# AbAgInterface

Ab-Ag Interface annotation of the Structural antibody database ([SAbDab](https://academic.oup.com/nar/article/42/D1/D1140/1044118)).

SabDab PDB structures present in the `structures` folder can be obtained [here](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/search/?all=true#downloads).

You may find a compressed file of `structures` folder (including both SadDab and Nanoshaper results) [here](https://u.pcloud.link/publink/show?code=XZCRiaXZrRr2fY5mUj89IwYOBqsaHLA00xfk).

Luca:
To reproduce data, first download structure folder and run get_exposed.py.
Then one can use explore.ipynb to get stats.
Surface patches used in the paper can be reproduced using the scripts in the patch_generation folder. Bear in mind that one has to select pertinent structures in the filename within the script (prior having generated those via the get_exposed script)
