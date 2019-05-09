# Dynamics of colistin resistance in chicken farms of the Mekong delta, Vietnam

## Data

The data is generated by the [make_data.Rmd](https://github.com/viparc/colistin_resistance/blob/master/make_data.Rmd)
Rmarkdown script (the HTML version of which can be seen on 
[rpubs.com/choisy/colistindata](http://rpubs.com/choisy/colistindata)) and is
stored in the
[data/viparc.csv](https://raw.githubusercontent.com/viparc/colistin_resistance/master/data/viparc.csv)
CSV file. It contains 5391 weeks of observation (rows) and 149 variables
(columns). The variables are

* **USUBJID:** farm ID;
* **FLOCKSEQUENCE:** flock ID (in a given farm);
* **WEEK:** week number (in a given flock of a given farm);
* **RESPIRATORY, ..., SUDDENDEATH:** presence/absence of 6 clinical signs in the flock;
* **amoxicillin_use, ..., unknown_use:** presence/absence of 44 antimicrobials
in the flock;
* **amoxicillin_g, ..., unknown_g:** mass, in g, of the antimicriobial used in
the flock;
* **amoxicillin_g.kg, ..., unknown_g.kg:** mass, in g per kg of chicken, of the antimicriobial used in the flock;
* **completed:** boolean informing whether the flock is done (`TRUE`) or still ongoing (`FALSE`);
* **CHICKENTOTAL:** total number of chicken in the farm;
* **arm_weight_kg:** mass, in kg, of chicken in the farm;
* **sampling:** boolean informing whether there is feces sampling during the
week.
