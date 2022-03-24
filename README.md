The main objective of this work is to replicate the model presented in "[A simple model for control of COVID-19 infections on an urban campus][article]" by Robert Brown. This article presents a modified SEIR model on COVID-19 transmission on an urban campus and is verified using data collected in Fall 2020. Specifically, by exctracting the data used for the suggested study, we will run the same SEIR model and generate figures 1 to 3 to verify the replicability of the article.

[article]: https://...

R. Brown suggested a SEIR - susceptible-exposed-infectious-recovered - models with 5 compartments;
▪	Susceptibles (S)
▪	Infectious and symptomatic I(t)
▪	Infectious and asymptomatic A(t)
▪	Recovered R(t)
▪	Quarantine Q(t)
Where susceptibles individuals belonging to the campus community S(t) become infectious and asymptomatic at infection rate β. Infection is transmitted to the susceptible population in 3 ways; from infectious and symptomatic individuals I(t), infectious and asymptomatic individuals A(t) or come from the community outside the campus, which are considered exogenous cases and represented by the variable E(t). Individuals from either group become recovered at rate fR while infectious and asymptomatic individuals only are tested at rate fT. Finally, infectious and symptomatic persons go to quarantine at rate fQ. 

### Intellectual property

#### Article template

(C) 2015-2020, Nicolas Rougier + co-authors GPL-3+, Apache v2+, SIL Open Font License

This set of template files is free-licensed. The files contained in
the sub-directories roboto/; source-code-pro/; source-sans-pro/;
source-serif-pro/; have their free licences indicated with a
"*License.txt" file. All other files, including this one, are licensed
under the GPL version 3 or later, at your choosing, by Nicolas Rougier
and co-authors, 2015-2020. See the file COPYING for details of the
GPL-3 licence.

#### Model

