---
title: 'OMEGA-Py: Python Tools for OMEGA Data'
tags:
 - Python
 - Astronomy
 - Planetary Science
 - OMEGA
 - Mars Express
 - Mars
 - ESA
authors:
 - name: Aurélien Stcherbinine
   orcid: 0000-0002-7086-5443
   affiliation: 1
 - name: Yves Langevin
   orcid: 0000-0002-4492-215X
   affiliation: 2
 - name: John Carter
   orcid: 0000-0002-2698-6926
   affiliation: "2, 3"
 - name: Mathieu Vincendon
   affiliation: 2
   orcid: 0000-0002-1072-7487
 - name: Yann Leseigneur
   orcid: 0000-0003-1787-2924
   affiliation: 2
 - name: Océane Barraud
   orcid: 0000-0002-9985-1109
   affiliation: 4
affiliations:
 - name: Department of Astronomy and Planetary Science, Northern Arizona University, Flagstaff, AZ USA
   index: 1
 - name: IAS, Université Paris-Saclay, CNRS, Orsay, France
   index: 2
 - name: LAM, Université Aix-Marseille, CNRS, CNES, Marseille, France
   index: 3
 - name: German Aerospace Center (DLR), Institute of Planetary Research, Berlin, Germany
   index: 4
date: 24 October 2023
bibliography: paper.bib
---

# Summary
`OMEGA-Py` is a Python 3 module dedicated to the scientific use of data provided by the 
Observatoire pour la Minéralogie, l'Eau, les Glaces et l'Activité (OMEGA) instrument onboard
the ESA Mars Express (MEx) orbiter [@bibring_2004].
It has been developed as an alternative to the IDL routines [@soft10]
of the OMEGA legacy software provided by the instrument team for the past 20 years.

The module notably includes a re-implementation of the most recent release of
the IDL OMEGA software (v10, `SOFT 10`),
which performs the reading, calibration and reduction of the level 1B data publicly
available on the ESA [PSA](https://archives.esac.esa.int/psa/#!Table%20View/OMEGA=instrument)
[@besse_2018] to produce level 2A data that can be used for the scientific analysis,
but also contains several additional
data reduction functions such as built-in atmospheric and thermal corrections
(using previously published methods) and graphics tools including interactive
visualization of the data or generation of composite OMEGA maps using the
`matplotlib` module, including geographic projection [@hunter_2007].

The objective of the module is to facilitate the scientific exploitation of OMEGA observations,
especially for the younger generation of planetary scientists who are more used to the Python
language than IDL. Plus, the presence of built-in correction and visualization functions
aims at making the huge and very complete OMEGA dataset (rich of 20 years of observations now) more 
easily accessible.

Since its first release in 2020, `OMEGA-Py` has been used in published studies
[@stcherbinine_2021b; @leseigneur_2023] as well as in currently ongoing projects
[e.g., @barraud_2022a].

`OMEGA-Py` can be installed from PyPI with `pip install omegapy`, and
is distributed as an official software by the OMEGA team
since the release of version 3.0 deployed in October 2023.


# Statement of need
The accessibility of data returned by space missions is a crucial point to ensure the
development of open science.
While the OMEGA dataset is public with yearly releases as of 2023, the legacy pipeline uses
a proprietary software and several crucial data reduction algorithms are not public,
thus severely hindering its use.
Since the beginning of the science phase of the OMEGA instrument in 2004, the instrument
team has provided 10 releases of the IDL software (`SOFT01` to `SOFT10`) to read the
level 1B binary files that can be downloaded from the ESA [PSA](https://archives.esac.esa.int/psa/#!Table%20View/OMEGA=instrument) 
and generate level 2A data with reflectance spectra.
However, the presence of an IDL solution only may raise some concerns:

 * The cost of an IDL license, as it is a proprietary language, makes it not accessible to everyone.
 * As the community (and especially the youngest generation) is moving to use mostly Python
   instead of IDL, the requirement to use the IDL language to access OMEGA data can limit its
   accessibility.

In addition, over the past years, the OMEGA dataset had a reputation in the community
for being challenging and requiring a lot of investment to use.
With `OMEGA-Py` we aim to tackle this reputation by providing a free all-in-one toolbox
to load, correct, analyze, and visualize the OMEGA data, and thus make the unique OMEGA dataset
rich of 20 years of observations easily accessible to the community and especially to
the younger generation of scientists and students.


# Acknowledgements
We thank all the people who helped with testing and improving the module.


# References

