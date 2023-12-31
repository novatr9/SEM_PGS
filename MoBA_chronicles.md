# MOBA Chronicles

This is a readme to document all information gathered from the start of the SEM_PGS project working on MoBa data. So far this readme is organised chronolically as an easy way for me to access information until completed, 
but it will one day be reorganised in sections at some point. \
Every person cited in this document has their email address linked to their name for ease of contact

Ps: Dates are european (ie DD.MM.YYYY)

# to do list:
 - ~~connect to the UiO server~~
 - ~~check if the NPSOL optimizer works on their server~~
 - find meta data about sequencing batches and batch effects of questionnaires gathering
 - get Rosa's scripts for calculating PGS
 - find partitioning script for hPGS
 - [google doc](https://docs.google.com/document/d/1VJSf6GZWipu-aCMFqXz3eeGK3aV1byOHL6vHYRvfMa0/edit?usp=sharing)

## 20.12.23

Singularities work on the cluster but not on TSD rooms
This [paper](https://www.biorxiv.org/content/10.1101/2022.06.23.496289v2) contains the QC steps for MoBa

## 18.12.23

Jeff made a singulrity to run the NPSOL optimizer in OpenMx. Dinka and I both imported it to our respective TDS rooms via the large file transfer link: [TSD import](https://data.tsd.usit.no/file-import/)

## 16.12.23

Dinka and I wrote a protocol for running SEM PGS on MoBa adding the tasks we’re working on at the moment. You can find it here: [google doc](https://docs.google.com/document/d/1VJSf6GZWipu-aCMFqXz3eeGK3aV1byOHL6vHYRvfMa0/edit?usp=sharing). It contains a step by step chronology of the files, packages/softwares, input and outputs we need for our analysis. This document is work in progress. 

## 12.12.23

Email from UiO's IT team: \
"Hi Noemie

Unfortunately we are not familiar with OpenMX. If you want to install openMX you need to contact your administrator and they will contact us. 

I am hereby closing this ticket."


## 11.12.23

Got access to the p805 teams chat

## 8.12.2023

Nice meeting with Eivind and Rosa. She felt a bit standoffish regarding sharing code, but it's most likely a matter of us giving appropriate reassurance that we will acknowledge her ownership. 
Xuanyu and Matt have access to the Microsoft teams channel of IuO but not yet me. [Clara Marie Fides Timpe](cmtimpe@uio.no) is admin and the person to request access to. \
We discussed possible phenotypes to start working on: 
 - **Height** 
 - **ADHD** is maternally reported and the mom self reports
 - **Depression** 8yo data for kids
Here the [link](https://github.com/rosacheesman/ADHD_schools) to Rosa's repo on analyzing MoBa and ADHD in schools \
\

**To connect to the Oslo server**: \
Go to the terminal, enter: \

`ssh p805-submit.tsd.usit.no` \
`cd /cluster/projects/p805/`

## 6.12.2023

Eivind could not make our monthly meeting. At least he let us know. \
This [link](https://www.fhi.no/en/ch/studies/moba/for-forskere-artikler/questionnaires-from-moba/#invitation-and-statements-of-consent) lists the questoinnaires included in MoBa.

## 4.12.23

We got access!!!!! [Dinka](dinka.smajlagic@psykologi.uio.no), a postdoc from UiO, has arrived to the lab for 10 months and she is very helpful navigating and understand how the data is organised in room p-805 (our project). 

## 25.11.2023

Still no access to the data, Matt and I have agreed that I should stop harassing Eivind (I even invited him to a zoom meeting but he never logged on), and I should focus on the AM project until we meet again.

## 23.11.2023

Found this github repo to [QC MoBa data](https://github.com/novatr9/MoBaPsychGen-QC-pipeline)

## 21.11.2023

Message from Rosa: \
"We have PRS for about 25 traits now, but these are just using normal LDpred and we haven’t distinguished non transmitted haplotypes. Here you can see where we’re at with making scores and 
what the sumstats we used are [Link text](https://docs.google.com/spreadsheets/d/1Jn_NgXWPQsHLCjCW-Pa1Bs21hIZtzQKdZUjsfkhMLZI/edit#gid=488070425)
For MDD we used the broad depression GWAS by Howard et al."
She does not know what is wrong with our access to data. I have contacted Eivind but have gotten no reponse. I email the [TSD help desk](tsd-drift@usit.uio.no), they are responsive but confirmed only Eivind can grant us the access.

## 20.11.2023

Eivind sent us our norsk credentials through WHatsapp and [Xuanyu](Xuanyu.Lyu@colorado.edu) and I got access to the TSD portal (through which we access sensitive data hosted by UiO), but no data to be seen anywhere. \
To connect remotly, download [VMware Horizon Client](https://my-horizon.vmware.com/portal/webclient/index.html). In the set up, the server to add is view.tsd.usit.no \
For more info, got to the UiO TSD help [page](https://www.uio.no/english/services/it/research/sensitive-data/help/project-access.html) \
We are access the project/ room p805 \
NOTE: the remote desktop has a scandinavian keyboard so when thinking of a password keep that in mind. For example, the sign "@" cannot be typed from a US keyboard.

## 19.11.2023

First meeting with [Eivind Ystom](eivind.ystrom@psykologi.uio.no) and [Rosa Cheeseman](r.c.g.cheesman@psykologi.uio.no). Eivind is a professor at UiO and Rosa his postdoc. He will grant us access to the data 
unformally as the legal way is so long and tedious. 



