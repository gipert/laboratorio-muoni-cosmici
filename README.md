Decadimento dei muoni cosmici in alluminio
=========================================

Contents:
--------

`README.md`
`env.sh`: shell script da eseguire per rendere visibili gli eseguibili ovunque
`Makefile`: Makefile generale per compilare tutto il progetto

`latex/`

* `muon.tex`:    file principale
* `analApp.tex`: appendici con i calcoli analitici e simulazioni
* `imgEff.tex`:  appendice con i grafici dell'efficienza
* `draw.tex`:    schema dell'apparato sperimentale
* `Scemi.tex`:   schemi dei circuiti
* `img`:         immagini

`code/`

* `efficiency/`
    * `eff.cc`:    macro ROOT per il plot delle efficienze
    * `*.txt`:     data files per ogni PMT
* `calibration/`
    * `cal.C` :    macro ROOT per calibrare i TAC
    * `TAC_437`:   dati calibrazione per il TAC nr. 437
    * `TAC_467`:   dati calibrazione per il TAC nr. 467
* `analysis/`
    * `lifetime_nocalib_bkgr_sub.cc` : macro ROOT per l'analisi (parziale 1° semestre)
    * `dati_semestre_1` : directory contenenti tutti i file usati nell'analisi del 1° semestre
    * `dati_semestre_2` : directory contenente tutti i file usati nell'analisi del 2° semestre
* `montecarlo/`
    * `montecarlo.cc` : codice simulazione MC


Collaborative Git:
-----------------

Ho creato due branches: `master` è la versione approvata e aggiornata della relazione, `modifiche` è riservata alle
proposte di modifica tramite pull request. Potete lavorare su `master` e caricare le modifiche direttamente sulla
versione finale oppure su `modifiche` per discuterne con gli altri.

La prima cosa da fare (dopo aver letto la man di `gittutorial` e configurato la propria identità con `git config --global --edit`) 
è copiare questa repository in locale tramite `git clone`; per cominciare a lavorare:

* spostatevi sulla branch `modifiche` con `git checkout modifiche`;
* fate le vostre modifiche e finalizzatele con `git commit -a`;
* caricate tutto in remoto con `git push bitbucket modifiche` (può essere che voi dobbiate ridefinire l'alias
  `bitbucket` con il comando `git remote add bitbucket git@bitbucket.org:luigipertoldi/relazionemuoni.git`);
* andate sulla pagina della repository e create una pull request. Chiunque può lavorare su `modifiche` e caricare le
  proprie modifiche lasciando la pull request in sospeso. Una volta che si è d'accordo si unisce il contenuto di
  `modifiche` con quello di `master`.

Nota Bene:

* Ricordatevi di dire a `git` di ignorare file output di compilazione aggiungendoli al file `.gitignore`
* Ricordatevi di mantenere sincronizzato il vostro lavoro con quello presente sulla repository tramite `git pull`

La procedura di cui sopra è un esempio dei tanti modi in cui si possono proporre modifiche ad un progetto, `git` è un
programma molto esteso, sbizzarritevi! Potete creare più branches, potete lavorare sulla vostra versione nella vostra
repo e poi proporne l'unione con il master ecc... Leggete la man di `git` per qualunque dubbio. 

Provate a modificare questo README.md per esercitarvi!

SSH config:
----------

Avete bisogno di una coppia di chiavi ssh, pubblica e privata (se non esiste il file `~/.ssh/id_rsa` allora non ce 
l'avete). Per generarle:

    ssh-keygen -t rsa

e seguite le istruzioni (potete lasciare tutto default, la passphrase non è necessaria). Infine dovete dare a Bitbucket
la vostra chiave pubblica, dimodochè la vostra macchina sia associata alla vostra identità, copiando il contenuto di
`~/.ssh/id_rsa.pub` nell'apposito form tra le impostazioni del vostro account.

`man` codici:
-------------

* `lifetime_nocalib_bkgr_sub.cc`: 
	`$ root -l .x lifetime_nocalib_bkgr_sub.cc( "[filelistName]" , [rebinFactor] = 12 , [midValue] = 1000 )`
