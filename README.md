Decadimento dei muoni cosmici in alluminio
=========================================

Contents:
--------

* `README.md`
* `env.sh`: shell script da eseguire per rendere visibili gli eseguibili ovunque
* `Makefile`: Makefile generale per compilare tutto il progetto

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
    * `semestre1` : macro e notebooks di Mathematica per la simulazione MonteCarlo delle coincidenze
    * `semestre2`
        * `baselineStart.cc` : studio della stima della baseline con fit "pol0" al variare del punto iniziale
                               nello spettro simulato exp+cost
        * `fitexpStart.cc` : studio della bontà del fit esponenziale al variare del punto iniziale nello spettro
                             simulato exp+exp+cost
        * `montecarlo_modifiedforbaseline.cc` : confronto metodi di fit della baseline
	* `montecarlo_SingleNsim.cc` : simulazione montecarlo completa
        * `montecarlo.cc` : funzione-simulazione montecarlo completa adatta all'iterazione in `main.cc`
	* `main.cc` : matrice di simulazioni

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
* `git pull` semplicemente aggiunge ciò che trova nel remote al vostro lavoro, se volete essere certi di avere una 
  lo stesso identico codice che si trova nella branch di riferimento nella repository dovete:
    * `$ git fetch bitbucket` : fa un download di quello che c'è nella repo senza fare il merge con il locale
    * `$ git reset --hard bitbucket/branch_name` : effettivamente rimpiazza ciò che avete in locale
* Quando c'è una nuova branch nel remote e volete averla in locale dovete prima crearla e associarla a quella remota.
  Solo dopo potete fare pull:
    * `$ git checkout -b newlocalbranchname bitbucket/branch-name`
    * `$ git pull`
* Se vi siete accorti che avete fatto una cazzata nell'ultimo commit e volete tornare indietro potete:
    * `$ git reset HEAD~` : annulla il commit
    * `$ git revert HEAD~` : crea un nuovo commit in cui eliminate le modifiche dell'ultimo

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
