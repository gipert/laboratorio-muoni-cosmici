Decadimento dei muoni cosmici in alluminio
=========================================

Contents:
--------

* `muon.tex`:    file principale
* `analApp.tex`: appendici con i calcoli analitici e simulazioni
* `imgEff.tex`:  appendice con i grafici dell'efficienza
* `draw.tex`:    schema dell'apparato sperimentale
* `Scemi.tex`:   schemi dei circuiti
* `img`:         immagini

Collaborative Git:
-----------------

Ho creato due branches: `master` è la versione approvata e aggiornata della relazione, `modifiche` è riservata alle
proposte di modifica tramite pull request. Potete lavorare su `master` e caricare le modifiche direttamente sulla
versione finale oppure su `modifiche` per discuterne con gli altri.

La prima cosa da fare (dopo aver letto la man di `gittutorial` e configurato la propria identità con `git config --global --edit`) 
è copiare questa repository in locale tramite `git clone`; per cominciare a lavorare:

* spostatevi sulla branch `modifiche` con `git checkout modifiche`;
* fate le vostre modifiche e finalizzatele con `git commit -a`;
* caricate tutto in remoto con `git push bitbucket modifiche`;
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