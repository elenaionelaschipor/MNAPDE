EDP = EXTERIOR DIRICHLET PROBLEM 



abbiamo omega un dominio lipschiz, k >0, gd \\in H^1/2(Gamma). u (u scatt) soddisfa il problema esterno (nel senso che ho Omega = Omega+ unito Omega - , dove omega + =  esterno, omega- = oggetto contro cui faccio lo scattering) SE



laplaciano u + k^2 u = 0 (Omega+)

traccia u = gd (bordo) => el caso di un ostacolo sound soft ho che gd = - traccia di u incidente

u uscente / radiating / sommerfeld





per risolvere EDP uso metodo BEM (boundary element method) che a differenza di FEM lavora solamente sul bordo gamma e torna sul resto alla fine. FEM è bellino, ma ho dimensione + alta. 

lavoro sul bordo!



𝒮 = single layer potential, è fatto in modo che:



(𝒮ψ)(x) = integrale su Gamma di soluzione fondamentale\*ψ dy, al variare di x in Omega + (parte esterna)



(Sψ)(x) = integrale su gamma di soluzione fondamentale\*ψ dy, al variare di x in Gamma!!! (bordo)



abbiamo che la traccia di (𝒮ψ) è uguale a Sψ quindi risolvere EDP è come risolvere



&#x09;						cerco ψ tale che Sψ = gd in Gamma.   (✨)

&#x20;

a quel punto, trovata ψ, ho u = (𝒮ψ) (in Omega! però uso ψ che è in Gamma, quindi semplifico i conti)



da ✨ posso proseguire cercando una forma variazionale del tipo:



&#x09;cerco ψ in H-1/2 (omega) tale che	A(ψ, xi) := <Sψ, xi> = <gd, xi> =: F(xi) (❤️)



dove A ed F sono la forma sesquilineare e antilineare del problema variazionale (va tutto bene per la teoria)



su (❤️) posso fare o Collocazione oppure Galerkin, ottengo una roba del tipo  **A \* ψ = F** dove cambiano le matrici in base al metodo. 



per entrambi i metodi fisso uno spazio V\_N contenuto in H -1/2 (Gamma) e cerco ψN che risolve la ✨



faccio una partizione di Gamma in una mesh T\_N(gamma) di N segmenti K1...KN disgiunti e prendo VN = funzioni costanti a tratti sulla mesh. Come funzioni di base prendiamo le indicatrici delle KN



* COLLOCAZIONE

  * devo risolvere (SψN)(xj) = gd(xj)

    * alla fine della fiera ci viene che **A\_col \* ψ = F\_col**
    * dove A\_Col\_j,m = integrale su Km di sol fondamntale(xj, y) ds (y), e F\_col, j = gd(xj)
* Galerkin

  * devo risolvere A(ψN, xiN) = F(xiN)

    * A Gal jm = integrale doppio su Kj e Km di sol fondamentale
    * F gal j = integrale su j di gd.













