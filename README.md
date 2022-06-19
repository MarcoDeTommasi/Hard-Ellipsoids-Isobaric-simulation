# Isobaric simulation of hard ellipsoids (C++)

In questa simulazione è studiato il comportamento di un gas di ellissoidi duri attraverso un algoritmo
Metropolis-Montecarlo in condizioni isobariche (NPT). Dall’analisi dei risultati ottenuti nella simulazione,
si è ricavata l’equazione di stato che è stata in seguito confrontata con quella ottenuta dallo sviluppo del
viriale al primo ordine.

 ![Eq_stato](https://user-images.githubusercontent.com/79840407/174475243-dcadaf98-fabb-44b9-b5fd-d0710b8d87fe.svg)
 
Il volume escluso ![CodeCogsEqn (1)](https://user-images.githubusercontent.com/79840407/174475411-7a5354ad-7050-47d3-bb68-569d0a3b8b24.svg) è stato stimato a sua volta attraverso un’ integrazione Montecarlo. Gli ellissoidi
utilizzati nella simulazione sono uniassiali con aspect ratio ![CodeCogsEqn](https://user-images.githubusercontent.com/79840407/174475390-eafbd1ce-45ea-4a2b-b32e-5f6e3f3a073a.svg) = 2. I semiassi sono b = c = 1 e a = ![CodeCogsEqn](https://user-images.githubusercontent.com/79840407/174475390-eafbd1ce-45ea-4a2b-b32e-5f6e3f3a073a.svg) .
Il volume degli ellissoidi è quindi ottenuto come ![CodeCogsEqn (2)](https://user-images.githubusercontent.com/79840407/174475428-9e8cdffa-fe60-451d-84c1-7a332115fbd6.svg) . Oltre alla stima dell’equazione di stato, è
stata ricavata la funzione di distribuzione radiale g(r).

## Alcuni Risultati:

![E0S400000](https://user-images.githubusercontent.com/79840407/174475543-f81abead-2c54-4c49-a015-c1057b042c2c.png)


![grUltima](https://user-images.githubusercontent.com/79840407/174475455-0f164ea8-919e-4ed1-aa0d-ab5fe4fa5f69.png)


### Analisi completa:
[Ellissoidi.pdf](https://github.com/MarcoDeTommasi/Hard-Ellipsoids-Isobaric-simulation/files/8935084/DeTommasiNPTsim.pdf)
