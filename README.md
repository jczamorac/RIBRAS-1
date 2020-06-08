# RIBRAS

Essa é uma simulação do RIBRAS (Radioactive Ion Beam In Brazil) feita com a ferramenta GEANT4. Foi o projeto de iniciação científica de alunos de graduação em física da USP. O programa visa simular a física de uma reação nuclear.

### Requisitos

Entre os requisitos estão as ferramentas:

* GEANT4 
* CLHEP
* ROOT 
* HepRApp

A maior parte do código foi feita em C++ e alguns scripts feitos em bash script.

### Inicialização

Antes de rodar a simulação, é necessário alterar o arquivo `vis.mac`, nele é possível alterar valores como a energia e posição do feixe primário de partículas, a posição do alvo, a reação nuclear, ligar e desligar o campo eletromagnético etc. Além disso, é essencial ficar atento com a conservação de energia da reação, caso contrário, a simulação não irá funcionar.

Para compilar o programa utilize:

```
bash compilar.sh
```

Esse comando criará um arquivo chamado `arquivobinario` usado para inicializar a simulação, basta utilizar:

```
./arquivobinario vis.mac
```

para rodar o programa. Após finalizado, é necessário abrir o arquivo `G4Data0.heprep` com o programa HepRApp para a visualização do evento. Os dados da simulação são guardados dentro da pasta `/ROOT` no formato `.root`.



