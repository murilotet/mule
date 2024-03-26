# mule
Encoder:
1) para compilar: navegar até o diretório Encoder/PPM/ e digitar make
2) Para codificar um lightfield: encoder -lf <nome do diretorio do lightfield> -o <nome do arquivo comprimido>
-nh <numero de vistas na direcao horizontal> -nv <numero de vistas na direcao vertical>
-off_h <numero da primeira vista na direcao horizontal> -off_v <numero da primeira vista na direcao vertical>
-t <tamanho do bloco 4D na direcao t> -s <tamanho do bloco 4D na direcao s>
-v <tamanho do bloco 4D na direcao v> -u <tamanho do bloco 4D na direcao u>
-min_t <tamanho minimo do bloco 4D na direcao t> -min_s <tamanho minimo do bloco 4D na direcao s>
-min_v <tamanho minimo do bloco 4D na direcao v> -min_u <tamanho minimo do bloco 4D na direcao u>
-lambda <valor do multiplicador de lagrange>

Parâmetros opcionais: -lenslet13x13 para ser usado com light fields lenslets,
-u_scale <ganho> multiplicador float opcional para controlar estouro de aritmética no caso de blocos muito grandes
(desnecessário para blocos de até 9x9x512x512)

Decoder:
1) para compilar: navegar até o diretório Decoder/PPM/ e digitar make
2) Para codificar um lightfield: decoder -lf <nome do diretorio do lightfield> -i <nome do arquivo comprimido>
-nh <numero de vistas na direcao horizontal> -nv <numero de vistas na direcao vertical>
-off_h <numero da primeira vista na direcao horizontal> -off_v <numero da primeira vista na direcao vertical>

Opcionalmente todos os parâmetros podem ser colocados em um arquivo de configuração. Nesse caso,
por exemplo para codificar: encoder -cf <nome do arquivo de configuração>
