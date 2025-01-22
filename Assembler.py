import sys
from collections import defaultdict, deque

# Importa o módulo sys para gerenciar argumentos da linha de comando.
# Importa defaultdict para criar dicionários com valores padrão automaticamente.
# Importa deque, uma estrutura de dados para manipulação eficiente de filas e pilhas.

def read_kmers_from_file(filename):
    """
    Lê os fragmentos de DNA (k-mers) do arquivo de entrada.
    O propósito é processar um arquivo de texto contendo k-mers separados por vírgulas.
    """
    with open(filename, 'r') as file:
        # Abre o arquivo especificado em modo de leitura.
        content = file.read()
        # Lê todo o conteúdo do arquivo como uma única string.
        kmers = content.strip().split(',')
        # Remove espaços em branco no início e fim da string e divide pelos separadores de vírgula.
    return kmers
    # Retorna uma lista contendo os k-mers extraídos do arquivo.

def build_de_bruijn_graph(kmers):
    """
    Constrói o grafo de De Bruijn a partir dos k-mers.
    Utiliza os prefixos e sufixos dos k-mers como nós para criar as conexões.
    """
    graph = defaultdict(list)
    # Cria um dicionário onde cada chave mapeia para uma lista de sufixos.
    in_degree = defaultdict(int)
    # Inicializa os graus de entrada dos nós como zero.
    out_degree = defaultdict(int)
    # Inicializa os graus de saída dos nós como zero.

    for kmer in kmers:
        # Itera sobre cada k-mer na lista fornecida.
        prefix = kmer[:-1]
        # Obtém o prefixo do k-mer (todos os caracteres menos o último).
        suffix = kmer[1:]
        # Obtém o sufixo do k-mer (todos os caracteres menos o primeiro).
        graph[prefix].append(suffix)
        # Adiciona o sufixo como um vizinho do prefixo no grafo.
        out_degree[prefix] += 1
        # Incrementa o grau de saída do nó prefixo.
        in_degree[suffix] += 1
        # Incrementa o grau de entrada do nó sufixo.

    return graph, in_degree, out_degree
    # Retorna o grafo, junto com os graus de entrada e saída.

def find_starting_node(graph, in_degree, out_degree):
    """
    Encontra o nó inicial para o passeio Euleriano com base nos graus de entrada e saída.
    """
    start_node = None
    # Inicializa a variável para armazenar o nó inicial.

    for node in graph:
        # Itera sobre todos os nós do grafo.
        if out_degree[node] - in_degree[node] == 1:
            # Se o nó tem grau de saída 1 maior que o de entrada:
            return node
            # Retorna esse nó imediatamente como o ponto de partida.
        if in_degree[node] == out_degree[node]:
            # Se o nó é equilibrado (grau de entrada igual ao grau de saída):
            start_node = node
            # Salva como uma possível escolha para o ponto de partida.

    return start_node
    # Retorna o nó inicial, ou None se nenhum nó válido for encontrado.

def reconstruct_dna(graph, start_node):
    """
    Reconstrói a sequência de DNA encontrando um caminho Euleriano no grafo.
    Utiliza uma pilha para explorar caminhos e uma fila para construir o caminho final.
    """
    stack = [start_node]
    # Inicializa a pilha com o nó inicial.
    path = deque()
    # Inicializa a fila para armazenar o caminho final.

    while stack:
        # Enquanto houver nós na pilha:
        current = stack[-1]
        # Obtém o nó no topo da pilha sem removê-lo.
        if graph[current]:
            # Se o nó atual tiver vizinhos não visitados:
            next_node = graph[current].pop()
            # Remove e obtém o próximo nó a ser visitado.
            stack.append(next_node)
            # Adiciona o próximo nó à pilha.
        else:
            # Se o nó atual não tiver vizinhos restantes:
            path.appendleft(stack.pop())
            # Remove o nó da pilha e adiciona ao início do caminho.

    return ''.join([path.popleft()] + [node[-1] for node in path])
    # Constrói a sequência de DNA combinando o primeiro nó com o último caractere de cada nó subsequente.

def main(input_file, output_file):
    """
    Função principal que gerencia o fluxo de execução do programa.
    Lê os k-mers, constrói o grafo de De Bruijn, encontra o nó inicial e reconstrói a sequência de DNA.
    """
    kmers = read_kmers_from_file(input_file)
    # Lê os k-mers do arquivo de entrada.
    graph, in_degree, out_degree = build_de_bruijn_graph(kmers)
    # Constrói o grafo de De Bruijn e calcula os graus de entrada e saída.
    start_node = find_starting_node(graph, in_degree, out_degree)
    # Encontra o nó inicial para o caminho Euleriano.

    if start_node is None:
        # Se nenhum nó inicial válido for encontrado:
        raise ValueError("Não foi possível determinar um ponto de partida para reconstrução.")
        # Lança um erro informando que a reconstrução falhou.

    reconstructed_dna = reconstruct_dna(graph, start_node)
    # Reconstrói a sequência de DNA utilizando o caminho Euleriano.

    with open(output_file, 'w') as file:
        # Abre o arquivo de saída em modo de escrita.
        file.write(reconstructed_dna)
        # Escreve a sequência reconstruída no arquivo.

if __name__ == "__main__":
    # Bloco que executa o programa quando o script é chamado diretamente.
    if len(sys.argv) != 2:
        # Verifica se o número de argumentos fornecidos é exatamente 2 (programa + arquivo de entrada).
        print("Uso: python dna_reconstruction.py <arquivo_de_entrada>")
        # Exibe instruções sobre como usar o programa.
        sys.exit(1)
        # Encerra o programa com código de erro.

    input_file = sys.argv[1]
    # Obtém o nome do arquivo de entrada fornecido na linha de comando.
    output_file = "KauaFarias.txt"
    # Define o nome do arquivo de saída onde a sequência será salva.

    main(input_file, output_file)
    # Chama a função principal para executar o programa.
    print(f"Sequência reconstruída salva em {output_file}")
    # Exibe uma mensagem informando que a sequência foi salva com sucesso.

