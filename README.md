# Estimativa de dependências gênicas em PDAC a partir do transcriptoma
 
Pipeline reprodutível que estima dependências gênicas em adenocarcinoma ductal
pancreático (PDAC) a partir da expressão gênica, aprendendo a relação entre estado
transcricional e dependência em linhagens celulares do **DepMap** e transferindo-a para
tumores de pacientes do **TCGA-PAAD**.
 
As linhagens celulares têm as duas medidas: expressão e dependência medida por CRISPR.
Os tumores têm apenas expressão. O pipeline cura as duas bases, alinha os perfis
transcricionais, treina um modelo *elastic net* multissaída nas linhagens, avalia-o sem
vazamento de dados e o aplica aos 150 tumores curados da coorte.
 
> Código e resultados do artigo *Predição de dependências gênicas em tumores de adenocarcinoma ductal pancreático a partir de dados públicos de transcriptoma*, submetido ao
> **Journal of Health Informatics**. <!-- Atualizar com DOI após o aceite. -->
 
---
 
## Estrutura do repositório
 
```
depmap-tcga-pdac/
├── depmap-tcga-pdac.Rmd     # relatório único: todas as etapas, da curadoria à aplicação
├── R/
│   ├── tcc_helpers.R        # funções do alinhamento, seleção de alvos e modelos
│   ├── tcc_simulate.R       # geração de dados simulados para testar o código
│   ├── download_tcga.R      # download do TCGA-PAAD via TCGAbiolinks
│   └── render.R             # script de renderização do relatório
├── figuras/                 # figuras geradas pelo relatório (.png)
├── tabelas/                 # tabelas geradas pelo relatório (.csv)
├── data/                    # dados brutos — NÃO versionado, ver "Dados" abaixo
│   ├── depmap/
│   └── tcga/
└── README.md
```
 
As pastas `figuras/` e `tabelas/` já vêm preenchidas com as saídas da execução que
originou os resultados do artigo, de modo que é possível conferir todos os números sem
baixar os dados brutos nem executar o pipeline.
 
---
 
## Requisitos
 
- **R 4.6.2** ou superior
- Pacotes: `TCGAbiolinks` (2.40.0), `glmnet`, `FNN`, `uwot`, `survival`, `survminer`,
  `ggplot2`, `dplyr`, `knitr`, `rmarkdown`
```r
install.packages(c("glmnet", "FNN", "uwot", "survival", "survminer",
                   "ggplot2", "dplyr", "knitr", "rmarkdown"))
 
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
```
 
Não é necessária infraestrutura de alto desempenho nem processamento gráfico. A execução
que gerou os resultados do artigo rodou em computador pessoal com processador de oito
núcleos e 16 GB de memória, e levou cerca de **duas horas**, incluindo as 50 permutações
do nulo.
 
---
 
## Dados
 
Os dados brutos **não estão versionados** neste repositório: somam vários gigabytes e são
distribuídos pelas fontes originais sob seus próprios termos de uso. Ambas as bases são
públicas e de acesso aberto.
 
### DepMap (versão 26Q3)
 
Baixe de <https://depmap.org/portal/download/> e coloque em `data/depmap/`:
 
| Arquivo | Conteúdo |
| --- | --- |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | expressão, log2(TPM+1) |
| `CRISPRGeneEffect.csv` | escores de dependência (Chronos) |
| `Model.csv` | descrição das linhagens |
| `OmicsCNGene.csv` | número de cópias gênicas |
 
O arquivo de número de cópias do 26Q3 é distribuído como `OmicsCNGeneWGS.csv`; renomeie
para `OmicsCNGene.csv`. Versões anteriores do DepMap usam outros nomes de arquivo — ver
os comentários no bloco `params` do Rmd.
 
### TCGA-PAAD
 
Baixado programaticamente pelo `TCGAbiolinks` a partir do Genomic Data Commons, com o
*workflow* **STAR Counts**. Para executar o download, defina `run_download: true` no
bloco `params` do Rmd, ou rode `R/download_tcga.R` diretamente. Os arquivos são gravados
em `data/tcga/` e mantidos em cache local.
 
---
 
## Como executar
 
```r
rmarkdown::render("depmap-tcga-pdac.Rmd")
```
 
ou, equivalentemente, `source("R/render.R")`.
 
O comportamento é controlado pelo bloco `params`, no cabeçalho do Rmd:
 
| Parâmetro | Padrão | Função |
| --- | --- | --- |
| `simulate` | `false` | `true` executa com dados simulados; útil para testar o código sem baixar nada |
| `run_download` | `false` | `true` baixa DepMap e TCGA |
| `use_cache` | `true` | reaproveita os dados já baixados |
| `run_models` | `true` | `false` pula treino e avaliação |
| `run_lpo` | `true` | treino sem linhagens pancreáticas |
| `run_eval` | `true` | reamostragem, permutação e correção de FDR |
 
**Para uma primeira execução sem baixar os dados**, use `simulate: true`. O pipeline roda
de ponta a ponta com dados sintéticos de mesma estrutura, o que permite verificar a
instalação e entender o fluxo em poucos minutos.
 
---
 
## Parâmetros da análise
 
Todos fixados antes da observação dos resultados.
 
| Parâmetro | Valor | Descrição |
| --- | --- | --- |
| `expr_thr` / `expr_frac` | 1 / 0,10 | filtro de expressão: log2(TPM+1) > 1 em ≥ 10% das amostras de cada fonte |
| `n_genes_align` | 5000 | genes de maior variância usados no alinhamento |
| `n_cpc` | 5 | componentes contrastivos removidos |
| `n_pc` | 50 | componentes principais usados como preditores |
| `center_anchor` | `"pdac"` | ancoragem declarada no centroide das linhagens de PDAC |
| `n_genes_y` | 1000 | genes-alvo por rodada, selecionados dentro do *fold* |
| `n_folds` / `inner_k` | 5 / 5 | validação cruzada externa e interna, agrupadas por paciente |
| `en_alpha` | 0,5 | proporção entre lasso e ridge no elastic net |
| `knn_k` | 20 | vizinhos do modelo de referência |
| `n_perm` | 50 | permutações do nulo |
| `n_boot` | 1000 | reamostragens para os intervalos de confiança |
| `r_thr` | 0,3 | limiar de predizibilidade, pré-especificado |
| `mtap_cn_rel_thr` | 0,3 | limiar de deleção de MTAP nas linhagens |
| `seed` | 20261115 | semente aleatória |
 
Um gene é considerado **previsível** quando satisfaz simultaneamente três critérios
definidos antes da análise: r ≥ 0,3, taxa de falsas descobertas < 0,10 e limite inferior
do intervalo de confiança de 95% acima de zero.
 
---
 
## Reprodutibilidade
 
- Versões **fixas** das duas bases (DepMap 26Q3; TCGA-PAAD via STAR Counts), mantidas em
  cache local.
- Semente aleatória única, declarada em `params`.
- Todas as etapas dependentes dos dados — alinhamento, seleção de genes-alvo,
  padronização e escolha de λ — são refeitas **dentro de cada rodada** da validação
  cruzada, usando apenas as linhagens de treino.
- Linhagens derivadas de um mesmo paciente são mantidas na mesma partição.
- O `sessionInfo()` completo é emitido ao final do relatório.
