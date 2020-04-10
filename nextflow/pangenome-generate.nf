#!/usr/bin/env nextflow

params.dir = "${baseDir}/../data"

fastaFiles = "${params.dir}/**.fasta"
sequences = Channel.fromPath(fastaFiles).map { path -> tuple(path.simpleName, path) }

process dedup {
  tag { sample }
  container "quay.io/biocontainers/seqkit:0.7.1--0"

  input:
    set sample, file(fasta) from sequences
  output:
    set sample, file("${sample}.dedup.fasta") into dedup

  """
  seqkit rmdup --by-seq --ignore-case -o ${sample}.dedup.fasta $fasta
  """
}

process overlapReads {
  tag { sample }
  publishDir "$sample", mode: 'copy'
  container "quay.io/biocontainers/minimap2:2.17--h8b12597_1"

  input:
    set sample, file(fasta) from dedup
  output:
    set sample, file(fasta), file("${sample}.paf") into alignments

  """
  minimap2 -cx asm20 -w 1 -t ${task.cpus} $fasta $fasta > ${sample}.paf
  """
}

process induceGraph {
  tag { sample }
  publishDir "$sample", mode: 'copy'
  container "quay.io/biocontainers/seqwish:0.4.1--h8b12597_0"

  input:
    set sample, file(fasta), file(alignment) from alignments
  output:
    set sample, file("${sample}.gfa") into graphs

  """
  seqwish -t ${task.cpus} -k 16 -s $fasta -p $alignment -g ${sample}.gfa
  """
}

(graphsToSort, graphsToVisualize) = graphs.into(2)

process traversePaths {
  tag { sample }
  publishDir "$sample", mode: 'copy'
  container "quay.io/biocontainers/dsh-bio:1.3.3--0"

  input:
    set sample, file(graph) from graphsToVisualize
  output:
    set sample, file("${sample}.withTraversals.gfa") into graphsWithTraversals

  """
  dsh-bio traverse-paths -i $graph | dsh-bio truncate-paths -o ${sample}.withTraversals.gfa
  """
}

process vizGraph {
  tag { sample }
  publishDir "$sample", mode: 'copy'
  container "ubuntu"

  input:
    set sample, file(graph) from graphsWithTraversals
  output:
    set sample, file("${sample}.forceDirected.svg"), file("${sample}.forceDirected.cys"), file("${sample}.hierarchical.svg"), file("${sample}.hierarchical.cys") into visualizations

  """
  touch ${sample}.forceDirected.svg ${sample}.forceDirected.cys ${sample}.hierarchical.svg ${sample}.hierarchical.cys
  """
  /*
  """
  run_cytoscape_assembly_viz.sh -i $graph -o ${sample}.forceDirected.svg ${sample}.forceDirected.cys ${sample}.hierarchical.svg ${sample}.hierarchical.cys
  """
  */
}

process buildSortedGraph {
  tag { sample }
  publishDir "$sample", mode: 'copy'
  container "quay.io/biocontainers/odgi:v0.3--py37h8b12597_0"

  input:
    set sample, file(graph) from graphsToSort
  output:
    set sample, file("${sample}.odgi") into sortedGraphs

  """
  odgi build -g $graph -s -o - |\
     odgi sort -i - -p s -o ${sample}.odgi
  """
}

process vizSortedGraph {
  tag { sample }
  publishDir "$sample", mode: 'copy'
  container "quay.io/biocontainers/odgi:v0.3--py37h8b12597_0"

  input:
    set sample, file(sortedGraph) from sortedGraphs
  output:
    set sample, file("${sample}.png") into sortedVisualizations

  """
  odgi viz -i $sortedGraph -o ${sample}.png -x 50000 -y 500 -R -P 4 -R
  """
}

