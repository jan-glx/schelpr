#source: https://github.com/lonsbio/fastqe
fastq_emoji_map_binned= c(
  '!'= 'no_entry_sign',
  '"'= 'no_entry_sign',

  #2–9 6
  '#'= 'skull',
  '$'= 'skull',
  '%'= 'skull',
  '&'= 'skull',
  '\''= 'skull',
  '('= 'skull',
  ')'= 'skull',
  '*'= 'skull',

  #10–19 15
  '+'= 'poop' ,
  ','= 'poop' ,
  '-'= 'poop' ,
  '.'= 'poop' ,
  '/'= 'poop' ,
  '0'= 'poop' ,
  '1'= 'poop' ,
  '2'= 'poop' ,
  '3'= 'poop' ,
  '4'= 'poop' ,

  #20–24 22
  '5'= 'warning',
  '6'= 'warning',
  '7'= 'warning',
  '8'= 'warning',
  '9'= 'warning',


  #25–29 27
  ':'= 'smile',
  ';'= 'smile',
  '<'= 'smile',
  '='= 'smile',
  '>'= 'smile',


  #30–34 33
  '?'= 'laughing',
  '@'= 'laughing',
  'A'= 'laughing',
  'B'= 'laughing',
  'C'= 'laughing',

  #35–39 37
  'D'= 'sunglasses',
  'E'= 'sunglasses',
  'F'= 'sunglasses',
  'G'= 'sunglasses',
  'H'= 'sunglasses',

  #≥ 40 40
  'I'= 'heart_eyes',
  'J'= 'heart_eyes'
)
fastq_emoji_map_binned <- sapply(fastq_emoji_map_binned, emo::ji)

#' @export
fastq2emoji <- function(x) do.call(paste0,lapply(tstrsplit(x, split=""), function(x) fastq_emoji_map_binned[x]))
