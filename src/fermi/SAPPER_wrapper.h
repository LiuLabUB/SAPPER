#ifdef __cplusplus
extern "C" {
#endif
  
// wrapper to assemble a fastq file and output anther fastq, given the overlap option.
  void assemble( char* infq, char* outfq, int unitig_k, int min_merge_len);

#ifdef __cplusplus
}
#endif
