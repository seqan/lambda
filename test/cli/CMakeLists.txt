cmake_minimum_required (VERSION 3.8)

add_cli_test (index_test.cpp)
target_use_datasources (index_test FILES db_nucl.fasta.gz
                                         db_nucl_bs.fasta.gz
                                         db_prot.fasta.gz
                                         db_nucl_bs_bifm.fasta.gz.lba
                                         db_nucl_bs_fm.fasta.gz.lba
                                         db_nucl_bifm.fasta.gz.lba
                                         db_nucl_fm.fasta.gz.lba
                                         db_prot_bifm.fasta.gz.lba
                                         db_prot_fm.fasta.gz.lba
                                         db_trans_bifm.fasta.gz.lba
                                         db_trans_fm.fasta.gz.lba)

add_cli_test (search_test.cpp)
target_use_datasources (search_test FILES db_nucl.fasta.gz
                                          db_nucl_bs.fasta.gz
                                          db_prot.fasta.gz
                                          queries_nucl.fasta.gz
                                          queries_nucl_bs.fasta.gz
                                          queries_prot.fasta.gz
                                          output_blastn_bs_fm.bam
                                          output_blastn_bs_fm.m0
                                          output_blastn_bs_fm.m8
                                          output_blastn_bs_fm.m9
                                          output_blastn_bs_fm.sam
                                          output_blastn_fm.bam
                                          output_blastn_fm.m0
                                          output_blastn_fm.m8
                                          output_blastn_fm.m9
                                          output_blastn_fm.sam
                                          output_blastp_fm.bam
                                          output_blastp_fm.m0
                                          output_blastp_fm.m8
                                          output_blastp_fm.m9
                                          output_blastp_fm.sam
                                          output_blastx_fm.bam
                                          output_blastx_fm.m0
                                          output_blastx_fm.m8
                                          output_blastx_fm.m9
                                          output_blastx_fm.sam
                                          output_tblastn_fm.bam
                                          output_tblastn_fm.m0
                                          output_tblastn_fm.m8
                                          output_tblastn_fm.m9
                                          output_tblastn_fm.sam
                                          output_tblastx_fm.bam
                                          output_tblastx_fm.m0
                                          output_tblastx_fm.m8
                                          output_tblastx_fm.m9
                                          output_tblastx_fm.sam
                                          output_blastn_bs_fm_fast.m8
                                          output_blastn_bs_fm_fast.sam
                                          output_blastn_fm_fast.m8
                                          output_blastn_fm_fast.sam
                                          output_blastp_fm_fast.m8
                                          output_blastp_fm_fast.sam
                                          output_blastx_fm_fast.m8
                                          output_blastx_fm_fast.sam
                                          output_tblastn_fm_fast.m8
                                          output_tblastn_fm_fast.sam
                                          output_tblastx_fm_fast.m8
                                          output_tblastx_fm_fast.sam
                                          output_blastn_bs_fm_sensitive.m8
                                          output_blastn_bs_fm_sensitive.sam
                                          output_blastn_fm_sensitive.m8
                                          output_blastn_fm_sensitive.sam
                                          output_blastp_fm_sensitive.m8
                                          output_blastp_fm_sensitive.sam
                                          output_blastx_fm_sensitive.m8
                                          output_blastx_fm_sensitive.sam
                                          output_tblastn_fm_sensitive.m8
                                          output_tblastn_fm_sensitive.sam
                                          output_tblastx_fm_sensitive.m8
                                          output_tblastx_fm_sensitive.sam)
