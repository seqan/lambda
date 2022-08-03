cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

## Index input files

declare_datasource (FILE db_nucl.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/input_files/db_nucl.fasta.gz
                    URL_HASH SHA256=614f8d7863c40facb7fffe666ce04341368a43ea90b8125cb675907f315bb0a2)

declare_datasource (FILE db_nucl_bs.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/input_files/db_nucl_bs.fasta.gz
                    URL_HASH SHA256=160375ac5ff4426f1768981a0215495efd6de0b847e3226c5c7c112370a599a1)

declare_datasource (FILE db_prot.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/input_files/db_prot.fasta.gz
                    URL_HASH SHA256=2c27f09f77e1f8ec0fea0aa1cae75f6634bab852a843a48033fc0f5251d57626)

## Index output files

declare_datasource (FILE db_nucl_fm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/db_nucl_fm.fasta.gz.lba
                    URL_HASH SHA256=8ad7fcaac2635b2c2aa9125e4ace062a0938a430828bb9de756bd4496decd75f)

declare_datasource (FILE db_nucl_bifm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/db_nucl_bifm.fasta.gz.lba
                    URL_HASH SHA256=b88be986b72c8eb137f191be30d6789aaab5d0126aae83029d537b130607e308)

declare_datasource (FILE db_nucl_bs_fm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/db_nucl_bs_fm.fasta.gz.lba
                    URL_HASH SHA256=12260d33026f2a90a1bd0a0e9db82aee306f931471b35b303c9b0caf7e28abf7)

declare_datasource (FILE db_nucl_bs_bifm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/db_nucl_bs_bifm.fasta.gz.lba
                    URL_HASH SHA256=a4fd46461db05f3688f3be332014ee337784a1a19c24ea92c4c9ab8c9f4720ca)

declare_datasource (FILE db_prot_fm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/db_prot_fm.fasta.gz.lba
                    URL_HASH SHA256=bf79914640873832e58bd68e8f4554cf3a44632fed6857517ccc7159c559cf02)

declare_datasource (FILE db_prot_bifm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/db_prot_bifm.fasta.gz.lba
                    URL_HASH SHA256=f60cb23bf2da10a705cd6cb3bd927a44a5acd0d86a1749fc820c396fbdf777b8)

declare_datasource (FILE db_trans_fm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/db_trans_fm.fasta.gz.lba
                    URL_HASH SHA256=1bc68a11e8bfaebbea758cd4d8712d6c361b5b03cea548d0bade9cf6609ddb9e)

declare_datasource (FILE db_trans_bifm.fasta.gz.lba
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/db_trans_bifm.fasta.gz.lba
                    URL_HASH SHA256=e56d8ba576ad313df4f7bae1cd2193e3c11699749a985ecc1807af65b9815fdc)

## Query input files

declare_datasource (FILE queries_nucl.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/input_files/queries_nucl.fasta.gz
                    URL_HASH SHA256=7a8adcc7ee5d967992a0624b9fb704118223abd112534c2c6456104b30ddad88)

declare_datasource (FILE queries_nucl_bs.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/input_files/queries_nucl_bs.fasta.gz
                    URL_HASH SHA256=a358ace3e6f35bd379854fd8870a63a6da62b6bd3e72d95f11e5386398881f0a)

declare_datasource (FILE queries_prot.fasta.gz
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/input_files/queries_prot.fasta.gz
                    URL_HASH SHA256=e21411f422c1dd844696c8ca8a08b93e8243d7f77307731c5f2ee8dcc60d0703)

## Search output files

declare_datasource (FILE output_blastn_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_fm.bam
                    URL_HASH SHA256=041916249e2ece6e3faf2b64dcadc8e4544a3c80ec2b1ef4e64d4435eedf478c)

declare_datasource (FILE output_blastn_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_fm.m0
                    URL_HASH SHA256=4e8034135def8cc9935e5c1bf615c9cf2e74dee59fe5a1a741bef3d7af11b449)

declare_datasource (FILE output_blastn_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_fm.m8
                    URL_HASH SHA256=7b0d76cd928e434389e185ba847ee716dd10e3403831b9cac4a85acfe49bfcda)

declare_datasource (FILE output_blastn_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_fm.m9
                    URL_HASH SHA256=6e143013d9185db57d6202edbd50f93748fd470c698c0bd583929302dfec06e2)

declare_datasource (FILE output_blastn_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_fm.sam
                    URL_HASH SHA256=b3811d72bbe57b7b1589adccfb7ee8e6707c2bc40c2f5f4fef920ce38980c8e9)

declare_datasource (FILE output_blastn_bs_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_bs_fm.bam
                    URL_HASH SHA256=d0fe045258c6210ce05a6a14d73fa8c0e115cb8e49f91a761a74251190e580ad)

declare_datasource (FILE output_blastn_bs_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_bs_fm.m0
                    URL_HASH SHA256=de5e4d989a43d9c25c1a135e7ca538c4c5d936b039e0d6e6650ae6e04a792da7)

declare_datasource (FILE output_blastn_bs_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_bs_fm.m8
                    URL_HASH SHA256=1fda67245d7659c1aca4952931904143914fd65a3f1f44ef7c55ed1a5fc038c6)

declare_datasource (FILE output_blastn_bs_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_bs_fm.m9
                    URL_HASH SHA256=d9f5e0aa346c755dae29cbebc72e2c8b8daab42e2387496360eb0ac503a9f822)

declare_datasource (FILE output_blastn_bs_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_bs_fm.sam
                    URL_HASH SHA256=a5f56977390f6faf4d8b147faa7ddcaad8bdd5544c943acbf0eeabe83b7181e5)

declare_datasource (FILE output_blastp_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastp_fm.bam
                    URL_HASH SHA256=a3c464f2f8151087e59c1cad119a2fcd04e34bd8c6bc91dec3bba61861eeb196)

declare_datasource (FILE output_blastp_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastp_fm.m0
                    URL_HASH SHA256=9299b623245baac372cb6c7c0a7ab966d9bc3d5f01578600912d9404a80e0280)

declare_datasource (FILE output_blastp_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastp_fm.m8
                    URL_HASH SHA256=527367af7a6a6482be01005176e78078498ec37398f945ea6ed32726b6761b10)

declare_datasource (FILE output_blastp_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastp_fm.m9
                    URL_HASH SHA256=e039aaca5c9564c8a79b5f8fcab98b4b38c521d737384034044eabe67010ae02)

declare_datasource (FILE output_blastp_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastp_fm.sam
                    URL_HASH SHA256=3835c081d6f6892b4619657feadbfb9974be15bbf8b639fbd37d6dfbcd06e78b)

declare_datasource (FILE output_blastx_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastx_fm.bam
                    URL_HASH SHA256=85b7370608b4df7ee57a26f6c6af68dce435f767c42d6d45e8064836792639c3)

declare_datasource (FILE output_blastx_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastx_fm.m0
                    URL_HASH SHA256=e7c2599620e18ce7c75e079e9357f75bd9928dcb273b6c860f6fe6c9c4b7cd87)

declare_datasource (FILE output_blastx_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastx_fm.m8
                    URL_HASH SHA256=d0558f47c6cc1b9682c49ff6eb10025f2eee3cc3573989e7e455e8686bfdbbf9)

declare_datasource (FILE output_blastx_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastx_fm.m9
                    URL_HASH SHA256=ce317b439daed2db0514335bc5477b158b62bfe97d0af7038c12ef7a3105bc52)

declare_datasource (FILE output_blastx_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastx_fm.sam
                    URL_HASH SHA256=acfc8b63a6e4c5547726caafd2a8c65b0e04f7bd67c64a7b724a6d7315637e5b)

declare_datasource (FILE output_tblastn_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastn_fm.bam
                    URL_HASH SHA256=25550f19f6d1add842ea7df63b308c3fa41a9261e794f4042c7cdaf194611395)

declare_datasource (FILE output_tblastn_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastn_fm.m0
                    URL_HASH SHA256=9d8bd728356f1627602021cd2c3ab70ae5c027584a10b7b23bb25cba549461fb)

declare_datasource (FILE output_tblastn_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastn_fm.m8
                    URL_HASH SHA256=82eb91515a223a992511dc5e88075fa3c94495464dc1d972eacc5afb9041a3c2)

declare_datasource (FILE output_tblastn_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastn_fm.m9
                    URL_HASH SHA256=0e29e58218cf7d154fb8ab865a04cd81bf8302c9337387254758f04f1e9ab5ab)

declare_datasource (FILE output_tblastn_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastn_fm.sam
                    URL_HASH SHA256=af55fb74957ecf1816bcdbcbc53777129193c0a4ac40a463327d777450d5f86c)

declare_datasource (FILE output_tblastx_fm.bam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastx_fm.bam
                    URL_HASH SHA256=69a3b82cc5ef84999cfcb560bf4122010f965300c7c157087fa2ea9909880547)

declare_datasource (FILE output_tblastx_fm.m0
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastx_fm.m0
                    URL_HASH SHA256=25de5ba2fef1a53582403381666179d7cac9f6f8d0e5299ff217cf879a93f7ed)

declare_datasource (FILE output_tblastx_fm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastx_fm.m8
                    URL_HASH SHA256=b15751f7f4686d9d3562a14ade26241046da046e6f341cf631ebbe1a9e0d27f4)

declare_datasource (FILE output_tblastx_fm.m9
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastx_fm.m9
                    URL_HASH SHA256=79981091202983e3d07c26876fc8bd6fe3cd214187db0d91606df17aade32577)

declare_datasource (FILE output_tblastx_fm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastx_fm.sam
                    URL_HASH SHA256=15b08fcc6186c99e8b5680a19708f9aba4f04bf0d5d5becac31ac501eaf6bd4e)

declare_datasource (FILE output_blastn_bifm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_bifm.sam
                    URL_HASH SHA256=b3811d72bbe57b7b1589adccfb7ee8e6707c2bc40c2f5f4fef920ce38980c8e9)

declare_datasource (FILE output_blastn_bifm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_bifm.m8
                    URL_HASH SHA256=7b0d76cd928e434389e185ba847ee716dd10e3403831b9cac4a85acfe49bfcda)

declare_datasource (FILE output_blastn_bs_bifm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_bs_bifm.sam
                    URL_HASH SHA256=a5f56977390f6faf4d8b147faa7ddcaad8bdd5544c943acbf0eeabe83b7181e5)

declare_datasource (FILE output_blastn_bs_bifm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastn_bs_bifm.m8
                    URL_HASH SHA256=1fda67245d7659c1aca4952931904143914fd65a3f1f44ef7c55ed1a5fc038c6)

declare_datasource (FILE output_blastp_bifm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastp_bifm.sam
                    URL_HASH SHA256=3835c081d6f6892b4619657feadbfb9974be15bbf8b639fbd37d6dfbcd06e78b)

declare_datasource (FILE output_blastp_bifm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastp_bifm.m8
                    URL_HASH SHA256=527367af7a6a6482be01005176e78078498ec37398f945ea6ed32726b6761b10)

declare_datasource (FILE output_blastx_bifm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastx_bifm.sam
                    URL_HASH SHA256=acfc8b63a6e4c5547726caafd2a8c65b0e04f7bd67c64a7b724a6d7315637e5b)

declare_datasource (FILE output_blastx_bifm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_blastx_bifm.m8
                    URL_HASH SHA256=d0558f47c6cc1b9682c49ff6eb10025f2eee3cc3573989e7e455e8686bfdbbf9)

declare_datasource (FILE output_tblastn_bifm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastn_bifm.sam
                    URL_HASH SHA256=af55fb74957ecf1816bcdbcbc53777129193c0a4ac40a463327d777450d5f86c)

declare_datasource (FILE output_tblastn_bifm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastn_bifm.m8
                    URL_HASH SHA256=82eb91515a223a992511dc5e88075fa3c94495464dc1d972eacc5afb9041a3c2)

declare_datasource (FILE output_tblastx_bifm.sam
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastx_bifm.sam
                    URL_HASH SHA256=15b08fcc6186c99e8b5680a19708f9aba4f04bf0d5d5becac31ac501eaf6bd4e)

declare_datasource (FILE output_tblastx_bifm.m8
                    URL https://raw.githubusercontent.com/h-2/lambda-testdata/main/output_files/output_tblastx_bifm.m8
                    URL_HASH SHA256=b15751f7f4686d9d3562a14ade26241046da046e6f341cf631ebbe1a9e0d27f4)
