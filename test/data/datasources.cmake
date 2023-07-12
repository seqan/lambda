cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

## Index input files

set(BASEURL "https://raw.githubusercontent.com/h-2/lambda-testdata/21e2a230deec072ddc6da4b67403ce684dddf2c2")

declare_datasource (FILE db_nucl.fasta.gz
                    URL ${BASEURL}/input_files/db_nucl.fasta.gz
                    URL_HASH SHA256=614f8d7863c40facb7fffe666ce04341368a43ea90b8125cb675907f315bb0a2)

declare_datasource (FILE db_nucl_bs.fasta.gz
                    URL ${BASEURL}/input_files/db_nucl_bs.fasta.gz
                    URL_HASH SHA256=160375ac5ff4426f1768981a0215495efd6de0b847e3226c5c7c112370a599a1)

declare_datasource (FILE db_prot.fasta.gz
                    URL ${BASEURL}/input_files/db_prot.fasta.gz
                    URL_HASH SHA256=2c27f09f77e1f8ec0fea0aa1cae75f6634bab852a843a48033fc0f5251d57626)

## Index output files

declare_datasource (FILE db_nucl_fm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_nucl_fm.fasta.gz.lba
                    URL_HASH SHA256=68de7b2822b3f47268efee8becd400697cb1581a7055d1b581af365dd27f7a89)

declare_datasource (FILE db_nucl_bifm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_nucl_bifm.fasta.gz.lba
                    URL_HASH SHA256=e812f75f1672b5dc4b8cb6078ce209ec809d21783424836ccbe0c9b78f37887c)

declare_datasource (FILE db_nucl_bs_fm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_nucl_bs_fm.fasta.gz.lba
                    URL_HASH SHA256=b17d4e29f424af54fa802344ae9d4af974dbfc16d0445741e08de6f3cfface0b)

declare_datasource (FILE db_nucl_bs_bifm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_nucl_bs_bifm.fasta.gz.lba
                    URL_HASH SHA256=6dc216251505f0a913dffb9dc97699898492d511934dbcd1c9961337a151a62b)

declare_datasource (FILE db_prot_fm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_prot_fm.fasta.gz.lba
                    URL_HASH SHA256=57fa7a4319c5eb4dded850a79beefa1c5eca6927c91541cd7d658fde0ab7d166)

declare_datasource (FILE db_prot_bifm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_prot_bifm.fasta.gz.lba
                    URL_HASH SHA256=76d0405235aff8145fe0d90399ffb7d3619ad2fedc37776771a1ca8a38b79a5d)

declare_datasource (FILE db_trans_fm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_trans_fm.fasta.gz.lba
                    URL_HASH SHA256=dc94a97d19a2ad224b995640bcacfc8fb6c7aa87ef134f3e61647c5b31286696)

declare_datasource (FILE db_trans_bifm.fasta.gz.lba
                    URL ${BASEURL}/output_files/db_trans_bifm.fasta.gz.lba
                    URL_HASH SHA256=3ffce0f1115888bce562460e7f4f1222af7052cd8a85a32db0813eac317d8be5)

## Query input files

declare_datasource (FILE queries_nucl.fasta.gz
                    URL ${BASEURL}/input_files/queries_nucl.fasta.gz
                    URL_HASH SHA256=7a8adcc7ee5d967992a0624b9fb704118223abd112534c2c6456104b30ddad88)

declare_datasource (FILE queries_nucl_bs.fasta.gz
                    URL ${BASEURL}/input_files/queries_nucl_bs.fasta.gz
                    URL_HASH SHA256=a358ace3e6f35bd379854fd8870a63a6da62b6bd3e72d95f11e5386398881f0a)

declare_datasource (FILE queries_prot.fasta.gz
                    URL ${BASEURL}/input_files/queries_prot.fasta.gz
                    URL_HASH SHA256=e21411f422c1dd844696c8ca8a08b93e8243d7f77307731c5f2ee8dcc60d0703)

## Search output files

declare_datasource (FILE output_blastn_fm.bam
                    URL ${BASEURL}/output_files/output_blastn_fm.bam
                    URL_HASH SHA256=f181cd5e083df6879044153eced4b545a6c392e3c3e55f733799c37f32be16fb)

declare_datasource (FILE output_blastn_fm.m0
                    URL ${BASEURL}/output_files/output_blastn_fm.m0
                    URL_HASH SHA256=8d4083e6c4f0979497df61fec5169c1e82a9111bf052b77980c987a5012acddf)

declare_datasource (FILE output_blastn_fm.m8
                    URL ${BASEURL}/output_files/output_blastn_fm.m8
                    URL_HASH SHA256=5ee0f2c00efa90aa564b4f3a1a05e5bce084421612d1776d4431377e296c9f5f)

declare_datasource (FILE output_blastn_fm.m9
                    URL ${BASEURL}/output_files/output_blastn_fm.m9
                    URL_HASH SHA256=46090b822625b64628c0b4bbb7c068a1224a77d5e6c6dd7f5dc4901943e94096)

declare_datasource (FILE output_blastn_fm.sam
                    URL ${BASEURL}/output_files/output_blastn_fm.sam
                    URL_HASH SHA256=a1de463d54ee72ff24cade4526b340aab9e48d6232e3e98b36b9f6560260ef0b)

declare_datasource (FILE output_blastn_bs_fm.bam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.bam
                    URL_HASH SHA256=8a73e5bca4b7a3225114170906c520c74bd066c2540898369a9c6ab24f517f5e)

declare_datasource (FILE output_blastn_bs_fm.m0
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.m0
                    URL_HASH SHA256=5be6f6e2e5e57a9a07d364d1fd5c376e13db9aa4093f1445cdca17ee7e143790)

declare_datasource (FILE output_blastn_bs_fm.m8
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.m8
                    URL_HASH SHA256=a03ec852936061867791a4dabffffe8cead43b24c0e0afba398e134841fc193d)

declare_datasource (FILE output_blastn_bs_fm.m9
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.m9
                    URL_HASH SHA256=2b71a95f02f5cc59037d6d1b34b4fd6141e2d36a0e9a8803a6d6c5c9e98482de)

declare_datasource (FILE output_blastn_bs_fm.sam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm.sam
                    URL_HASH SHA256=3fe84a825c55004a87b34ee24312d64961a08737e91a948bfe0b5b3bcb2f3e38)

declare_datasource (FILE output_blastp_fm.bam
                    URL ${BASEURL}/output_files/output_blastp_fm.bam
                    URL_HASH SHA256=322932054405c802211f9fd1f9e86a96c8b19ea85d2eaa8cbd18f4bcfb1d5131)

declare_datasource (FILE output_blastp_fm.m0
                    URL ${BASEURL}/output_files/output_blastp_fm.m0
                    URL_HASH SHA256=56e29329f35f1868d6c5a2ab0c88e7ea1446d433f5d7bb70308ceb97599fcd5f)

declare_datasource (FILE output_blastp_fm.m8
                    URL ${BASEURL}/output_files/output_blastp_fm.m8
                    URL_HASH SHA256=99f520bb55f5c1b371ae5ba41b881c21360a84e7fcb8ff097e01036beb801d4c)

declare_datasource (FILE output_blastp_fm.m9
                    URL ${BASEURL}/output_files/output_blastp_fm.m9
                    URL_HASH SHA256=3d7b519b26419361c37f5512b28220f3004620aafa2aae62f0e13bbee3ee92f2)

declare_datasource (FILE output_blastp_fm.sam
                    URL ${BASEURL}/output_files/output_blastp_fm.sam
                    URL_HASH SHA256=d7267490e2c6c1ebf838fb7f08a5f370bf0f2b23024ca76db6bfed440851fc91)

declare_datasource (FILE output_blastx_fm.bam
                    URL ${BASEURL}/output_files/output_blastx_fm.bam
                    URL_HASH SHA256=490b13b75528c87530a3f3028bb825ddb96b46bcbe6dbea76ee782ec17abb5a5)

declare_datasource (FILE output_blastx_fm.m0
                    URL ${BASEURL}/output_files/output_blastx_fm.m0
                    URL_HASH SHA256=08b1500f717cb596bd1be7919cd72b6ac2c69caccae2716072c046c45122fbdc)

declare_datasource (FILE output_blastx_fm.m8
                    URL ${BASEURL}/output_files/output_blastx_fm.m8
                    URL_HASH SHA256=f6345d4856a143b9a425a1cca75889e49d2dbc834a817b8ef261d02ac33d2d0d)

declare_datasource (FILE output_blastx_fm.m9
                    URL ${BASEURL}/output_files/output_blastx_fm.m9
                    URL_HASH SHA256=4965194ec39831d4b9343eebc49cbc0871b0cee3fb7219c5626846477dc0307c)

declare_datasource (FILE output_blastx_fm.sam
                    URL ${BASEURL}/output_files/output_blastx_fm.sam
                    URL_HASH SHA256=2447e15abfcfca8b7b91ec788e46b96b363d9d7b125a95471898956da58fa80c)

declare_datasource (FILE output_tblastn_fm.bam
                    URL ${BASEURL}/output_files/output_tblastn_fm.bam
                    URL_HASH SHA256=0d646463d379eeadca4af6705081854034d2e07bce870dee1165f11dc50f9336)

declare_datasource (FILE output_tblastn_fm.m0
                    URL ${BASEURL}/output_files/output_tblastn_fm.m0
                    URL_HASH SHA256=d3922d4db0245d3c797b4ca5a982e73c3331be7b7f3282754e31daf2820691aa)

declare_datasource (FILE output_tblastn_fm.m8
                    URL ${BASEURL}/output_files/output_tblastn_fm.m8
                    URL_HASH SHA256=316ff6d392d1553bf3115f41e14e5cfced72f6b22459e487b0def80c37123b44)

declare_datasource (FILE output_tblastn_fm.m9
                    URL ${BASEURL}/output_files/output_tblastn_fm.m9
                    URL_HASH SHA256=720895bca93a35543fc58f2dbc490ba9ce531195fe71d0dbbc0ed41e6392002f)

declare_datasource (FILE output_tblastn_fm.sam
                    URL ${BASEURL}/output_files/output_tblastn_fm.sam
                    URL_HASH SHA256=b2ca3d794d03154ed956d095e531858d71aac8d81b4fdee5ce18ed585300e034)

declare_datasource (FILE output_tblastx_fm.bam
                    URL ${BASEURL}/output_files/output_tblastx_fm.bam
                    URL_HASH SHA256=ab3d60dc3beaa6d5ac78aaa7c6bd212297c9e7ab982fc9c9f6d3c69c6b819b87)

declare_datasource (FILE output_tblastx_fm.m0
                    URL ${BASEURL}/output_files/output_tblastx_fm.m0
                    URL_HASH SHA256=fcb43c8084305372ebaf931736bc1dc61a1589021aa3d36ff1d0610e7a423cbb)

declare_datasource (FILE output_tblastx_fm.m8
                    URL ${BASEURL}/output_files/output_tblastx_fm.m8
                    URL_HASH SHA256=d040e8debb4e5cff6443e676e77449a65f4ffdb83924660f339d640b121964f1)

declare_datasource (FILE output_tblastx_fm.m9
                    URL ${BASEURL}/output_files/output_tblastx_fm.m9
                    URL_HASH SHA256=970a2a628e1816243645e6eabc26ee8f20f2f2ae447e5d18b432fb125d25e867)

declare_datasource (FILE output_tblastx_fm.sam
                    URL ${BASEURL}/output_files/output_tblastx_fm.sam
                    URL_HASH SHA256=6f218b39e8785e1bebbe8be486237e03e06d9fccd41a5e57e575539cfe4eab6e)

declare_datasource (FILE output_blastn_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastn_fm_fast.sam
                    URL_HASH SHA256=a8121e6138ab3c1d5c0e6b32a38782f07806ff4806a46b706444ebf8cdb318b8)

declare_datasource (FILE output_blastn_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastn_fm_fast.m8
                    URL_HASH SHA256=b6dd323c73944d78367cf4f8b6f22f292d89876c4e6a6020fde7461a5115cfad)

declare_datasource (FILE output_blastn_bs_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_fast.sam
                    URL_HASH SHA256=3fe84a825c55004a87b34ee24312d64961a08737e91a948bfe0b5b3bcb2f3e38)

declare_datasource (FILE output_blastn_bs_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_fast.m8
                    URL_HASH SHA256=a03ec852936061867791a4dabffffe8cead43b24c0e0afba398e134841fc193d)

declare_datasource (FILE output_blastp_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastp_fm_fast.sam
                    URL_HASH SHA256=36bca3f2ee56d2234ed44f8351ad83b5fe267c9e75a9eea5c13d34565f20f2ac)

declare_datasource (FILE output_blastp_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastp_fm_fast.m8
                    URL_HASH SHA256=753acc7afccbcacc8e585e7a83b7b295548ea536b23800ca90a9566bcb7a48cb)

declare_datasource (FILE output_blastx_fm_fast.sam
                    URL ${BASEURL}/output_files/output_blastx_fm_fast.sam
                    URL_HASH SHA256=557394a3c75c2da49d3d15c34f10f1eca2e4a934102cfbe92081df0f704b1828)

declare_datasource (FILE output_blastx_fm_fast.m8
                    URL ${BASEURL}/output_files/output_blastx_fm_fast.m8
                    URL_HASH SHA256=5554e8dcd4de9ef02818e18ad8fb185594762a425aded7929a0b2a574137f46e)

declare_datasource (FILE output_tblastn_fm_fast.sam
                    URL ${BASEURL}/output_files/output_tblastn_fm_fast.sam
                    URL_HASH SHA256=e416ddf4d79ea72f8240d291c17c5e41fc2a2d7b4560f1c9ad600b498f9e4d1a)

declare_datasource (FILE output_tblastn_fm_fast.m8
                    URL ${BASEURL}/output_files/output_tblastn_fm_fast.m8
                    URL_HASH SHA256=e7e34710b246347dfdc12619b064da1ec0f1d248316901209d89672e92d9bbb2)

declare_datasource (FILE output_tblastx_fm_fast.sam
                    URL ${BASEURL}/output_files/output_tblastx_fm_fast.sam
                    URL_HASH SHA256=3c23be9c52c1f05acef72acbfc2ea316e18a9e3fa9a9ad1f59903118783c4129)

declare_datasource (FILE output_tblastx_fm_fast.m8
                    URL ${BASEURL}/output_files/output_tblastx_fm_fast.m8
                    URL_HASH SHA256=dea0901a96f8dafd9e3f7518223257ac6ad0b8e4bf2ec42fe05701707ed72b4a)

declare_datasource (FILE output_blastn_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastn_fm_sensitive.sam
                    URL_HASH SHA256=40042463f6a551e27cf2820650be691f39417bce1d6f24782a1718ff9cbfc6c5)

declare_datasource (FILE output_blastn_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastn_fm_sensitive.m8
                    URL_HASH SHA256=1656959b59c448c919c3e03398eece936fb196f141e1ffbc21185859a144ebd7)

declare_datasource (FILE output_blastn_bs_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_sensitive.sam
                    URL_HASH SHA256=5f4bb4f9b8f43563856c04998069972e01882b694151bc235ead384f97aa3780)

declare_datasource (FILE output_blastn_bs_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastn_bs_fm_sensitive.m8
                    URL_HASH SHA256=947860f9fb506eb0626de5c18bf47b6ba8efaa75fd765a5adc66b95df7958cfd)

declare_datasource (FILE output_blastp_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastp_fm_sensitive.sam
                    URL_HASH SHA256=83fc52cfe429c8ac4d790f69f32856f1426ed003c128d6dc4c8767e977e5e9c1)

declare_datasource (FILE output_blastp_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastp_fm_sensitive.m8
                    URL_HASH SHA256=25b8fa3de4eda51b016baa90ade64f5ff4f409b37b32b16226e46ca350673823)

declare_datasource (FILE output_blastx_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_blastx_fm_sensitive.sam
                    URL_HASH SHA256=d7bf59418964c5b7c5e7d10a173d07d2b4fdaaaa7d684a1a2df1c323fe8620c6)

declare_datasource (FILE output_blastx_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_blastx_fm_sensitive.m8
                    URL_HASH SHA256=039b2929a78ea1c4768a215dd7c6f060bea8c11feab6831d767ca4d7d9cdc0fa)

declare_datasource (FILE output_tblastn_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_tblastn_fm_sensitive.sam
                    URL_HASH SHA256=b9b54e58f158746e1d69fc24d841ba637c859849bc9fe50e2c75cf3bc634a3e8)

declare_datasource (FILE output_tblastn_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_tblastn_fm_sensitive.m8
                    URL_HASH SHA256=c07e11add4b48d1d910c76f95a007de452c16a8fe8ea83f695a1c75ee4152db1)

declare_datasource (FILE output_tblastx_fm_sensitive.sam
                    URL ${BASEURL}/output_files/output_tblastx_fm_sensitive.sam
                    URL_HASH SHA256=8b400bf75ba657dd29983fcf0e5f0f8de2243d8b68287af8b884aab61524ba2e)

declare_datasource (FILE output_tblastx_fm_sensitive.m8
                    URL ${BASEURL}/output_files/output_tblastx_fm_sensitive.m8
                    URL_HASH SHA256=a590cff9711dbbad654a7317c081af5ffacda9c934ece3bf7b59eec223c2cf1f)
