#include <fstream>
#include <string>

#include "cli_test.hpp"

struct index_test : public cli_test
{
    std::string md5_cmd;

    void determine_md5_cmd()
    {
        if (std::system("echo foo | md5sum > /dev/null 2>&1") == 0)
            md5_cmd = "md5sum";
        else if (std::system("echo foo | md5 -r > /dev/null 2>&1") == 0)
            md5_cmd = "md5 -r";

        ASSERT_FALSE(md5_cmd.empty());
    }

    void run_index_test(std::string const & index_command,
                        std::string const & db_file,
                        std::string const & index_file,
                        std::string const & index_type,
                        std::string const & reduction,
                        std::string const & control_file)
    {
        cli_test_result result = execute_app("lambda3", index_command,
                                             "-d", data(db_file),
                                             "-i", index_file,
                                             "--db-index-type", index_type,
                                             "-r", reduction);

        ASSERT_EQ(result.exit_code, 0);

        if (md5_cmd.empty())
            determine_md5_cmd();

        int ret = std::system((md5_cmd + " " + index_file + " > md5sum_test.txt").c_str());
        ASSERT_TRUE(ret == 0);
        ret = std::system((md5_cmd + " " + (std::string) data(control_file) + " > md5sum_control.txt").c_str());
        ASSERT_TRUE(ret == 0);

        std::ifstream test_output ("md5sum_test.txt");
        std::ifstream control_output ("md5sum_control.txt");

        std::string test_line;
        std::string control_line;

        getline(test_output, test_line);
        getline(control_output, control_line);

        std::string test_sum = test_line.substr(0, test_line.find(' '));
        std::string control_sum = control_line.substr(0, control_line.find(' '));

        ASSERT_EQ(test_sum, control_sum);
    }
};

TEST_F(index_test, no_options)
{
    cli_test_result result = execute_app("lambda3");
    std::string expected
    {
        "lambda3 - Lambda, the Local Aligner for Massive Biological DatA.\n"
        "================================================================\n"
        "    [OPTIONS] COMMAND [COMMAND-OPTIONS]\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(index_test, mkindexn_no_options)
{
    cli_test_result result = execute_app("lambda3", "mkindexn");
    std::string expected
    {
        "lambda3-mkindexn - the Local Aligner for Massive Biological DatA\n"
        "================================================================\n"
        "    [OPTIONS] -d DATABASE.fasta [-i INDEX.lba]\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(index_test, mkindexp_no_options)
{
    cli_test_result result = execute_app("lambda3", "mkindexp");
    std::string expected
    {
        "lambda3-mkindexp - the Local Aligner for Massive Biological DatA\n"
        "================================================================\n"
        "    [OPTIONS] -d DATABASE.fasta [-i INDEX.lba]\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(index_test, nucl_fm)
{
    run_index_test("mkindexn", "db_nucl.fasta.gz", "db_nucl_fm_test.fasta.gz.lba", "fm", "dna4",
                   "db_nucl_fm.fasta.gz.lba");
}

TEST_F(index_test, nucl_bifm)
{
    run_index_test("mkindexn", "db_nucl.fasta.gz", "db_nucl_bifm_test.fasta.gz.lba", "bifm", "dna4",
                   "db_nucl_bifm.fasta.gz.lba");
}

TEST_F(index_test, nucl_bs_fm)
{
    run_index_test("mkindexn", "db_nucl_bs.fasta.gz", "db_nucl_bs_fm_test.fasta.gz.lba", "fm", "dna3bs",
                   "db_nucl_bs_fm.fasta.gz.lba");
}

TEST_F(index_test, nucl_bs_bifm)
{
    run_index_test("mkindexn", "db_nucl_bs.fasta.gz", "db_nucl_bs_bifm_test.fasta.gz.lba", "bifm", "dna3bs",
                   "db_nucl_bs_bifm.fasta.gz.lba");
}

TEST_F(index_test, prot_fm)
{
    run_index_test("mkindexp", "db_prot.fasta.gz", "db_prot_fm_test.fasta.gz.lba", "fm", "li10",
                   "db_prot_fm.fasta.gz.lba");
}

TEST_F(index_test, prot_bifm)
{
    run_index_test("mkindexp", "db_prot.fasta.gz", "db_prot_bifm_test.fasta.gz.lba", "bifm", "li10",
                   "db_prot_bifm.fasta.gz.lba");
}

TEST_F(index_test, trans_fm)
{
    run_index_test("mkindexp", "db_nucl.fasta.gz", "db_trans_fm_test.fasta.gz.lba", "fm", "li10",
                   "db_trans_fm.fasta.gz.lba");
}

TEST_F(index_test, trans_bifm)
{
    run_index_test("mkindexp", "db_nucl.fasta.gz", "db_trans_bifm_test.fasta.gz.lba", "bifm", "li10",
                   "db_trans_bifm.fasta.gz.lba");
}
