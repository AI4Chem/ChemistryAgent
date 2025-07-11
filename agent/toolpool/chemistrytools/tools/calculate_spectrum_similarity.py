from chemistry_tools.spectrum_similarity import SpectrumSimilarity, create_array


def calculate_spectrum_similarity(
    mz_top: list, intensities_top: list, mz_bottom: list, intensities_bottom: list
) -> float:
    """
    Name: calculate_spectrum_similarity
    Description: Calculate spectrum similarity score between two spectra.
    Parameters:
        mz_top : list, List of m/z values for the top spectrum.
        intensities_top : list, List of intensity values corresponding to mz_top.
        mz_bottom : list, List of m/z values for the bottom spectrum.
        intensities_bottom : list, List of intensity values corresponding to mz_bottom.
    Returns:
        score_1: float, Similarity score of spectra_top to spectra_bottom.
        score_2: float, Similarity score of spectra_bottom to spectra_top.
    """
    top_spec = create_array(mz=mz_top, intensities=intensities_top)
    bottom_spec = create_array(mz=mz_bottom, intensities=intensities_bottom)
    spec_sim = SpectrumSimilarity(top_spec, bottom_spec)
    score = spec_sim.score()
    return float(score[0]), float(score[1])


if __name__ == "__main__":
    print(calculate_spectrum_similarity(
        [
            27,
            28,
            32,
            37,
            38,
            39,
            40,
            41,
            50,
            51,
            52,
            53,
            57,
            59,
            61,
            62,
            63,
            64,
            65,
            66,
            67,
            70,
            71,
            72,
            73,
            74,
            75,
            76,
            77,
            78,
            84,
            87,
            88,
            89,
            90,
            91,
            92,
            93,
            101,
            102,
            103,
            104,
            113,
            114,
            115,
            116,
            117,
            126,
            127,
            128,
            129,
            130,
            139,
            140,
            141,
            142,
            143,
            152,
            153,
            154,
            166,
            167,
            168,
            169,
            170,
            171,
        ],
        [
            138,
            210,
            59,
            70,
            273,
            895,
            141,
            82,
            710,
            2151,
            434,
            49,
            41,
            121,
            73,
            229,
            703,
            490,
            1106,
            932,
            68,
            159,
            266,
            297,
            44,
            263,
            233,
            330,
            1636,
            294,
            1732,
            70,
            86,
            311,
            155,
            219,
            160,
            107,
            65,
            111,
            99,
            188,
            107,
            120,
            686,
            150,
            91,
            46,
            137,
            201,
            73,
            69,
            447,
            364,
            584,
            279,
            182,
            37,
            60,
            286,
            718,
            3770,
            6825,
            9999,
            1210,
            85,
        ],
        [
            15,
            26,
            27,
            29,
            30,
            37,
            38,
            39,
            40,
            41,
            42,
            43,
            44,
            50,
            51,
            52,
            53,
            54,
            56,
            62,
            63,
            64,
            65,
            66,
            67,
            68,
            70,
            73,
            74,
            75,
            76,
            77,
            78,
            79,
            80,
            88,
            89,
            90,
            91,
            92,
            93,
            94,
            95,
            102,
            103,
            104,
            105,
            106,
            107,
            115,
            116,
            117,
            118,
            119,
            120,
            121,
            122,
            123,
            127,
            128,
            130,
            131,
            132,
            133,
            134,
            135,
            136,
            139,
            140,
            141,
            146,
            148,
            149,
            150,
            151,
            152,
            153,
            154,
            163,
            164,
            165,
            166,
            167,
            168,
            169,
            175,
            176,
            177,
            179,
            180,
            181,
            182,
            183,
            193,
            194,
            195,
            196,
            197,
            198,
            209,
            210,
            211,
            212,
            223,
            237,
            238,
            239,
            240,
            241,
            251,
            253,
            254,
            267,
            268,
            269,
            270,
        ],
        [
            10,
            10,
            240,
            980,
            30,
            10,
            30,
            180,
            30,
            50,
            90,
            10,
            90,
            100,
            700,
            110,
            30,
            20,
            10,
            20,
            140,
            210,
            400,
            110,
            20,
            10,
            10,
            20,
            20,
            30,
            90,
            3830,
            410,
            160,
            30,
            10,
            20,
            70,
            470,
            1140,
            700,
            90,
            10,
            20,
            200,
            1120,
            650,
            480,
            40,
            10,
            10,
            50,
            340,
            410,
            9999,
            1300,
            100,
            10,
            10,
            10,
            10,
            10,
            30,
            10,
            60,
            20,
            10,
            10,
            10,
            10,
            20,
            7889,
            940,
            70,
            10,
            20,
            10,
            10,
            180,
            1010,
            110,
            30,
            60,
            40,
            10,
            110,
            20,
            10,
            10,
            30,
            30,
            190,
            30,
            10,
            10,
            30,
            130,
            60,
            10,
            10,
            20,
            40,
            10,
            20,
            10,
            10,
            170,
            30,
            10,
            10,
            70,
            10,
            290,
            2940,
            590,
            60,
        ],
    ))
    calculate_spectrum_similarity(
        [
            53,
            55,
            56,
            57,
            66,
            70,
            72,
            79,
            81,
            82,
            85,
            86,
            88,
            96,
            98,
            99,
            100,
            101,
            105,
            109,
            111,
            113,
            114,
            162,
            233,
            250,
            252,
        ],
        [
            24.708788,
            92.449566,
            80.622542,
            100.000000,
            19.276350,
            48.954177,
            21.797781,
            6.331112,
            40.232018,
            36.033970,
            16.470557,
            16.631925,
            6.472862,
            26.350812,
            31.059507,
            8.616107,
            2.344028,
            13.626328,
            72.847200,
            6.071102,
            59.792807,
            23.845567,
            6.604971,
            1.647800,
            2.208517,
            1.070046,
            2.117789,
        ],
        [
            50,
            51,
            52,
            63,
            65,
            74,
            75,
            76,
            77,
            78,
            93,
            107,
            123,
            124,
        ],
        [
            13.113113,
            43.143143,
            2.202202,
            1.301301,
            11.811812,
            4.604605,
            2.902903,
            2.602603,
            100.000000,
            6.606607,
            10.810811,
            1.001001,
            52.652653,
            3.803804,
        ],
    )
    calculate_spectrum_similarity(
        [
            54,
            55,
            56,
            57,
            66,
            70,
            72,
            79,
            81,
            82,
            85,
            86,
            88,
            96,
            98,
            99,
            100,
            101,
            105,
            109,
            111,
            113,
            114,
            162,
            233,
            250,
            252,
        ],
        [
            24.708788,
            92.449566,
            80.622542,
            100.000000,
            19.276350,
            48.954177,
            21.797781,
            6.331112,
            40.232018,
            36.033970,
            16.470557,
            16.631925,
            6.472862,
            26.350812,
            31.059507,
            8.616107,
            2.344028,
            13.626328,
            72.847200,
            6.071102,
            59.792807,
            23.845567,
            6.604971,
            1.647800,
            2.208517,
            1.070046,
            2.117789,
        ],
        [
            53,
            51,
            52,
            63,
            65,
            74,
            75,
            76,
            77,
            78,
            93,
            107,
            123,
            124,
        ],
        [
            13.113113,
            43.143143,
            2.202202,
            1.301301,
            11.811812,
            4.604605,
            2.902903,
            2.602603,
            100.000000,
            6.606607,
            10.810811,
            1.001001,
            52.652653,
            3.803804,
        ],
    )
