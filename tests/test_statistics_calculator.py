from collections import Counter

from astair2.caller import add_to_statistics, statistics_rows


def record(snv, specific="CAT", context="CHH", modified=8, unmodified=2):
    return (
        "some_reference_genome",
        1,
        2,
        0.8,
        modified,
        unmodified,
        "C",
        "T",
        specific,
        context,
        snv,
        10,
    )


def test_modified_position_is_counted():
    statistics = (Counter(), Counter(), Counter())
    add_to_statistics(statistics, record("No"), None)
    covered, modified, unmodified = statistics
    assert covered == {"CHH": 1, "CAT": 1}
    assert modified == {"CHH": 8, "CAT": 8}
    assert unmodified == {"CHH": 2, "CAT": 2}


def test_snv_is_covered_but_not_counted_as_modification():
    statistics = (Counter(), Counter(), Counter())
    add_to_statistics(statistics, record("homozygous"), None)
    covered, modified, unmodified = statistics
    assert covered == {"CHH": 1, "CAT": 1}
    assert not modified and not unmodified


def test_user_context_counts_towards_its_own_row():
    statistics = (Counter(), Counter(), Counter())
    add_to_statistics(
        statistics,
        record(
            "No",
            specific="CAG",
            context="user defined context",
            modified=3,
            unmodified=1,
        ),
        "CAG",
    )
    assert statistics[1] == {"CAG": 3, "user defined context": 3}


def test_statistics_rows():
    statistics = (Counter(), Counter(), Counter())
    add_to_statistics(
        statistics,
        record("No", specific="CGA", context="CpG", modified=3, unmodified=1),
        None,
    )
    rows = list(
        statistics_rows(
            statistics, Counter({"CG": 5, "CGb": 5, "CGA": 4}), "CpG", None, "*"
        )
    )
    assert rows[1] == ("CpG", "*", 75.0, 10, 1, 3, 1)
    assert rows[2] == ("*", "CGA", 75.0, 4, 1, 3, 1)
    assert rows[3] == ("*", "CGC", "*", 0, 0, 0, 0)
