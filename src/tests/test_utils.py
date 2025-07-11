# tests/test_utils.py

from src.utils import is_non_academic_author, extract_email

def test_is_non_academic_author_positive():
    assert is_non_academic_author("Incyte Corporation, Wilmington, USA") is True
    assert is_non_academic_author("Roche Pharma (Schweiz) AG") is True
    assert is_non_academic_author("Biotech Solutions Pvt Ltd") is True

def test_is_non_academic_author_negative():
    assert is_non_academic_author("Department of Oncology, Stanford University") is False
    assert is_non_academic_author("Massachusetts General Hospital") is False
    assert is_non_academic_author("Harvard Medical School") is False

def test_extract_email_single_string():
    aff = "Roche Pharma (Schweiz) AG, Basel, Switzerland. Email: john.doe@roche.com"
    assert extract_email(aff) == "john.doe@roche.com"

def test_extract_email_list_of_affiliations():
    affs = [
        "Stanford University, CA, USA.",
        "Genentech Inc., San Francisco, USA. Contact: jane.smith@genentech.com"
    ]
    assert extract_email(affs) == "jane.smith@genentech.com"

def test_extract_email_not_found():
    assert extract_email("No email here at all") is None
