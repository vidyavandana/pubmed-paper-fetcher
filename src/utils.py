import re
from typing import Union, Optional, List
NON_ACADEMIC_KEYWORDS = [
    "inc", "ltd", "llc", "biotech", "pharma", "technologies", "solutions", "corporation", "company", "gmbh", "s.a.", "co.", "pvt"
]
ACADEMIC_KEYWORDS = [
    "university", "institute", "hospital", "college", "school", "center", "centre", "department", "faculty"
]

def is_non_academic_author(affiliation: str) -> bool:
    """
    Returns True if the affiliation looks like a company (non-academic), based on keywords.
    """
    affil = affiliation.lower()
    return (
        any(word in affil for word in NON_ACADEMIC_KEYWORDS)
        and not any(word in affil for word in ACADEMIC_KEYWORDS)
    )

def extract_email(affiliations: Union[str, List[str]]) -> Optional[str]:
    """
    Tries to extract an email address from the affiliation string(s).
    """
    email_pattern = r"[\w\.-]+@[\w\.-]+\.\w+"

    if isinstance(affiliations, str):
        matches = re.findall(email_pattern, affiliations)
        return matches[0] if matches else None

    if isinstance(affiliations, list):
        for aff in affiliations:
            matches = re.findall(email_pattern, aff)
            if matches:
                return matches[0]

    return None
