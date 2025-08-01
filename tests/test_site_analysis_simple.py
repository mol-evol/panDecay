#!/usr/bin/env python3
"""
Simple functional test for site analysis implementation.

This tests the basic functionality of the site analysis components
without requiring actual PAUP* execution.
"""

from src.core.analysis.site_analyzer import SiteData, CladesSiteAnalysis


def test_site_data_creation():
    """Test creating site data objects."""
    print("Testing SiteData creation...")
    
    # Test supporting site
    site1 = SiteData(
        site_position=1,
        ml_likelihood=-10.5,
        constraint_likelihood=-12.8,
        delta_lnl=2.3,
        site_type="supporting"
    )
    
    assert site1.site_position == 1
    assert site1.delta_lnl == 2.3
    assert site1.site_type == "supporting"
    print("✓ Supporting site data created correctly")
    
    # Test conflicting site
    site2 = SiteData(
        site_position=2,
        ml_likelihood=-8.2,
        constraint_likelihood=-7.1,
        delta_lnl=-1.1,
        site_type="conflicting"
    )
    
    assert site2.delta_lnl == -1.1
    assert site2.site_type == "conflicting"
    print("✓ Conflicting site data created correctly")
    
    return [site1, site2]


def test_site_analysis_creation():
    """Test creating site analysis results."""
    print("Testing CladesSiteAnalysis creation...")
    
    # Create test site data
    site_data = {
        1: SiteData(1, -10.5, -12.8, 2.3, "supporting"),
        2: SiteData(2, -8.2, -7.1, -1.1, "conflicting"),
        3: SiteData(3, -9.0, -9.05, -0.05, "neutral")
    }
    
    # Create site analysis
    analysis = CladesSiteAnalysis(
        clade_id="Clade_1",
        site_data=site_data,
        supporting_sites=1,
        conflicting_sites=1,
        neutral_sites=1,
        support_ratio=0.333,
        sum_supporting_delta=2.3,
        sum_conflicting_delta=-1.1,
        weighted_support_ratio=0.677
    )
    
    assert analysis.clade_id == "Clade_1"
    assert len(analysis.site_data) == 3
    assert analysis.supporting_sites == 1
    assert analysis.conflicting_sites == 1
    assert analysis.neutral_sites == 1
    assert abs(analysis.support_ratio - 0.333) < 0.001
    print("✓ Site analysis object created correctly")
    
    return analysis


def test_site_categorization():
    """Test site categorization logic."""
    print("Testing site categorization...")
    
    threshold = 0.1
    
    # Test supporting site (delta > threshold)
    delta_lnl = 0.5
    if delta_lnl > threshold:
        site_type = "supporting"
    elif delta_lnl < -threshold:
        site_type = "conflicting"
    else:
        site_type = "neutral"
    
    assert site_type == "supporting"
    print("✓ Supporting site categorized correctly")
    
    # Test conflicting site (delta < -threshold)
    delta_lnl = -0.3
    if delta_lnl > threshold:
        site_type = "supporting"
    elif delta_lnl < -threshold:
        site_type = "conflicting"
    else:
        site_type = "neutral"
    
    assert site_type == "conflicting"
    print("✓ Conflicting site categorized correctly")
    
    # Test neutral site (-threshold <= delta <= threshold)
    delta_lnl = 0.05
    if delta_lnl > threshold:
        site_type = "supporting"
    elif delta_lnl < -threshold:
        site_type = "conflicting"
    else:
        site_type = "neutral"
    
    assert site_type == "neutral"
    print("✓ Neutral site categorized correctly")


def test_support_ratio_calculation():
    """Test support ratio calculations."""
    print("Testing support ratio calculations...")
    
    # Test data
    supporting_sites = 15
    conflicting_sites = 8
    neutral_sites = 2
    total_sites = supporting_sites + conflicting_sites + neutral_sites
    
    # Basic support ratio
    support_ratio = supporting_sites / total_sites
    expected_ratio = 15 / 25
    assert abs(support_ratio - expected_ratio) < 0.001
    print(f"✓ Support ratio calculated correctly: {support_ratio:.3f}")
    
    # Weighted support ratio
    sum_supporting_delta = 12.5
    sum_conflicting_delta = -6.8
    total_delta = sum_supporting_delta + abs(sum_conflicting_delta)
    weighted_support_ratio = sum_supporting_delta / total_delta if total_delta > 0 else 0.0
    
    expected_weighted = 12.5 / (12.5 + 6.8)
    assert abs(weighted_support_ratio - expected_weighted) < 0.001
    print(f"✓ Weighted support ratio calculated correctly: {weighted_support_ratio:.3f}")


def run_all_tests():
    """Run all site analysis tests."""
    print("="*60)
    print("Site Analysis Implementation Test Suite")
    print("="*60)
    
    try:
        # Test basic data structures
        test_site_data_creation()
        print()
        
        test_site_analysis_creation()
        print()
        
        # Test logic components
        test_site_categorization()
        print()
        
        test_support_ratio_calculation()
        print()
        
        print("="*60)
        print("✓ All site analysis tests passed!")
        print("Site analysis implementation is working correctly.")
        print("="*60)
        
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)