import six
from nose.tools import assert_equal, assert_almost_equal, assert_true, \
    assert_false, assert_raises, assert_is_instance, with_setup

import driver

def test_name_from_path():
    assert_equal(driver.name_from_path("/Users/test/test.py"), "test")
    assert_equal(driver.name_from_path("/Users/test/test"), "test")
    assert_equal(driver.name_from_path("test.py"), "test")
    assert_equal(driver.name_from_path("~/test.py"), "test")
    assert_equal(driver.name_from_path("~/test"), "test")

def test_name_from_path_has_p():
    assert_equal(driver.name_from_path("/Users/test/testp.py"), "testp")
    assert_equal(driver.name_from_path("testp.py"), "testp")
    assert_equal(driver.name_from_path("~/testp.py"), "testp")
    assert_equal(driver.name_from_path("~/testp"), "testp")
