/**
 * \file 
 * \brief Source file for testing Vectors class using googletest.
 * \author Mukund Yedunuthala
 */

// Must include the gtest header to use the testing library
#include "../include/linalg.h"
#include <gtest/gtest.h>

/**
 * \brief Tests if the size of a Vectors object is correctly set after initialization.
 *
 * This test checks that after creating a Vectors object with size 3,
 * the `getSize()` method correctly returns 3.
 */
TEST(TestVectors, TestSizeAfterAssignment) {
    unsigned int size {3};
    Vectors aVec(size);
    ASSERT_EQ(aVec.getSize(), 3);
}

/**
 * \brief Tests if all elements of a Vectors object are correctly assigned a value.
 *
 * This test initializes a Vectors object of size 3, sets all elements to 1.0 
 * using `setToValue(1.0)`, and verifies that the first element is 1.0.
 */
TEST(TestVectors, TestValueAfterAssignment) {
    unsigned int size {3};
    Vectors aVec(size);
    aVec.setToValue(1.0);
    EXPECT_DOUBLE_EQ(aVec.getDataAtIndex(0), 1.0);
}
