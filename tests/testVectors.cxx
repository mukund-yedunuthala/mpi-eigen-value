// Must include the gtest header to use the testing library
#include "../include/linalg.h"
#include <gtest/gtest.h>

TEST(TestVectors, TestSizeAfterAssignment) {
    unsigned int size {3};
    Vectors aVec(size);
    ASSERT_EQ(aVec.getSize(), 3);
}

TEST(TestVectors, TestValueAfterAssignment) {
    unsigned int size {3};
    Vectors aVec(size);
    aVec.setToValue(1.0);
    EXPECT_DOUBLE_EQ(aVec.getDataAtIndex(0), 1.0);
}