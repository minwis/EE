# EE
RQ: To what extent does SIMD vectorization using a striped layout improve the computational efficiency of the Smith-Waterman?

Refer to these journals:
1. https://academic.oup.com/bioinformatics/article/23/2/156/205631
2. https://www.researchgate.net/publication/354390311_A_Review_of_Parallel_Implementations_for_the_Smith-Waterman_Algorithm
3. https://arxiv.org/pdf/1909.00899

How to run program(IntelliJ)
When creating application to run the program, 
1. Browse MAIN class 
2. Click "Modify options"
3. Click "Add VM options"
4. Copy and pase: "--add-modules jdk.incubator.vector"
5. Save and run