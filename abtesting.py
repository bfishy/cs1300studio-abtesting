from scipy import stats
from scipy.stats import t as t_dist
from scipy.stats import chi2

from abtesting_test import *

# # You can comment out these lines! They are just here to help follow along to the tutorial.
# print(t_dist.cdf(-2, 20)) # should print .02963
# print(t_dist.cdf(2, 20)) # positive t-score (bad), should print .97036 (= 1 - .2963)

# print(chi2.cdf(23.6, 12)) # prints 0.976
# print(1 - chi2.cdf(23.6, 12)) # prints 1 - 0.976 = 0.023 (yay!)

# TODO: Fill in the following functions! Be sure to delete "pass" when you want to use/run a function!
# NOTE: You should not be using any outside libraries or functions other than the simple operators (+, **, etc)
# and the specifically mentioned functions (i.e. round, cdf functions...)

def slice_2D(list_2D, start_row, end_row, start_col, end_col):
    '''
    Splices a the 2D list via start_row:end_row and start_col:end_col
    :param list: list of list of numbers
    :param nums: start_row, end_row, start_col, end_col
    :return: the spliced 2D list (ending indices are exclsive)
    '''
    to_append = []
    for l in range(start_row, end_row):
        to_append.append(list_2D[l][start_col:end_col])

    return to_append

def get_avg(nums):
    '''
    Helper function for calculating the average of a sample.
    :param nums: list of numbers
    :return: average of list
    '''
    s = 0
    n = len(nums)
    for num in nums:
        s += num
    return s / n

def get_stdev(nums):
    '''
    Helper function for calculating the standard deviation of a sample.
    :param nums: list of numbers
    :return: standard deviation of list
    '''
    n = len(nums)
    avg = get_avg(nums)
    s = 0
    for num in nums:
        diff = (num - avg) ** 2
        s += diff
    std = (s/(n-1))**0.5
    return std

def get_standard_error(a, b):
    '''
    Helper function for calculating the standard error, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: standard error of a and b (see studio 6 guide for this equation!)
    '''
    na = len(a)
    nb = len(b)
    std_a = get_stdev(a)
    std_b = get_stdev(b)
    return ((std_a**2/na) + (std_b**2/nb))**0.5

def get_2_sample_df(a, b):
    '''
    Calculates the combined degrees of freedom between two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: integer representing the degrees of freedom between a and b (see studio 6 guide for this equation!)
    HINT: you can use Math.round() to help you round!
    '''
    se = get_standard_error(a, b)
    std_a = get_stdev(a)
    std_b = get_stdev(b)
    n_a = len(a)
    n_b = len(b)
    df = round(se**4 / (
        (((std_a**2)/n_a)**2 / (n_a-1))
        + (((std_b**2)/n_b)**2 / (n_b-1))))
    return df

def get_t_score(a, b):
    '''
    Calculates the t-score, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: number representing the t-score given lists a and b (see studio 6 guide for this equation!)
    '''
    se = get_standard_error(a, b)
    t_score = ((get_avg(a) - get_avg(b)) / se)
    t_score = t_score if t_score < 0 else t_score * (-1)
    return t_score

def perform_2_sample_t_test(a, b):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates a p-value by performing a 2-sample t-test, given two lists of numbers.
    :param a: list of numbers
    :param b: list of numbers
    :return: calculated p-value
    HINT: the t_dist.cdf() function might come in handy!
    '''
    p_value = t_dist.cdf(get_t_score(a, b), get_2_sample_df(a, b))
    return p_value


# [OPTIONAL] Some helper functions that might be helpful in get_expected_grid().
def row_sum(observed_grid, ele_row):
    s = 0
    m = len(observed_grid[ele_row])
    for i in range(m):
        s += observed_grid[ele_row][i]
    return s
def col_sum(observed_grid, ele_col):
    s = 0
    n = len(observed_grid)
    for i in range(n):
        s += observed_grid[i][ele_col]
    return s
def total_sum(observed_grid):
    s = 0
    n = len(observed_grid)
    for i in range(n):
        s += row_sum(observed_grid, i)
    return s
def calculate_expected(row_sum, col_sum, tot_sum):
    expected = (row_sum * col_sum) / tot_sum
    return expected

def get_expected_grid(observed_grid):
    '''
    Calculates the expected counts, given the observed counts.
    ** DO NOT modify the parameter, observed_grid. **
    :param observed_grid: 2D list of observed counts
    :return: 2D list of expected counts
    HINT: To clean up this calculation, consider filling in the optional helper functions below!
    '''
    n = len(observed_grid)
    m = len(observed_grid[0])
    result = [[0 for i in range(m)] for j in range(n)]
    for i in range(n):
        row_s = row_sum(observed_grid, i)
        for j in range(m):
            col_s = col_sum(observed_grid, j)
            total_s = total_sum(observed_grid)
            result[i][j] = calculate_expected(row_s, col_s, total_s)
    return result
# print(get_expected_grid([[207, 282], [231, 242]]))
def df_chi2(observed_grid):
    '''
    Calculates the degrees of freedom of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: degrees of freedom of expected counts (see studio 6 guide for this equation!)
    '''
    n = len(observed_grid)
    m = len(observed_grid[0])
    return (n-1) * (m-1)

def chi2_value(observed_grid):
    '''
    Calculates the chi^2 value of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: associated chi^2 value of expected counts (see studio 6 guide for this equation!)
    '''
    s = 0
    n = len(observed_grid)
    m = len(observed_grid[0])
    expected = get_expected_grid(observed_grid)
    for i in range(n):
        for j in range(m):
            s += (observed_grid[i][j] - expected[i][j]) ** 2 / expected[i][j]
    return s  

def perform_chi2_homogeneity_test(observed_grid):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates the p-value by performing a chi^2 test, given a list of observed counts
    :param observed_grid: 2D list of observed counts
    :return: calculated p-value
    HINT: the chi2.cdf() function might come in handy!
    '''
    return 1 - chi2.cdf(chi2_value(observed_grid), df_chi2(observed_grid))

# These commented out lines are for testing your main functions. 
# Please uncomment them when finished with your implementation and confirm you get the same values :)
def data_to_num_list(s):
  '''
    Takes a copy and pasted row/col from a spreadsheet and produces a usable list of nums. 
    This will be useful when you need to run your tests on your cleaned log data!
    :param str: string holding data
    :return: the spliced list of numbers
    '''
  return list(map(float, s.split()))


# # t_test 1:
# a_t1_list = data_to_num_list(a1) 
# b_t1_list = data_to_num_list(b1)
# print(get_t_score(a_t1_list, b_t1_list)) # this should be -129.500
# print(perform_2_sample_t_test(a_t1_list, b_t1_list)) # this should be 0.0000
# # why do you think this is? Take a peek at a1 and b1 in abtesting_test.py :)

# # t_test 2:
# a_t2_list = data_to_num_list(a2) 
# b_t2_list = data_to_num_list(b2)
# print(get_t_score(a_t2_list, b_t2_list)) # this should be -1.48834
# print(perform_2_sample_t_test(a_t2_list, b_t2_list)) # this should be .082379

# # t_test 3:
# a_t3_list = data_to_num_list(a3) 
# b_t3_list = data_to_num_list(b3)
# print(get_t_score(a_t3_list, b_t3_list)) # this should be -2.88969
# print(perform_2_sample_t_test(a_t3_list, b_t3_list)) # this should be .005091

print('#######################################')

# chi2_test 1:
a_c1_list = data_to_num_list(a_count_1) 
b_c1_list = data_to_num_list(b_count_1)
c1_observed_grid = [a_c1_list, b_c1_list]
print(chi2_value(c1_observed_grid)) # this should be 4.103536
print(perform_chi2_homogeneity_test(c1_observed_grid)) # this should be .0427939

# chi2_test 2:
a_c2_list = data_to_num_list(a_count_2) 
b_c2_list = data_to_num_list(b_count_2)
c2_observed_grid = [a_c2_list, b_c2_list]
print(chi2_value(c2_observed_grid)) # this should be 33.86444
print(perform_chi2_homogeneity_test(c2_observed_grid)) # this should be 0.0000
# Again, why do you think this is? Take a peek at a_count_2 and b_count_2 in abtesting_test.py :)

# chi2_test 3:
a_c3_list = data_to_num_list(a_count_3) 
b_c3_list = data_to_num_list(b_count_3)
c3_observed_grid = [a_c3_list, b_c3_list]
print(chi2_value(c3_observed_grid)) # this should be .3119402
print(perform_chi2_homogeneity_test(c3_observed_grid)) # this should be .57649202


