//  arr_ = [0, 1, 2, 3, 4, 5, 6, 7,  8,  9, 10,11];
var arr1 = [5, 3, 2, 6, 7, 9, 1, 0, -1, -7, 10, 8];
var arr2 = [-5, 3, 2, -6, 7, 9, 1, -10, -1, -7, 10, 8];
var arr3 = [5, 3, 2, 6, 7, 9, 1, 0, -1, -7, 10, 8];
var arr5 = [5, 3, 2, 6, 7, 9, 1, 0, -1, -7, 10, 8];
var table;
var data = {
    // A labels array that can contain any sort of values
    labels: [],
    // Our series array that contains series objects or in this case series data arrays
    series: []
};
var matrix1 = [
    [3,-2,5],
    [0,8,4],
    [3,5,1]
];

var matrix2 = [
    [1,2,3],
    [1,2,3],
    [1,2,3]
];
window.onload = () => {
    // let res = testForDifferentInputSizes(
    //     [maxSubArrayRecursive,maxSubArrayLinear]
    //     ,100,10000000,100000);
    // data.labels = res.labels;
    // data.series = res.results;
    // console.log(res.values);
    // new Chartist.Line('.ct-chart', data, {showPoint: false});
    // table = $('#tbody');
    // let lcs = getAllSebsquences('ABCDEF'.split(''),4);
    // printArray(lcs);
    // console.log(longestCommonSubsequenceBruteforce('ACDCABDBADCACDCDBCCD'.split(''),
    //                                             'CBDDDDDBCCAAAABBCDAC'.split('')));
}

function testForDifferentInputSizes (algorithms, minN=10, maxN=20000, step=100,
    rangeMin=-1000, rangeMax=1000, floats=false) {
    let arr = [];
    let labels = []; // holds sample size labels
    let results = []; // holds test results
    let values = []; // holds test return values
    for (let size=minN; size<=maxN; size += step) {
        // console.log('Started testing ' + algo.name + ' with sample size ' + size);
        
        // fill the array to the next sample size
        for (let i=arr.length; i<size; i++) {
            const rand = Math.random() * (rangeMax - rangeMin);
            let val = rangeMin + (floats ? rand : Math.floor(rand));
            arr.push(val);
        }
        // labels.push(size);
        // test the given algorithms on this sample size
        for (let i=0; i<algorithms.length; i++) {
            let start = performance.now();
            let val = (algorithms[i])(arr.slice());
            let elapsed = performance.now() - start;
            // console.log('Finished testing ' + algo.name + ' with sample size ' + size);
            if (!results[i]) results[i] = [];
            results[i].push(elapsed);
            if (!values[i]) values[i] = [];
            values[i].push(val);
        }
    }
    
    return {labels: labels, results: results, values: values};
}

function generateRandomArray (size, array, offset=0, rangeMin=-1000, rangeMax=1000, floats=false) {
    let arr = array || [];
    offset = array ? array.length : 0;
    // fill the array to the next sample size
    for (let i=offset; i<size; i++) {
        const rand = Math.random() * (rangeMax - rangeMin);
        let val = rangeMin + (floats ? rand : Math.floor(rand));
        arr.push(val);
    }
    return arr;
}

function generateRandomMatrix(n, rangeMin=-1000, rangeMax=1000, floats=false) {
    let res = [];
    // let log = '';
    for (let i=0; i<n; i++) {
        res[i] = [];
        for (let j=0; j<n; j++) {
            const rand = Math.random() * (rangeMax - rangeMin);
            res[i][j] = rangeMin + (floats ? rand : Math.floor(rand));
            // log += res[i][j] + ' ';
        }
        // log += '\n';
    }
    // console.log(log);
    return res;
}

function printArray(arr, _log=true, printSize=true) {
    let log = '[';
    for (let i=0; i<arr.length; i++) {
        if (arr[i] instanceof Array)
            log += printArray(arr[i], false, false) + ',';
        else
            log += (i === 0) ? arr[i] : ',' + arr[i];
    }
    log += (printSize) ?  '] => ' + arr.length + ' elements' : ']';
    if (_log)
        console.log(log);
    else
        return log;
}

function print2DArray(arr) {
    let log = '[\n  [';
    for (let i=0; i<arr.length; i++) {
        log += arr[i].join(', ');
        log += (i === arr.length-1) ? ']' : ']\n  [';
    }
    log += '\n]';
    console.log(log);
}

/**
 * sorting an array is the process of removing inversions
 * so most sorting algorithms can be easily modified to return 
 * the number of inversions in an array
 */

// ϴ(n^2)
function bubbleSort(arr) {
    let inversions = 0;
    for (let i = 0; i < arr.length - 1; i++) {
        let noSwaps = true;
        for (let j = arr.length; j > i; j--)
            if (arr[j] < arr[j - 1]) {
                arr[j - 1] = [arr[j], arr[j] = arr[j - 1]][0];
                noSwaps = false;
                inversions++;
            }

        if (noSwaps)
            break;
    }
    return inversions;
}


// ϴ(n^2)
function selectionSort(arr) {
    for (let i = 0; i < arr.length; i++) {
        let min = i;
        for (let j = i + 1; j < arr.length; j++) {
            if (arr[j] < arr[min]) {
                min = j;
            }
        }
        arr[i] = [arr[min], arr[min] = arr[i]][0];
    }
}

// ϴ(n^2)
function insertionSort(arr, decreasing = false) {
    let inversions = 0;
    for (let i = 1; i < arr.length; i++) {
        const ele = arr[i];
        let j = i - 1;
        while (j >= 0 && ((arr[j] > ele) !== decreasing)) {
            arr[j + 1] = arr[j];
            --j;
        }
        arr[j + 1] = ele;
        inversions += i-j-1;
    }
    return inversions;
}

function insertionSortSubArr(arr, b, e) {
    for (let i = b + 1; i < e + 1; i++) {
        const ele = arr[i];
        let j = i - 1;
        while (j >= b && arr[j] > ele) {
            arr[j + 1] = arr[j];
            --j;
        }
        arr[j + 1] = ele;
    }
}

// convenience function
function recurInsertionSort(arr) {
    recur_insertion_sort(arr, 0, arr.length - 1);
}

function recur_insertion_sort(arr, b, e) {
    if (b < e) {
        //recursively sort arr[b.. e-1]
        recur_insertion_sort(arr, b, e - 1);
        //insert arr[e] into the already sorted arr[b.. e-1]
        const ele = arr[e];
        let j = e - 1;
        while (j >= 0 && arr[j] > ele) {
            arr[j + 1] = arr[j];
            --j;
        }
        arr[j + 1] = ele;
    }
}

// convenience function
function mergeSort(arr) {
    merge_sort(arr, 0, arr.length - 1);
}

// ϴ(n)
function merge(arr, a, b, c) {
    let arr1 = arr.slice(a, b + 1);
    arr1.push(Infinity); // sentinal card
    let arr2 = arr.slice(b + 1, c + 1);
    arr2.push(Infinity); // sentinal card

    let j = 0,
        k = 0,
        inversions = 0; // not required

    for (let i = a; i <= c; i++) {
        if (arr1[j] < arr2[k])
            arr[i] = arr1[j++];
        else {
            arr[i] = arr2[k++];
            // count it as an inversion only if left and right arrays are not empty
            // (if sentinal cards are not exposed)
            if (j<b && k<c)
                inversions += b - j;
        }
    }
    return inversions;
}

// ϴ(nlogn)
function merge_sort(arr, b, e) {
    if (b < e) {
        let q = Math.floor((b + e) / 2)
        merge_sort(arr, b, q);
        merge_sort(arr, q + 1, e);
        merge(arr, b, q, e);
    }
}

// ϴ(nk+nlog(n/k))
function hybridMergeInsertionSort(arr, k = 20) {
    hybrid_merge_ins_sort(arr, 0, arr.length - 1, k);
}

function hybrid_merge_ins_sort(arr, b, e, k) {
    // add 1 because we start from 0
    if ((e - b + 1) > k) {
        let q = Math.floor((b + e) / 2)
        hybrid_merge_ins_sort(arr, b, q, k);
        hybrid_merge_ins_sort(arr, q + 1, e, k);
        merge(arr, b, q, e);
    } else {
        insertionSortSubArr(arr, b, e);
    }
}

// ϴ(n)
function linearSearch(val, arr) {
    let j = 0;
    while (j++ < arr.length && arr[j] != val);
    return (j < arr.length) ? j : -1;
}

function recursiveLinearSearch(val, arr) {
    // base cases
    if (!arr.length) return -1;
    if (arr[arr.length-1] === val) return arr.length-1;
    return recursiveLinearSearch(val, arr.slice(0,arr.length - 1))
}

// ϴ(logn)
function binarySearch(val, arr) {
    let l = 0,
        r = arr.length - 1;
    while (l <= r) {
        let mid = Math.floor((l + r) / 2);
        if (arr[mid] == val)
            return mid;
        else if (arr[mid] < val)
            l = mid + 1;
        else
            r = mid - 1;
    }
    return -1;
}

// ϴ(nlogn + n)
function areThere2ValWithSum(arr, val) {
    mergeSort(arr); //nlogn
    for (let i = 0; i < arr.length; i++) //n
        if (binarySearch(val - arr[i], arr) != -1) //logn
            return true;
    return false;
}


// adds two binary numbers in form of arrays
function add2BinaryNumbers(a, b) {
    const longest = (a.length > b.length) ? a.length : b.length;
    let result = [];
    a.reverse();
    b.reverse();
    let carriage = 0;
    for (let i = 0; i < longest; i++) {
        const A = (i < a.length) ? a[i] : 0;
        const B = (i < b.length) ? b[i] : 0;
        result.push(A ^ B ^ carriage);
        if (A + B + carriage > 1) carriage = 1;
        else carriage = 0;
    }
    if (carriage == 1) result.push(carriage);
    result.reverse();
    return result;
}

// Horner’s rule ϴ(n)
// a_0 + x * ( a_1 + x * (a_2 + ... + x * (a_(n-1) + x*a_n) ...))
function calcPolynomial(x, a) {
    let y = 0;
    for (let i = 0; i < a.length; i++)
        y = a[i] + x * y;
    return y;
}

// ϴ(n^2)
// Inversion : i < j and per[i] > per[j]
// Note :  i+1 < j < per.length
function numberOfInversions(per) {
    let inversions = 0;
    for (let i = 0; i < per.length; i++)
        for (let j = i + 1; j < per.length; j++)
            if (per[i] > per[j]) inversions++;
    return inversions;
}

// ϴ(n logn)
// number of inversions with modified merge sort
// convenience function
function numberOfInversionsModifiedMergeSort(arr) {
    return _numberOfInversionsModifiedMergeSort(arr, 0, arr.length - 1);
}

function _numberOfInversionsModifiedMergeSort(arr, b, e) {
    if (b < e) {
        let q = Math.floor((b + e) / 2)
        let i = _numberOfInversionsModifiedMergeSort(arr, b, q);
        let j = _numberOfInversionsModifiedMergeSort(arr, q + 1, e);
        return i + j + merge(arr, b, q, e);
    }
    return 0;
}

// find the subarray with the highest sum
// ϴ(n^2)
function maxSubArrayBrutforce(arr) {
    let maxSubArr = {
        start: 0,
        end: 0,
        sum: -Infinity
    };
    for (let i=0; i<arr.length; i++) {
        let sum = 0;
        for (let j=i; j< arr.length; j++) {
            sum += arr[j];
            if(sum >= maxSubArr.sum)
                maxSubArr = {start: i, end: j, sum: sum};
        }
    }
    return maxSubArr;
}

// ϴ(n logn)
// convenience funtion
function maxSubArrayRecursive(arr) {
    return _maxSubArrayRecursive(arr,0,arr.length - 1);
}

function _maxSubArrayRecursive(arr, start, end) {
    // base case (1 element in array)
    if (start === end) return {start: start, end: end, sum: arr[start]}
    /**
     * else solve recursively
     * the maximum subarray may fall in totally in the left half,
     * totally in the right half, or in between
     *  */
    let mid = Math.floor((start+end)/2);
    let right = _maxSubArrayRecursive(arr,mid+1,end);
    let left = _maxSubArrayRecursive(arr,start,mid);
    let crossing = maxCrossingSubArray(arr,start,mid,end);
    // for debugging 
    // let row = '<tr><td>' + start + '</td><td>' + mid + '</td><td>' + end + '</td><td>'
    //             + '[' + arr.slice(start,end+1).join(' ') + ']' + '</td><td>'
    //             + '(' + left.start + ',' + left.end + ',' + left.sum + ')</td><td>' 
    //             + '(' + crossing.start + ',' + crossing.end + ',' + crossing.sum + ')</td><td>'
    //             + '(' + right.start + ',' + right.end + ',' + right.sum + ')</td></tr>';
    // table.prepend(row);
    // this can be replacesd by an array sort
    if (left.sum >= right.sum && left.sum >= crossing.sum)
        return left;
    else if (right.sum >= left.sum && right.sum >= crossing.sum)
        return right;
    else return crossing;
}

// for array of 6 elements (arr.length) mid should be 2
// for array of 5 elements (arr.length) mid should be 1
function maxCrossingSubArray(arr, start, mid, end) {
    // find maximum left subarray
    let maxSubArrLeft = {
        start: 0,
        end: mid,
        sum: -Infinity
    };
    let tempSum = 0;
    // for array of 6 this loops through INDEXES 0,1,2
    for (let i=mid; i>=start; i--) {
        tempSum += arr[i];
        if (tempSum > maxSubArrLeft.sum)
            maxSubArrLeft = {start: i, end: mid, sum: tempSum};
    }
    // find maximum right subarray
    let maxSubArrRight = {
        start: 0,
        end: mid,
        sum: -Infinity
    };
    tempSum = 0;
    // for array of 6 this loops through INDEXES 3,4,5
    for (let i=mid+1; i<=end; i++) {
        tempSum += arr[i];
        if (tempSum > maxSubArrRight.sum)
            maxSubArrRight = {start: mid, end: i, sum: tempSum};
    }

    return {start: maxSubArrLeft.start,
        end: maxSubArrRight.end,
        sum: maxSubArrLeft.sum + maxSubArrRight.sum}
}

/**
 * Kadane's Algorithm
 * explanation : at each iteration `i` find the longest subarray that ends
 * at the current index (the current index MUST be included in the subarray)
 * so at each iteration one has two choices, either extend the last maximum
 * subarray, ending at index `i-1`, or start a new maximum subarray starting
 * at index `i`, the compare the sum of this maiximum "local" subarray to the
 * global maximum subarray
 * ϴ(n)
 */
function maxSubArrayLinear(arr) {
    let msa = {start: 0, end:0, sum: arr[0]};
    let general = {start: 0, end:0, sum: arr[0]};
    for (let i=1; i<arr.length; i++) {
        // extend current maximum subarray
        // or start a new one
        if (msa.sum + arr[i] > arr[i]) {
            msa.end = i;
            msa.sum += arr[i];
        } else {
            msa.start = msa.end = i;
            msa.sum = arr[i];
        }

        if (msa.sum > general.sum) {
            general = {...msa} // shallow copy;
        }
    }

    return general;
}

// ϴ(n^3)
function multiplySquareMatrixNormalMehtod(matrix1, matrix2) {
    let res = [];
    let n = matrix1.length;
    for (let i = 0; i < n; i++) {
        res[i] = [];
        for (let j = 0; j < n; j++) {
            res[i][j] = 0;
            for (let k=0; k<n; k++)
                res[i][j] += matrix1[i][k] * matrix2[k][j];
        }
    }
    return res;
}

function getAllSebsquences(arr, size, offset = 0) {
    // base cases
    if (size === 1)
        return arr.slice(offset).map(val => [val]);
    else if (size === 0)
        return [arr];
    /**
     * for each position list all the possible subsequences
     * by calling the function recursively then prepending the
     * current element, like listing all possible combinations
     * with two nested for loops
     */
    let subs = [];
    for (let i = offset; i < arr.length; i++) {
        let _subs = getAllSebsquences(arr, size - 1, i + 1);
        _subs.forEach(subArray => subs.push([arr[i], ...subArray]));
    }
    return subs;
}

// ϴ(n*2^n)
// this will break easily for large inputs
function longestCommonSubsequenceBruteforce(arr1, arr2) {
    let lcs = '';
    // let arr be the shorter and secondArr the longer array
    // because the sequence length can be a maximum of arr length
    let arr, secondArr;
    if (arr1.length <= arr2.length) {
        arr = arr1; secondArr = arr2;
    } else {
        arr = arr2; secondArr = arr1;
    }

    /**
     * list all possible subsequences of arr of all sizes
     * and for each one check if it exists in secondArr
     * if it does update the longest common subsequence
     */
    for (let i = 1; i <= arr.length; i++) {
        let subs = getAllSebsquences(arr, i);
        console.log(`Number of Subsequences of size ${i} : ${subs.length}`);
        
        for (let j = 0; j < subs.length; j++) {
            let subsFound = true;
            let lastIndex = -1;
            for (let k = 0; k < subs[j].length; k++) {
                let index = secondArr.indexOf(subs[j][k], lastIndex + 1);
                if (index === -1) {
                    subsFound = false; break;
                } else
                    lastIndex = index;
            }
            if (subsFound && subs[j].length > lcs.length)
                lcs = subs[j].join('');
        }
    }
    return lcs;
}

// ϴ(n^2)
// Dynamic Programming 
function longestCommonSubsequence(arr1, arr2) {
    let dpTable = [];
    let pathTable = [];
    dpTable[0] = new Array(arr2.length+1).fill(0);
    // build the dp table and the path table
    for (let i=0; i<arr1.length; i++) {
        dpTable[i+1] = [0];
        pathTable[i] = [];
        for (let j=0; j<arr2.length; j++) {
            if ( arr1[i] === arr2[j]) {
                dpTable[i+1][j+1] = 1 + dpTable[i][j];
                pathTable[i][j] = '↖';
            } else
                if (dpTable[i][j+1] >= dpTable[i+1][j]) {
                    dpTable[i+1][j+1] = dpTable[i][j+1];
                    pathTable[i][j] = '↑';
                } else {
                    dpTable[i+1][j+1] = dpTable[i+1][j];
                    pathTable[i][j] = '←';
                }
        }
    }
    // follow the path to find the subsequence
    let i=arr1.length-1, j=arr2.length-1;
    let lcs = [];
    while (i>=0 && j>=0) {
        if (pathTable[i][j] === '↖') {
            lcs.push(arr1[i]);
            i--; j--; // go diagnol
        } else if (pathTable[i][j] === '↑') {
            i--; // go up
        } else {
            j--; // go left
        }
    }
    return lcs.reverse().join('');
    // return dpTable[arr1.length][arr2.length];
}