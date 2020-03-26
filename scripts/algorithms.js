var arr1 = [5, 3, 2, 6, 7, 9, 1, 0, -1, -7, 10, 8];
var arr2 = [5, 3, 2, 6, 7, 9, 1, 0, -1, -7, 10, 8];
var arr3 = [5, 3, 2, 6, 7, 9, 1, 0, -1, -7, 10, 8];
var arr5 = [5, 3, 2, 6, 7, 9, 1, 0, -1, -7, 10, 8];

window.onload = () => {
    // [2,3,8,6,1] = 5 inversions
    console.log(selectionSort(arr3));
}

/**
 * sorting an array is the process of removing inversions
 * so most sorting algorithms can be easily modified to return 
 * the number of inversions in an array
 */

// O(n^2)
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


// O(n^2)
function selectionSort(arr) {
    let inversions = 0;
    for (let i = 0; i < arr.length; i++) {
        let min = i;
        for (let j = i + 1; j < arr.length; j++)
            if (arr[j] < arr[min]) {
                min = j;
                inversions++;
            }
        arr[i] = [arr[min], arr[min] = arr[i]][0];
    }
    return inversions;
}

// O(n^2)
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

// O(n)
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

// O(nlogn)
function merge_sort(arr, b, e) {
    if (b < e) {
        let q = Math.floor((b + e) / 2)
        merge_sort(arr, b, q);
        merge_sort(arr, q + 1, e);
        merge(arr, b, q, e);
    }
}

// O(nk+nlog(n/k))
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

// O(n)
function linearSearch(val, arr) {
    let j = 0;
    while (j++ < arr.length && arr[j] != val);
    return (j < arr.length) ? j : -1;
}

// O(logn)
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

// O(nlogn + n)
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

// Horner’s rule O(n)
// a_0 + x * ( a_1 + x * (a_2 + ... + x * (a_(n-1) + x*a_n) ...))
function calcPolynomial(x, a) {
    let y = 0;
    for (let i = 0; i < a.length; i++)
        y = a[i] + x * y;
    return y;
}

// O(n^2)
// Inversion : i < j and per[i] > per[j]
// Note :  i+1 < j < per.length
function numberOfInversions(per) {
    let inversions = 0;
    for (let i = 0; i < per.length; i++)
        for (let j = i + 1; j < per.length; j++)
            if (per[i] > per[j]) inversions++;
    return inversions;
}

// O(n logn)
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