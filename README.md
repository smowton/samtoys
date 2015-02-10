# samtoys
Small htslib-based toys

* **intersect**: Take two ?AM files in qname-sorted order and compute the intersection and/or set difference, regarding records as equal if their qnames and sequences match.
* **subset**: Extract records whose qnames match a specified list, where the input is not assumed qname-sorted.
* **seektest**: Test that seek functionality still appears to work, for developers.

More toys coming as I need them :) Note that some of these tools use my fork of htslib to improve I/O efficiency. You can use that fork to build them, or else just comment out the incompatible changes, such as using hts_set_opt to configure I/O buffer sizes.
