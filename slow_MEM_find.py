def find_mems(pattern, text):
    mems = []
    pattern_length = len(pattern)
    text_length = len(text)

    # Generate all substrings of the pattern
    for start in range(pattern_length):
        for end in range(start + 1, pattern_length + 1):
            substring = pattern[start:end]

            # Find all occurrences of the substring in the text
            pos = text.find(substring)
            if pos != -1:
                # Check for left and right extensions
                left_extendable = (start > 0 and pattern[start-1:end]in text)
                right_extendable = (end + 1 <= pattern_length and pattern[start: end+1] in text)

                if not (left_extendable or right_extendable):
                    mems.append((substring, start))  # Store (substring, start_index)

                # Find next occurrence
    return mems

# Example usage
pattern = "CATGCTATGACTGACCACTTAACTTCTATACCAGTAAACCGACATGCTCTCATATAGGCCCAACAAATTGGCTCGCAGGTCTAATAG"
text = "TTGCATGTTCTCTTGCCGCGATCGTTTGTGGACTCCGCAGATCTGCATGCTATGACTGACGCTTGGAGACGTCCATGTGTACCTCGGTTGACTTAACTTGTATACCAGTAAACCGACATGCTCTCATATAGGCCCAACAAATTGGGTCGCAGGTCTAATAGGCCGGCCACGCTCGTAGTAGCGACGTCATGTAGACGCAAGGCATGATGCCCGCACGGCCAGCGACTTTTTTGCAATCTCAATCGCACCGGGGGAGTGATCAACTCCTTCTCGAGCGCAACGTTTAGTAATGTAGTAGGCACGGGGAGTAAGAAGGGGTCGGCCGGGGGCGGTGAGG$"

# Find MEMs
mems = find_mems(pattern, text)
print("Maximal Exact Matches (MEMs):")
for mem in mems:
    print(f"Match: {mem[0]}, Start Index: {mem[1]}")
