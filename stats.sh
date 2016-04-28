cat out.txt | grep "ENTER 0 as 0 for 1" | wc -l
cat out.txt | grep "ENTER 1 as 1 for 0" | wc -l

cat out.txt | grep "ENTER 0 as 0 for 2" | wc -l
cat out.txt | grep "ENTER 2 as 1 for 0" | wc -l

cat out.txt | grep "ENTER 0 as 0 for 3" | wc -l
cat out.txt | grep "ENTER 3 as 1 for 0" | wc -l

cat out.txt | grep "ENTER 1 as 0 for 0" | wc -l
cat out.txt | grep "ENTER 0 as 1 for 1" | wc -l

cat out.txt | grep "ENTER 1 as 0 for 2" | wc -l
cat out.txt | grep "ENTER 2 as 1 for 1" | wc -l

cat out.txt | grep "ENTER 1 as 0 for 3" | wc -l
cat out.txt | grep "ENTER 3 as 1 for 1" | wc -l

cat out.txt | grep "ENTER 2 as 0 for 0" | wc -l
cat out.txt | grep "ENTER 0 as 1 for 2" | wc -l

cat out.txt | grep "ENTER 2 as 0 for 1" | wc -l
cat out.txt | grep "ENTER 1 as 1 for 2" | wc -l

cat out.txt | grep "ENTER 2 as 0 for 3" | wc -l
cat out.txt | grep "ENTER 3 as 1 for 2" | wc -l

cat out.txt | grep "ENTER 3 as 0 for 0" | wc -l
cat out.txt | grep "ENTER 0 as 1 for 3" | wc -l

cat out.txt | grep "ENTER 3 as 0 for 1" | wc -l
cat out.txt | grep "ENTER 1 as 1 for 3" | wc -l

cat out.txt | grep "ENTER 3 as 0 for 2" | wc -l
cat out.txt | grep "ENTER 2 as 1 for 3" | wc -l

#
echo "RA"
cat out.txt | grep "ENTER 1 as Ra 0" | wc -l
cat out.txt | grep "EXIT 1 as Ra 0" | wc -l
cat out.txt | grep "ENTER 1 as Ra 3" | wc -l
cat out.txt | grep "EXIT 1 as Ra 3" | wc -l
cat out.txt | grep "ENTER 2 as Ra 0" | wc -l
cat out.txt | grep "EXIT 2 as Ra 0" | wc -l
cat out.txt | grep "ENTER 2 as Ra 3" | wc -l
cat out.txt | grep "EXIT 2 as Ra 3" | wc -l
cat out.txt | grep "ENTER 3 as Ra 1" | wc -l
cat out.txt | grep "EXIT 3 as Ra 1" | wc -l
cat out.txt | grep "ENTER 3 as Ra 2" | wc -l
cat out.txt | grep "EXIT 3 as Ra 2" | wc -l
cat out.txt | grep "ENTER 0 as Ra 2" | wc -l
cat out.txt | grep "EXIT 0 as Ra 2" | wc -l
cat out.txt | grep "ENTER 0 as Ra 1" | wc -l
cat out.txt | grep "EXIT 0 as Ra 1" | wc -l

# cat out.txt | grep "ENTER 0 as 1" | wc -l
# cat out.txt | grep "ENTER 1 as 1" | wc -l
# cat out.txt | grep "ENTER 2 as 1" | wc -l
# cat out.txt | grep "ENTER 3 as 1" | wc -l
echo "S"
cat out.txt | grep "ENTER 0 as S 1" | wc -l
cat out.txt | grep "EXIT 0 as S 1" | wc -l
cat out.txt | grep "ENTER 0 as S 2" | wc -l
cat out.txt | grep "EXIT 0 as S 2" | wc -l

cat out.txt | grep "ENTER 1 as S 0" | wc -l
cat out.txt | grep "EXIT 1 as S 0" | wc -l
cat out.txt | grep "ENTER 1 as S 3" | wc -l
cat out.txt | grep "EXIT 1 as S 3" | wc -l


cat out.txt | grep "ENTER 2 as S 0" | wc -l
cat out.txt | grep "EXIT 2 as S 0" | wc -l
cat out.txt | grep "ENTER 2 as S 3" | wc -l
cat out.txt | grep "EXIT 2 as S 3" | wc -l


cat out.txt | grep "ENTER 3 as S 1" | wc -l
cat out.txt | grep "EXIT 3 as S 1" | wc -l
cat out.txt | grep "ENTER 3 as S 2" | wc -l
cat out.txt | grep "EXIT 3 as S 2" | wc -l

