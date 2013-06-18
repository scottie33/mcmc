#!/bin/bash

if [ -f "allfig.html" ]; then
	rm -fr allfig.html
	echo " delete allfig.html done!"
fi

format="png"

if [ `ls *.$format | wc -l` -ne 0 ]; then
	rm -f *.$format
	echo " delete *.$format done!"
fi

echo "<html>" > allfig.html
echo "	<body>" >> allfig.html
for i in `ls *bootstrap*.eps r*2*.eps`; do
	if [ -f $i.$format ]; then 
		echo " $i.$format exists, do nothing..."
	else
		convert -rotate 90 $i $i.$format
		echo " $i was converted into $i.$format"
	fi
	echo "	<p>$i.$format</p><br>" >> allfig.html
	echo "	<img src=\"$PWD/$i.$format\"><br>" >> allfig.html
  	#echo "	</a>" >> allfig.html
  	#echo "	</p>" >> allfig.html
done
echo "	</body>" >> allfig.html
echo "</html>" >> allfig.html

echo " allfig.html created, please double-click to check it out."
open allfig.html

exit
