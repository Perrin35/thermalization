OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(-0.17012574) q[0];
sx q[0];
rz(2.3556019) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(-2.5075066) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2848628) q[0];
sx q[0];
rz(-0.75548178) q[0];
sx q[0];
rz(-1.5289686) q[0];
rz(-pi) q[1];
rz(2.2251031) q[2];
sx q[2];
rz(-1.620703) q[2];
sx q[2];
rz(0.33545845) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12373911) q[1];
sx q[1];
rz(-2.1542319) q[1];
sx q[1];
rz(0.53637335) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5473458) q[3];
sx q[3];
rz(-1.1707077) q[3];
sx q[3];
rz(1.401702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9156076) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(1.263164) q[2];
rz(1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(-0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(2.7040226) q[0];
rz(-0.63105398) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(2.9204869) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3723345) q[0];
sx q[0];
rz(-2.6466469) q[0];
sx q[0];
rz(1.9702205) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2220076) q[2];
sx q[2];
rz(-1.4663327) q[2];
sx q[2];
rz(-1.8609357) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21048966) q[1];
sx q[1];
rz(-1.1138798) q[1];
sx q[1];
rz(-1.7060075) q[1];
x q[2];
rz(2.9651871) q[3];
sx q[3];
rz(-0.50900148) q[3];
sx q[3];
rz(2.727946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9235886) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(2.5533) q[2];
rz(-2.6925987) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(-2.0935521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5730729) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(-2.1411238) q[0];
rz(-0.72552848) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(-2.3838938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43514566) q[0];
sx q[0];
rz(-0.62951127) q[0];
sx q[0];
rz(-2.5802617) q[0];
rz(0.20747848) q[2];
sx q[2];
rz(-1.4138655) q[2];
sx q[2];
rz(-2.4776138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.12049) q[1];
sx q[1];
rz(-2.0334525) q[1];
sx q[1];
rz(0.91336577) q[1];
rz(-1.682231) q[3];
sx q[3];
rz(-2.206344) q[3];
sx q[3];
rz(-3.0524658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2086601) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(-2.3051252) q[2];
rz(1.6992735) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02012415) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(-2.6304723) q[0];
rz(3.064149) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(1.3285332) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7610361) q[0];
sx q[0];
rz(-1.7839) q[0];
sx q[0];
rz(-2.9831432) q[0];
rz(-pi) q[1];
rz(-2.0650234) q[2];
sx q[2];
rz(-0.87072125) q[2];
sx q[2];
rz(1.1717403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.37175825) q[1];
sx q[1];
rz(-1.3039939) q[1];
sx q[1];
rz(-0.00053243551) q[1];
x q[2];
rz(-1.1174326) q[3];
sx q[3];
rz(-2.8957267) q[3];
sx q[3];
rz(-2.9256431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0478583) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(1.1070586) q[2];
rz(-2.6691061) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(-0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040314019) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(-2.8033946) q[0];
rz(1.8473373) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(0.50450528) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.550068) q[0];
sx q[0];
rz(-1.803777) q[0];
sx q[0];
rz(-2.1165119) q[0];
rz(2.4733876) q[2];
sx q[2];
rz(-2.5720111) q[2];
sx q[2];
rz(2.6936206) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.63197631) q[1];
sx q[1];
rz(-2.7512433) q[1];
sx q[1];
rz(0.901464) q[1];
rz(-pi) q[2];
rz(-0.65977804) q[3];
sx q[3];
rz(-1.0228844) q[3];
sx q[3];
rz(-1.9632615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.111104) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(-1.4939235) q[2];
rz(1.6882287) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21022739) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(1.0908303) q[0];
rz(2.6018654) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(-2.9528023) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15121962) q[0];
sx q[0];
rz(-1.2114721) q[0];
sx q[0];
rz(-1.8381401) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0374755) q[2];
sx q[2];
rz(-1.3008442) q[2];
sx q[2];
rz(1.6402668) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.61154762) q[1];
sx q[1];
rz(-1.5049184) q[1];
sx q[1];
rz(-2.8987315) q[1];
rz(-2.4981899) q[3];
sx q[3];
rz(-2.8426369) q[3];
sx q[3];
rz(1.9770196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.171689) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(-1.8264654) q[2];
rz(-1.6361489) q[3];
sx q[3];
rz(-0.35741487) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0514907) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(-2.8175957) q[0];
rz(-1.8404768) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(1.3791929) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.798511) q[0];
sx q[0];
rz(-1.8146975) q[0];
sx q[0];
rz(0.084246158) q[0];
rz(1.9786644) q[2];
sx q[2];
rz(-1.1424354) q[2];
sx q[2];
rz(-0.63901627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.56475793) q[1];
sx q[1];
rz(-1.733629) q[1];
sx q[1];
rz(-2.8388883) q[1];
rz(-pi) q[2];
rz(0.93382436) q[3];
sx q[3];
rz(-2.056042) q[3];
sx q[3];
rz(-2.3683734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1088915) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(-0.94318715) q[2];
rz(2.8105248) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(-0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11809764) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(0.76422894) q[0];
rz(3.0006192) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.2932628) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4828669) q[0];
sx q[0];
rz(-1.7687161) q[0];
sx q[0];
rz(-2.0107962) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4752611) q[2];
sx q[2];
rz(-0.97180688) q[2];
sx q[2];
rz(2.4107188) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5871208) q[1];
sx q[1];
rz(-1.9738102) q[1];
sx q[1];
rz(-2.9906669) q[1];
rz(-pi) q[2];
rz(-0.39956283) q[3];
sx q[3];
rz(-1.385653) q[3];
sx q[3];
rz(-3.1294587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3747037) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-0.58132201) q[2];
rz(-0.86822048) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54810369) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(1.8909489) q[0];
rz(0.99682322) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.8539799) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33289136) q[0];
sx q[0];
rz(-1.5960072) q[0];
sx q[0];
rz(-1.0813792) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31918819) q[2];
sx q[2];
rz(-0.50779283) q[2];
sx q[2];
rz(2.9920981) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77546706) q[1];
sx q[1];
rz(-0.93535103) q[1];
sx q[1];
rz(0.25823621) q[1];
x q[2];
rz(-0.052863315) q[3];
sx q[3];
rz(-2.3164146) q[3];
sx q[3];
rz(2.4406274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.320497) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(1.8851177) q[2];
rz(-0.16658941) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(-2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2607516) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(-0.32521954) q[0];
rz(1.1351769) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(2.7744055) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8242278) q[0];
sx q[0];
rz(-1.5396848) q[0];
sx q[0];
rz(-0.048006417) q[0];
rz(-pi) q[1];
rz(-2.7634002) q[2];
sx q[2];
rz(-0.66236712) q[2];
sx q[2];
rz(-2.3364002) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.90875188) q[1];
sx q[1];
rz(-2.4526261) q[1];
sx q[1];
rz(0.31713756) q[1];
rz(-0.96079798) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8964768) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(1.0661351) q[2];
rz(-3.0623479) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8463678) q[0];
sx q[0];
rz(-1.5548779) q[0];
sx q[0];
rz(1.5665733) q[0];
rz(-0.29905839) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(2.3921641) q[2];
sx q[2];
rz(-1.6288169) q[2];
sx q[2];
rz(-1.5734869) q[2];
rz(-2.9885837) q[3];
sx q[3];
rz(-0.85481337) q[3];
sx q[3];
rz(-0.92683642) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];