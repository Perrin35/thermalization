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
rz(-0.78599077) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(0.63408607) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91416042) q[0];
sx q[0];
rz(-2.3254546) q[0];
sx q[0];
rz(0.039365191) q[0];
x q[1];
rz(0.91648957) q[2];
sx q[2];
rz(-1.620703) q[2];
sx q[2];
rz(-0.33545845) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0178535) q[1];
sx q[1];
rz(-2.1542319) q[1];
sx q[1];
rz(-2.6052193) q[1];
rz(-pi) q[2];
rz(-0.4001873) q[3];
sx q[3];
rz(-1.549198) q[3];
sx q[3];
rz(-0.17822972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2259851) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(-1.263164) q[2];
rz(-1.8566711) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(-2.7040226) q[0];
rz(-0.63105398) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(-2.9204869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3723345) q[0];
sx q[0];
rz(-0.49494574) q[0];
sx q[0];
rz(1.1713722) q[0];
rz(-pi) q[1];
rz(1.8684623) q[2];
sx q[2];
rz(-2.7781099) q[2];
sx q[2];
rz(0.01089451) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21048966) q[1];
sx q[1];
rz(-1.1138798) q[1];
sx q[1];
rz(-1.4355852) q[1];
rz(-pi) q[2];
rz(1.6684181) q[3];
sx q[3];
rz(-1.0704346) q[3];
sx q[3];
rz(-0.2122768) q[3];
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
rz(-0.58829266) q[2];
rz(-2.6925987) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(-2.0935521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5730729) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(1.0004689) q[0];
rz(-2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-2.3838938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.706447) q[0];
sx q[0];
rz(-0.62951127) q[0];
sx q[0];
rz(-2.5802617) q[0];
rz(-pi) q[1];
rz(1.4104841) q[2];
sx q[2];
rz(-1.3659039) q[2];
sx q[2];
rz(-0.93970539) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.11848395) q[1];
sx q[1];
rz(-2.1493836) q[1];
sx q[1];
rz(-2.5793377) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5030701) q[3];
sx q[3];
rz(-1.6604074) q[3];
sx q[3];
rz(-1.5935957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9329325) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(-2.3051252) q[2];
rz(1.4423192) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(-2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1214685) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(0.51112038) q[0];
rz(-3.064149) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(-1.8130594) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0275832) q[0];
sx q[0];
rz(-0.264835) q[0];
sx q[0];
rz(2.2008405) q[0];
rz(-pi) q[1];
rz(-0.76339108) q[2];
sx q[2];
rz(-1.1995458) q[2];
sx q[2];
rz(0.064918092) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.37175825) q[1];
sx q[1];
rz(-1.3039939) q[1];
sx q[1];
rz(-3.1410602) q[1];
rz(-pi) q[2];
rz(-1.3489181) q[3];
sx q[3];
rz(-1.6776049) q[3];
sx q[3];
rz(2.2282003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.093734309) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(2.034534) q[2];
rz(0.47248653) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012786) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(2.8033946) q[0];
rz(1.2942554) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(0.50450528) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.550068) q[0];
sx q[0];
rz(-1.3378157) q[0];
sx q[0];
rz(1.0250807) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.675823) q[2];
sx q[2];
rz(-1.911474) q[2];
sx q[2];
rz(-0.53616947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0658768) q[1];
sx q[1];
rz(-1.8738235) q[1];
sx q[1];
rz(-2.8916343) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65977804) q[3];
sx q[3];
rz(-2.1187083) q[3];
sx q[3];
rz(-1.1783311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.111104) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(1.4939235) q[2];
rz(-1.6882287) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21022739) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(1.0908303) q[0];
rz(-0.53972721) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(2.9528023) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8180346) q[0];
sx q[0];
rz(-1.8206882) q[0];
sx q[0];
rz(0.37139335) q[0];
rz(-1.2994453) q[2];
sx q[2];
rz(-1.6711298) q[2];
sx q[2];
rz(0.041610418) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1660359) q[1];
sx q[1];
rz(-1.3284725) q[1];
sx q[1];
rz(1.5029328) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24174989) q[3];
sx q[3];
rz(-1.7484192) q[3];
sx q[3];
rz(-2.1135981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9699036) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(1.3151273) q[2];
rz(-1.6361489) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.090102) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(-2.8175957) q[0];
rz(1.3011159) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(-1.7623998) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1349072) q[0];
sx q[0];
rz(-2.8838257) q[0];
sx q[0];
rz(-1.8968614) q[0];
x q[1];
rz(2.6799913) q[2];
sx q[2];
rz(-1.201655) q[2];
sx q[2];
rz(-2.387407) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.56475793) q[1];
sx q[1];
rz(-1.4079637) q[1];
sx q[1];
rz(0.30270438) q[1];
rz(0.58052766) q[3];
sx q[3];
rz(-2.1248098) q[3];
sx q[3];
rz(2.011727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1088915) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(-2.1984055) q[2];
rz(-0.33106783) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11809764) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(-0.76422894) q[0];
rz(0.14097342) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(1.8483298) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8338776) q[0];
sx q[0];
rz(-2.6617962) q[0];
sx q[0];
rz(2.0108372) q[0];
x q[1];
rz(-1.6663315) q[2];
sx q[2];
rz(-2.1697858) q[2];
sx q[2];
rz(0.73087382) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1845491) q[1];
sx q[1];
rz(-0.42889412) q[1];
sx q[1];
rz(1.2317608) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7420298) q[3];
sx q[3];
rz(-1.7559397) q[3];
sx q[3];
rz(3.1294587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3747037) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-2.5602706) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54810369) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(-1.2506437) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(1.8539799) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8087013) q[0];
sx q[0];
rz(-1.5960072) q[0];
sx q[0];
rz(-2.0602134) q[0];
x q[1];
rz(-0.31918819) q[2];
sx q[2];
rz(-2.6337998) q[2];
sx q[2];
rz(-2.9920981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3661256) q[1];
sx q[1];
rz(-2.2062416) q[1];
sx q[1];
rz(2.8833564) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.052863315) q[3];
sx q[3];
rz(-0.82517805) q[3];
sx q[3];
rz(0.70096522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.320497) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(1.256475) q[2];
rz(0.16658941) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88084108) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(0.32521954) q[0];
rz(1.1351769) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(2.7744055) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31736483) q[0];
sx q[0];
rz(-1.6019078) q[0];
sx q[0];
rz(3.0935862) q[0];
rz(1.2904097) q[2];
sx q[2];
rz(-2.1791611) q[2];
sx q[2];
rz(2.8031363) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.91017427) q[1];
sx q[1];
rz(-1.3712198) q[1];
sx q[1];
rz(0.66399666) q[1];
rz(-0.48093421) q[3];
sx q[3];
rz(-2.1554865) q[3];
sx q[3];
rz(-1.1993053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.24511589) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(1.0661351) q[2];
rz(0.079244763) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8463678) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(2.8425343) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(-2.3921641) q[2];
sx q[2];
rz(-1.5127758) q[2];
sx q[2];
rz(1.5681058) q[2];
rz(-0.15300898) q[3];
sx q[3];
rz(-2.2867793) q[3];
sx q[3];
rz(2.2147562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
