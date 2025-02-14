OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.4975118) q[0];
sx q[0];
rz(-1.8026135) q[0];
sx q[0];
rz(-2.8259377) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(2.7538731) q[1];
sx q[1];
rz(8.8207689) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7496628) q[0];
sx q[0];
rz(-1.6214341) q[0];
sx q[0];
rz(-1.2681566) q[0];
x q[1];
rz(2.1748105) q[2];
sx q[2];
rz(-0.76672115) q[2];
sx q[2];
rz(-1.3441835) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4872514) q[1];
sx q[1];
rz(-2.5120334) q[1];
sx q[1];
rz(0.97626026) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8995057) q[3];
sx q[3];
rz(-0.32530537) q[3];
sx q[3];
rz(-1.6360375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.151256) q[2];
sx q[2];
rz(-1.1948816) q[2];
sx q[2];
rz(-1.7830431) q[2];
rz(1.547706) q[3];
sx q[3];
rz(-0.93021506) q[3];
sx q[3];
rz(-1.6319298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58594054) q[0];
sx q[0];
rz(-0.63783115) q[0];
sx q[0];
rz(1.9441388) q[0];
rz(-2.6400631) q[1];
sx q[1];
rz(-1.5634368) q[1];
sx q[1];
rz(-1.4362358) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.013151) q[0];
sx q[0];
rz(-1.801898) q[0];
sx q[0];
rz(2.1559155) q[0];
x q[1];
rz(0.11709638) q[2];
sx q[2];
rz(-0.70636049) q[2];
sx q[2];
rz(1.7603086) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9402071) q[1];
sx q[1];
rz(-1.6734945) q[1];
sx q[1];
rz(2.9443355) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6858176) q[3];
sx q[3];
rz(-1.3823798) q[3];
sx q[3];
rz(-0.51793232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2531835) q[2];
sx q[2];
rz(-1.9056355) q[2];
sx q[2];
rz(0.02296981) q[2];
rz(-0.90211558) q[3];
sx q[3];
rz(-1.2227367) q[3];
sx q[3];
rz(2.7149916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4334634) q[0];
sx q[0];
rz(-1.7925649) q[0];
sx q[0];
rz(0.9915114) q[0];
rz(1.0333215) q[1];
sx q[1];
rz(-2.8585377) q[1];
sx q[1];
rz(-1.906377) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2184705) q[0];
sx q[0];
rz(-1.8184796) q[0];
sx q[0];
rz(-0.60486233) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60977817) q[2];
sx q[2];
rz(-1.1705361) q[2];
sx q[2];
rz(-2.5230809) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.24656235) q[1];
sx q[1];
rz(-0.8324648) q[1];
sx q[1];
rz(-2.8042996) q[1];
rz(-pi) q[2];
rz(-2.1622873) q[3];
sx q[3];
rz(-1.1011754) q[3];
sx q[3];
rz(-0.33794935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.11423763) q[2];
sx q[2];
rz(-2.3894252) q[2];
sx q[2];
rz(-0.99119622) q[2];
rz(1.8441955) q[3];
sx q[3];
rz(-1.2605366) q[3];
sx q[3];
rz(1.7245002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7436413) q[0];
sx q[0];
rz(-0.93248168) q[0];
sx q[0];
rz(-2.3241296) q[0];
rz(-2.815333) q[1];
sx q[1];
rz(-1.9605109) q[1];
sx q[1];
rz(2.356333) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9049587) q[0];
sx q[0];
rz(-1.3044622) q[0];
sx q[0];
rz(-1.8840709) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26008028) q[2];
sx q[2];
rz(-2.4818015) q[2];
sx q[2];
rz(-1.5957956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2438812) q[1];
sx q[1];
rz(-0.63279018) q[1];
sx q[1];
rz(-0.80676078) q[1];
x q[2];
rz(2.7866632) q[3];
sx q[3];
rz(-1.6454434) q[3];
sx q[3];
rz(1.1433034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.12104812) q[2];
sx q[2];
rz(-0.64457568) q[2];
sx q[2];
rz(-2.5784967) q[2];
rz(0.8935039) q[3];
sx q[3];
rz(-1.1734791) q[3];
sx q[3];
rz(1.9068498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23984443) q[0];
sx q[0];
rz(-0.32308602) q[0];
sx q[0];
rz(1.5455986) q[0];
rz(-2.1196938) q[1];
sx q[1];
rz(-1.1232168) q[1];
sx q[1];
rz(2.9603069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.793042) q[0];
sx q[0];
rz(-0.53213813) q[0];
sx q[0];
rz(-1.9572958) q[0];
rz(-pi) q[1];
rz(-3.0490169) q[2];
sx q[2];
rz(-2.6642054) q[2];
sx q[2];
rz(-0.30711781) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1039155) q[1];
sx q[1];
rz(-1.5427151) q[1];
sx q[1];
rz(-1.7930536) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4138172) q[3];
sx q[3];
rz(-2.008956) q[3];
sx q[3];
rz(2.639132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1229317) q[2];
sx q[2];
rz(-0.67492008) q[2];
sx q[2];
rz(-1.2459416) q[2];
rz(0.59257007) q[3];
sx q[3];
rz(-1.4453245) q[3];
sx q[3];
rz(-2.3004801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5331921) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(-2.8327508) q[0];
rz(-2.8187075) q[1];
sx q[1];
rz(-2.7643118) q[1];
sx q[1];
rz(0.13993941) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9062146) q[0];
sx q[0];
rz(-1.8188634) q[0];
sx q[0];
rz(2.1885022) q[0];
rz(-1.1148648) q[2];
sx q[2];
rz(-1.4856292) q[2];
sx q[2];
rz(0.87907253) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9819239) q[1];
sx q[1];
rz(-0.86729151) q[1];
sx q[1];
rz(2.7419006) q[1];
x q[2];
rz(1.9512964) q[3];
sx q[3];
rz(-0.69529136) q[3];
sx q[3];
rz(-2.8075346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66190019) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(-0.24197401) q[2];
rz(0.74609977) q[3];
sx q[3];
rz(-1.7753764) q[3];
sx q[3];
rz(1.7297176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.022843) q[0];
sx q[0];
rz(-1.043909) q[0];
sx q[0];
rz(1.2710849) q[0];
rz(0.56308833) q[1];
sx q[1];
rz(-0.92910281) q[1];
sx q[1];
rz(-1.7005327) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6800425) q[0];
sx q[0];
rz(-0.6667887) q[0];
sx q[0];
rz(2.0186485) q[0];
rz(-0.38132847) q[2];
sx q[2];
rz(-1.4421808) q[2];
sx q[2];
rz(2.6113457) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5477834) q[1];
sx q[1];
rz(-1.9873706) q[1];
sx q[1];
rz(-1.3273456) q[1];
x q[2];
rz(0.24347476) q[3];
sx q[3];
rz(-2.4015275) q[3];
sx q[3];
rz(-2.8748663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9219804) q[2];
sx q[2];
rz(-1.2954804) q[2];
sx q[2];
rz(-1.8249576) q[2];
rz(0.5611788) q[3];
sx q[3];
rz(-0.96446529) q[3];
sx q[3];
rz(1.5285899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.104326) q[0];
sx q[0];
rz(-0.21793652) q[0];
sx q[0];
rz(2.2547146) q[0];
rz(0.2001702) q[1];
sx q[1];
rz(-1.6693516) q[1];
sx q[1];
rz(2.1194469) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5258785) q[0];
sx q[0];
rz(-1.7019042) q[0];
sx q[0];
rz(-0.025630533) q[0];
rz(-pi) q[1];
rz(-2.3273507) q[2];
sx q[2];
rz(-2.1541284) q[2];
sx q[2];
rz(0.63177201) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31412582) q[1];
sx q[1];
rz(-0.74605251) q[1];
sx q[1];
rz(1.5681727) q[1];
rz(-0.57557801) q[3];
sx q[3];
rz(-1.8691571) q[3];
sx q[3];
rz(0.23250599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6518121) q[2];
sx q[2];
rz(-0.45280364) q[2];
sx q[2];
rz(-2.32302) q[2];
rz(-2.1583648) q[3];
sx q[3];
rz(-0.95663095) q[3];
sx q[3];
rz(-0.53808588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035456903) q[0];
sx q[0];
rz(-1.1470733) q[0];
sx q[0];
rz(-1.5863093) q[0];
rz(-2.4447794) q[1];
sx q[1];
rz(-2.9882444) q[1];
sx q[1];
rz(-2.3568025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9826384) q[0];
sx q[0];
rz(-1.0305911) q[0];
sx q[0];
rz(2.0870952) q[0];
rz(-pi) q[1];
x q[1];
rz(1.773514) q[2];
sx q[2];
rz(-1.7200816) q[2];
sx q[2];
rz(1.9307856) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1363298) q[1];
sx q[1];
rz(-1.6044093) q[1];
sx q[1];
rz(-0.43515794) q[1];
x q[2];
rz(2.574563) q[3];
sx q[3];
rz(-2.4483557) q[3];
sx q[3];
rz(2.6854613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67184225) q[2];
sx q[2];
rz(-1.2776813) q[2];
sx q[2];
rz(-0.81076852) q[2];
rz(1.9226711) q[3];
sx q[3];
rz(-0.54088497) q[3];
sx q[3];
rz(2.0764652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78275457) q[0];
sx q[0];
rz(-0.29958075) q[0];
sx q[0];
rz(1.7175571) q[0];
rz(2.3639823) q[1];
sx q[1];
rz(-2.168226) q[1];
sx q[1];
rz(-0.75540677) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41235218) q[0];
sx q[0];
rz(-1.5662114) q[0];
sx q[0];
rz(1.3608749) q[0];
rz(-pi) q[1];
rz(-2.7451186) q[2];
sx q[2];
rz(-1.5546397) q[2];
sx q[2];
rz(1.1262058) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3302119) q[1];
sx q[1];
rz(-0.63618681) q[1];
sx q[1];
rz(0.61417555) q[1];
rz(-pi) q[2];
rz(-0.67556583) q[3];
sx q[3];
rz(-1.8295975) q[3];
sx q[3];
rz(1.344556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5767673) q[2];
sx q[2];
rz(-0.29594031) q[2];
sx q[2];
rz(-0.15288615) q[2];
rz(1.2711924) q[3];
sx q[3];
rz(-1.8891687) q[3];
sx q[3];
rz(-0.66711867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85970238) q[0];
sx q[0];
rz(-0.49260456) q[0];
sx q[0];
rz(2.3717666) q[0];
rz(1.9181171) q[1];
sx q[1];
rz(-1.2212831) q[1];
sx q[1];
rz(2.4881359) q[1];
rz(-0.6003004) q[2];
sx q[2];
rz(-0.68927286) q[2];
sx q[2];
rz(2.2130028) q[2];
rz(2.2822904) q[3];
sx q[3];
rz(-2.6283384) q[3];
sx q[3];
rz(2.8859861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
