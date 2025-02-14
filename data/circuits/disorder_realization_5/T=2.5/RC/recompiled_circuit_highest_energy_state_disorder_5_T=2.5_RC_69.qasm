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
rz(-0.37762168) q[0];
sx q[0];
rz(3.5701516) q[0];
sx q[0];
rz(9.188434) q[0];
rz(-1.6726681) q[1];
sx q[1];
rz(-1.5863215) q[1];
sx q[1];
rz(2.9798887) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2387213) q[0];
sx q[0];
rz(-1.6797025) q[0];
sx q[0];
rz(-2.1542473) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5125844) q[2];
sx q[2];
rz(-1.3752642) q[2];
sx q[2];
rz(-0.19042507) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75085109) q[1];
sx q[1];
rz(-3.1138732) q[1];
sx q[1];
rz(-1.3185459) q[1];
x q[2];
rz(2.6246965) q[3];
sx q[3];
rz(-0.42498838) q[3];
sx q[3];
rz(1.4864511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24254313) q[2];
sx q[2];
rz(-0.87729064) q[2];
sx q[2];
rz(-2.7522932) q[2];
rz(0.23046514) q[3];
sx q[3];
rz(-3.1233628) q[3];
sx q[3];
rz(0.61483312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.56735754) q[0];
sx q[0];
rz(-2.1901972) q[0];
sx q[0];
rz(-1.4943328) q[0];
rz(1.5859454) q[1];
sx q[1];
rz(-2.9289398) q[1];
sx q[1];
rz(1.1344604) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0527549) q[0];
sx q[0];
rz(-2.2688534) q[0];
sx q[0];
rz(1.237117) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6642642) q[2];
sx q[2];
rz(-2.3466718) q[2];
sx q[2];
rz(1.3238465) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8061773) q[1];
sx q[1];
rz(-1.9796438) q[1];
sx q[1];
rz(-2.2061636) q[1];
x q[2];
rz(-1.0026917) q[3];
sx q[3];
rz(-0.8352931) q[3];
sx q[3];
rz(-0.42772528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0906585) q[2];
sx q[2];
rz(-2.2218573) q[2];
sx q[2];
rz(1.8402137) q[2];
rz(-2.0672412) q[3];
sx q[3];
rz(-2.8315872) q[3];
sx q[3];
rz(1.5955135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.018547) q[0];
sx q[0];
rz(-0.30236852) q[0];
sx q[0];
rz(2.5240335) q[0];
rz(-1.1005719) q[1];
sx q[1];
rz(-3.1220084) q[1];
sx q[1];
rz(2.7380131) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7779008) q[0];
sx q[0];
rz(-0.11523031) q[0];
sx q[0];
rz(-1.4794502) q[0];
rz(-pi) q[1];
rz(1.4694655) q[2];
sx q[2];
rz(-2.0943542) q[2];
sx q[2];
rz(-0.061274139) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.61768764) q[1];
sx q[1];
rz(-2.9343961) q[1];
sx q[1];
rz(2.2788384) q[1];
rz(-2.5837901) q[3];
sx q[3];
rz(-1.7770168) q[3];
sx q[3];
rz(0.6129515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5201716) q[2];
sx q[2];
rz(-1.6941864) q[2];
sx q[2];
rz(-2.3604895) q[2];
rz(2.805294) q[3];
sx q[3];
rz(-1.7016442) q[3];
sx q[3];
rz(0.70359105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4562562) q[0];
sx q[0];
rz(-2.6254613) q[0];
sx q[0];
rz(-1.2180895) q[0];
rz(-1.7958027) q[1];
sx q[1];
rz(-0.0082052611) q[1];
sx q[1];
rz(1.2340612) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06935519) q[0];
sx q[0];
rz(-1.580582) q[0];
sx q[0];
rz(0.18982498) q[0];
x q[1];
rz(-2.7379964) q[2];
sx q[2];
rz(-2.4880313) q[2];
sx q[2];
rz(0.60774481) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7750596) q[1];
sx q[1];
rz(-1.9294191) q[1];
sx q[1];
rz(2.8721362) q[1];
rz(-0.0099700516) q[3];
sx q[3];
rz(-1.5775351) q[3];
sx q[3];
rz(-0.88047319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9464843) q[2];
sx q[2];
rz(-0.4230963) q[2];
sx q[2];
rz(1.1484324) q[2];
rz(0.75442433) q[3];
sx q[3];
rz(-1.1634588) q[3];
sx q[3];
rz(-2.8783126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36963439) q[0];
sx q[0];
rz(-0.088005528) q[0];
sx q[0];
rz(-0.60246402) q[0];
rz(0.98214904) q[1];
sx q[1];
rz(-0.0063449675) q[1];
sx q[1];
rz(2.6314661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0786297) q[0];
sx q[0];
rz(-1.3853243) q[0];
sx q[0];
rz(-1.7679035) q[0];
x q[1];
rz(3.0544364) q[2];
sx q[2];
rz(-1.6917602) q[2];
sx q[2];
rz(2.6270514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.229631) q[1];
sx q[1];
rz(-1.0385286) q[1];
sx q[1];
rz(2.0144573) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3272133) q[3];
sx q[3];
rz(-1.571221) q[3];
sx q[3];
rz(-2.9546776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.043857418) q[2];
sx q[2];
rz(-1.1172224) q[2];
sx q[2];
rz(-1.0350234) q[2];
rz(-1.6739316) q[3];
sx q[3];
rz(-0.36164713) q[3];
sx q[3];
rz(-2.5586832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3007091) q[0];
sx q[0];
rz(-2.1568334) q[0];
sx q[0];
rz(-2.050052) q[0];
rz(-2.7409399) q[1];
sx q[1];
rz(-3.1392097) q[1];
sx q[1];
rz(0.56652743) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5015903) q[0];
sx q[0];
rz(-0.8336319) q[0];
sx q[0];
rz(-2.9596674) q[0];
rz(0.80711694) q[2];
sx q[2];
rz(-0.87701488) q[2];
sx q[2];
rz(-1.4556231) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8174584) q[1];
sx q[1];
rz(-2.5310881) q[1];
sx q[1];
rz(0.39892674) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9683365) q[3];
sx q[3];
rz(-0.57954407) q[3];
sx q[3];
rz(0.15523191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.94622508) q[2];
sx q[2];
rz(-0.75104284) q[2];
sx q[2];
rz(1.6132149) q[2];
rz(-2.1402806) q[3];
sx q[3];
rz(-0.87037218) q[3];
sx q[3];
rz(2.9010469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2840851) q[0];
sx q[0];
rz(-0.79428285) q[0];
sx q[0];
rz(1.9965782) q[0];
rz(-1.6192294) q[1];
sx q[1];
rz(-0.018773627) q[1];
sx q[1];
rz(-1.9782664) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49646851) q[0];
sx q[0];
rz(-1.2815406) q[0];
sx q[0];
rz(0.048286322) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3848035) q[2];
sx q[2];
rz(-2.459068) q[2];
sx q[2];
rz(-0.96772099) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.39966291) q[1];
sx q[1];
rz(-0.96782875) q[1];
sx q[1];
rz(-0.16027995) q[1];
rz(-2.9807253) q[3];
sx q[3];
rz(-0.64023521) q[3];
sx q[3];
rz(-1.36659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5160949) q[2];
sx q[2];
rz(-0.98573804) q[2];
sx q[2];
rz(0.37334785) q[2];
rz(0.18800023) q[3];
sx q[3];
rz(-1.5550273) q[3];
sx q[3];
rz(-1.9787623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2242551) q[0];
sx q[0];
rz(-0.52977109) q[0];
sx q[0];
rz(-0.24714558) q[0];
rz(-2.9180134) q[1];
sx q[1];
rz(-3.1378742) q[1];
sx q[1];
rz(-1.5162969) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7784268) q[0];
sx q[0];
rz(-1.4545621) q[0];
sx q[0];
rz(-3.0768125) q[0];
rz(-pi) q[1];
rz(0.76659492) q[2];
sx q[2];
rz(-1.8803036) q[2];
sx q[2];
rz(-0.37878894) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.030236472) q[1];
sx q[1];
rz(-1.9647536) q[1];
sx q[1];
rz(1.5880821) q[1];
rz(0.070930158) q[3];
sx q[3];
rz(-2.5462357) q[3];
sx q[3];
rz(-1.2687253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.94459263) q[2];
sx q[2];
rz(-2.4304515) q[2];
sx q[2];
rz(1.5947343) q[2];
rz(2.366015) q[3];
sx q[3];
rz(-1.2838793) q[3];
sx q[3];
rz(-1.4625134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.813886) q[0];
sx q[0];
rz(-2.1729108) q[0];
sx q[0];
rz(1.1073329) q[0];
rz(2.2164717) q[1];
sx q[1];
rz(-0.0019625891) q[1];
sx q[1];
rz(-0.75606871) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30413142) q[0];
sx q[0];
rz(-0.63429773) q[0];
sx q[0];
rz(-2.1055806) q[0];
x q[1];
rz(-1.0399516) q[2];
sx q[2];
rz(-1.5762946) q[2];
sx q[2];
rz(-2.0844242) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1518095) q[1];
sx q[1];
rz(-2.1266471) q[1];
sx q[1];
rz(-1.9152867) q[1];
rz(-pi) q[2];
rz(-0.57188923) q[3];
sx q[3];
rz(-1.6731305) q[3];
sx q[3];
rz(-2.1646647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5209311) q[2];
sx q[2];
rz(-2.2563917) q[2];
sx q[2];
rz(-2.1405641) q[2];
rz(1.3641317) q[3];
sx q[3];
rz(-0.93327156) q[3];
sx q[3];
rz(2.6076243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022024632) q[0];
sx q[0];
rz(-1.7689995) q[0];
sx q[0];
rz(0.46510988) q[0];
rz(1.8067092) q[1];
sx q[1];
rz(-2.7692134) q[1];
sx q[1];
rz(-1.575527) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4441285) q[0];
sx q[0];
rz(-0.23374548) q[0];
sx q[0];
rz(1.759619) q[0];
rz(-0.31120531) q[2];
sx q[2];
rz(-1.0619783) q[2];
sx q[2];
rz(-0.52647639) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.59908341) q[1];
sx q[1];
rz(-1.570382) q[1];
sx q[1];
rz(3.1387657) q[1];
rz(0.09327831) q[3];
sx q[3];
rz(-2.7607548) q[3];
sx q[3];
rz(-0.24254984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2472725) q[2];
sx q[2];
rz(-0.044581052) q[2];
sx q[2];
rz(-1.0856005) q[2];
rz(1.2224489) q[3];
sx q[3];
rz(-0.51210755) q[3];
sx q[3];
rz(1.2781757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4492252) q[0];
sx q[0];
rz(-1.7074371) q[0];
sx q[0];
rz(-1.3771124) q[0];
rz(-1.5832681) q[1];
sx q[1];
rz(-0.91455864) q[1];
sx q[1];
rz(0.22462489) q[1];
rz(1.5660607) q[2];
sx q[2];
rz(-1.4700459) q[2];
sx q[2];
rz(-2.8503304) q[2];
rz(2.0153981) q[3];
sx q[3];
rz(-1.7845407) q[3];
sx q[3];
rz(0.15872628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
