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
rz(-1.4671833) q[0];
sx q[0];
rz(-1.6011342) q[0];
sx q[0];
rz(1.849548) q[0];
rz(2.0436824) q[1];
sx q[1];
rz(-0.75610375) q[1];
sx q[1];
rz(2.4296711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8956937) q[0];
sx q[0];
rz(-0.14764365) q[0];
sx q[0];
rz(2.9013322) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8991382) q[2];
sx q[2];
rz(-2.3676845) q[2];
sx q[2];
rz(2.7737308) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8375875) q[1];
sx q[1];
rz(-1.3437573) q[1];
sx q[1];
rz(-0.63804896) q[1];
rz(-pi) q[2];
rz(-0.96613066) q[3];
sx q[3];
rz(-1.1388766) q[3];
sx q[3];
rz(-2.8063065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0277675) q[2];
sx q[2];
rz(-1.0372459) q[2];
sx q[2];
rz(2.2859331) q[2];
rz(-1.8850108) q[3];
sx q[3];
rz(-0.62658566) q[3];
sx q[3];
rz(1.9470866) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41334316) q[0];
sx q[0];
rz(-3.0276868) q[0];
sx q[0];
rz(0.86694992) q[0];
rz(1.4504245) q[1];
sx q[1];
rz(-1.9536641) q[1];
sx q[1];
rz(-0.11122045) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040505458) q[0];
sx q[0];
rz(-1.2049647) q[0];
sx q[0];
rz(2.5554197) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1028981) q[2];
sx q[2];
rz(-1.1490029) q[2];
sx q[2];
rz(-1.8820764) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3635369) q[1];
sx q[1];
rz(-0.95603525) q[1];
sx q[1];
rz(0.78301439) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7584959) q[3];
sx q[3];
rz(-0.71569136) q[3];
sx q[3];
rz(2.7445499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0337246) q[2];
sx q[2];
rz(-2.3794231) q[2];
sx q[2];
rz(-2.6503358) q[2];
rz(-1.2563502) q[3];
sx q[3];
rz(-1.4278853) q[3];
sx q[3];
rz(0.45891941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1678109) q[0];
sx q[0];
rz(-1.351492) q[0];
sx q[0];
rz(-3.1381881) q[0];
rz(-2.7963474) q[1];
sx q[1];
rz(-2.7216941) q[1];
sx q[1];
rz(-0.87510625) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7686979) q[0];
sx q[0];
rz(-0.56124306) q[0];
sx q[0];
rz(0.13613693) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98710635) q[2];
sx q[2];
rz(-2.9265227) q[2];
sx q[2];
rz(-1.849086) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.7516084) q[1];
sx q[1];
rz(-0.64629236) q[1];
sx q[1];
rz(2.3306952) q[1];
rz(-0.77121021) q[3];
sx q[3];
rz(-1.5516087) q[3];
sx q[3];
rz(2.8018708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.7615937) q[2];
sx q[2];
rz(-0.21445175) q[2];
sx q[2];
rz(-2.8144515) q[2];
rz(1.3867779) q[3];
sx q[3];
rz(-1.1969457) q[3];
sx q[3];
rz(3.0676214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73579329) q[0];
sx q[0];
rz(-2.1402335) q[0];
sx q[0];
rz(-1.5843947) q[0];
rz(0.34218511) q[1];
sx q[1];
rz(-1.0289959) q[1];
sx q[1];
rz(2.7142966) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34278579) q[0];
sx q[0];
rz(-1.6248934) q[0];
sx q[0];
rz(-1.537772) q[0];
x q[1];
rz(0.9924381) q[2];
sx q[2];
rz(-2.522744) q[2];
sx q[2];
rz(-0.98471314) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7088741) q[1];
sx q[1];
rz(-2.1418835) q[1];
sx q[1];
rz(2.4576748) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5749608) q[3];
sx q[3];
rz(-1.8265836) q[3];
sx q[3];
rz(-1.2454288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7067318) q[2];
sx q[2];
rz(-2.2123983) q[2];
sx q[2];
rz(-0.84563196) q[2];
rz(-3.1137858) q[3];
sx q[3];
rz(-1.3549201) q[3];
sx q[3];
rz(-1.6644679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
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
rz(2.528275) q[0];
sx q[0];
rz(-2.1149825) q[0];
sx q[0];
rz(3.0897621) q[0];
rz(1.1593646) q[1];
sx q[1];
rz(-2.4720188) q[1];
sx q[1];
rz(2.5313098) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51175115) q[0];
sx q[0];
rz(-1.399222) q[0];
sx q[0];
rz(-2.9338994) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29783037) q[2];
sx q[2];
rz(-1.1874842) q[2];
sx q[2];
rz(-0.9155067) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7597639) q[1];
sx q[1];
rz(-2.0501158) q[1];
sx q[1];
rz(-1.3717531) q[1];
x q[2];
rz(1.5984437) q[3];
sx q[3];
rz(-0.96284722) q[3];
sx q[3];
rz(0.77800084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.74756885) q[2];
sx q[2];
rz(-1.3538066) q[2];
sx q[2];
rz(2.5547011) q[2];
rz(-1.3823357) q[3];
sx q[3];
rz(-1.047225) q[3];
sx q[3];
rz(-0.97572485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0368283) q[0];
sx q[0];
rz(-2.6703175) q[0];
sx q[0];
rz(-0.18071827) q[0];
rz(-1.7432927) q[1];
sx q[1];
rz(-1.9075874) q[1];
sx q[1];
rz(-1.3999375) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0796868) q[0];
sx q[0];
rz(-0.8041412) q[0];
sx q[0];
rz(1.2595533) q[0];
x q[1];
rz(0.39269523) q[2];
sx q[2];
rz(-1.942607) q[2];
sx q[2];
rz(-1.2758881) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9638979) q[1];
sx q[1];
rz(-0.70353466) q[1];
sx q[1];
rz(1.4369318) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65989699) q[3];
sx q[3];
rz(-2.0273034) q[3];
sx q[3];
rz(2.5259301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3900628) q[2];
sx q[2];
rz(-2.2921102) q[2];
sx q[2];
rz(2.0709822) q[2];
rz(1.6603598) q[3];
sx q[3];
rz(-2.345572) q[3];
sx q[3];
rz(-1.3720366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2299131) q[0];
sx q[0];
rz(-0.80753082) q[0];
sx q[0];
rz(2.5813778) q[0];
rz(-2.920976) q[1];
sx q[1];
rz(-0.95181528) q[1];
sx q[1];
rz(-2.2645948) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3426039) q[0];
sx q[0];
rz(-1.5734125) q[0];
sx q[0];
rz(-1.553234) q[0];
x q[1];
rz(-2.5602402) q[2];
sx q[2];
rz(-0.46581163) q[2];
sx q[2];
rz(-1.1752446) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0045053) q[1];
sx q[1];
rz(-1.3045038) q[1];
sx q[1];
rz(-3.0354663) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0362066) q[3];
sx q[3];
rz(-0.8871173) q[3];
sx q[3];
rz(-0.87382853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.47473869) q[2];
sx q[2];
rz(-2.2467504) q[2];
sx q[2];
rz(-1.8043) q[2];
rz(-0.21833359) q[3];
sx q[3];
rz(-2.1209013) q[3];
sx q[3];
rz(-1.7201503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6364994) q[0];
sx q[0];
rz(-1.8890843) q[0];
sx q[0];
rz(0.1388347) q[0];
rz(-0.88169634) q[1];
sx q[1];
rz(-1.9813184) q[1];
sx q[1];
rz(-0.82463157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3809842) q[0];
sx q[0];
rz(-3.0467817) q[0];
sx q[0];
rz(-1.8065971) q[0];
rz(-pi) q[1];
rz(3.0443848) q[2];
sx q[2];
rz(-1.4667454) q[2];
sx q[2];
rz(1.0417633) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6528553) q[1];
sx q[1];
rz(-1.0751546) q[1];
sx q[1];
rz(1.7486217) q[1];
rz(-pi) q[2];
rz(-2.2240198) q[3];
sx q[3];
rz(-1.6503449) q[3];
sx q[3];
rz(0.099640007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3489939) q[2];
sx q[2];
rz(-0.94453347) q[2];
sx q[2];
rz(-0.83354956) q[2];
rz(2.870657) q[3];
sx q[3];
rz(-2.0201594) q[3];
sx q[3];
rz(2.2942395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18868294) q[0];
sx q[0];
rz(-1.2485349) q[0];
sx q[0];
rz(-0.46440014) q[0];
rz(2.4480827) q[1];
sx q[1];
rz(-0.92718569) q[1];
sx q[1];
rz(-0.69673353) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9342142) q[0];
sx q[0];
rz(-1.6262282) q[0];
sx q[0];
rz(1.4595281) q[0];
rz(-pi) q[1];
rz(-1.7510173) q[2];
sx q[2];
rz(-1.5308793) q[2];
sx q[2];
rz(-0.60933622) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2285959) q[1];
sx q[1];
rz(-1.3218398) q[1];
sx q[1];
rz(3.1133786) q[1];
x q[2];
rz(0.018432004) q[3];
sx q[3];
rz(-0.64041172) q[3];
sx q[3];
rz(0.7225534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.36737475) q[2];
sx q[2];
rz(-3.0493224) q[2];
sx q[2];
rz(-2.3027159) q[2];
rz(-1.687626) q[3];
sx q[3];
rz(-1.5435012) q[3];
sx q[3];
rz(-1.4723697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0812747) q[0];
sx q[0];
rz(-1.769279) q[0];
sx q[0];
rz(1.6171612) q[0];
rz(-2.7853107) q[1];
sx q[1];
rz(-1.7226912) q[1];
sx q[1];
rz(1.3391395) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.08218) q[0];
sx q[0];
rz(-2.1689183) q[0];
sx q[0];
rz(-0.077738387) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21276591) q[2];
sx q[2];
rz(-0.32684718) q[2];
sx q[2];
rz(-0.67732993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.061627541) q[1];
sx q[1];
rz(-1.6929757) q[1];
sx q[1];
rz(-2.4604718) q[1];
rz(-pi) q[2];
x q[2];
rz(0.036063866) q[3];
sx q[3];
rz(-1.4286005) q[3];
sx q[3];
rz(-2.9309931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39849207) q[2];
sx q[2];
rz(-2.3091381) q[2];
sx q[2];
rz(-1.2644281) q[2];
rz(-2.8049331) q[3];
sx q[3];
rz(-2.2008937) q[3];
sx q[3];
rz(-1.8077067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3867415) q[0];
sx q[0];
rz(-1.6254397) q[0];
sx q[0];
rz(2.0425015) q[0];
rz(-2.5418538) q[1];
sx q[1];
rz(-0.44282423) q[1];
sx q[1];
rz(2.5437358) q[1];
rz(-2.2221634) q[2];
sx q[2];
rz(-2.6391657) q[2];
sx q[2];
rz(0.68965672) q[2];
rz(-0.029765263) q[3];
sx q[3];
rz(-0.53277848) q[3];
sx q[3];
rz(3.0938704) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
