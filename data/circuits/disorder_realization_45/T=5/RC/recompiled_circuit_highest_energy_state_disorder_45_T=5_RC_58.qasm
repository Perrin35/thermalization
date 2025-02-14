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
rz(-1.4650605) q[0];
sx q[0];
rz(-0.2413916) q[0];
sx q[0];
rz(0.13225947) q[0];
rz(-0.013068696) q[1];
sx q[1];
rz(-0.71813923) q[1];
sx q[1];
rz(-0.018996039) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8762704) q[0];
sx q[0];
rz(-0.92060584) q[0];
sx q[0];
rz(0.27133915) q[0];
x q[1];
rz(-1.7573518) q[2];
sx q[2];
rz(-1.6636408) q[2];
sx q[2];
rz(1.2684737) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2680696) q[1];
sx q[1];
rz(-1.2482845) q[1];
sx q[1];
rz(1.7544062) q[1];
rz(-pi) q[2];
rz(1.8034987) q[3];
sx q[3];
rz(-1.4178172) q[3];
sx q[3];
rz(0.97319727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.12048177) q[2];
sx q[2];
rz(-2.3349473) q[2];
sx q[2];
rz(-2.3497904) q[2];
rz(1.0857438) q[3];
sx q[3];
rz(-1.2191685) q[3];
sx q[3];
rz(2.048548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2708112) q[0];
sx q[0];
rz(-2.9830611) q[0];
sx q[0];
rz(-0.32919163) q[0];
rz(-1.6393433) q[1];
sx q[1];
rz(-2.2396125) q[1];
sx q[1];
rz(0.42218581) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5188816) q[0];
sx q[0];
rz(-2.4269613) q[0];
sx q[0];
rz(-0.67784061) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5936826) q[2];
sx q[2];
rz(-1.7284677) q[2];
sx q[2];
rz(1.6239177) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.77743769) q[1];
sx q[1];
rz(-1.6171439) q[1];
sx q[1];
rz(-0.32407659) q[1];
rz(-1.6738191) q[3];
sx q[3];
rz(-2.3405082) q[3];
sx q[3];
rz(-2.5689606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4210356) q[2];
sx q[2];
rz(-1.1853508) q[2];
sx q[2];
rz(1.1400918) q[2];
rz(-0.0044048443) q[3];
sx q[3];
rz(-1.5811698) q[3];
sx q[3];
rz(-0.13979039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36188257) q[0];
sx q[0];
rz(-0.10040586) q[0];
sx q[0];
rz(-2.7959339) q[0];
rz(2.0480305) q[1];
sx q[1];
rz(-2.3193017) q[1];
sx q[1];
rz(2.0872769) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7436126) q[0];
sx q[0];
rz(-2.6184887) q[0];
sx q[0];
rz(2.4260957) q[0];
rz(-pi) q[1];
rz(-2.3286211) q[2];
sx q[2];
rz(-1.6457498) q[2];
sx q[2];
rz(-0.048150657) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.24040996) q[1];
sx q[1];
rz(-2.6844271) q[1];
sx q[1];
rz(2.9967876) q[1];
x q[2];
rz(-0.90462138) q[3];
sx q[3];
rz(-1.5525177) q[3];
sx q[3];
rz(-1.3410904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67636079) q[2];
sx q[2];
rz(-0.54764843) q[2];
sx q[2];
rz(2.5993627) q[2];
rz(0.17255653) q[3];
sx q[3];
rz(-1.4901284) q[3];
sx q[3];
rz(1.1057314) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3365823) q[0];
sx q[0];
rz(-1.7886826) q[0];
sx q[0];
rz(-1.2611058) q[0];
rz(1.1359967) q[1];
sx q[1];
rz(-2.0295862) q[1];
sx q[1];
rz(0.13066185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4379556) q[0];
sx q[0];
rz(-1.9398296) q[0];
sx q[0];
rz(-0.28965182) q[0];
rz(-pi) q[1];
rz(2.7515) q[2];
sx q[2];
rz(-2.60139) q[2];
sx q[2];
rz(-0.15946968) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4266721) q[1];
sx q[1];
rz(-1.0627278) q[1];
sx q[1];
rz(-3.1053468) q[1];
x q[2];
rz(0.66417905) q[3];
sx q[3];
rz(-1.8669493) q[3];
sx q[3];
rz(0.68063762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38264349) q[2];
sx q[2];
rz(-3.068277) q[2];
sx q[2];
rz(-2.9200413) q[2];
rz(-0.70212901) q[3];
sx q[3];
rz(-1.5626855) q[3];
sx q[3];
rz(0.73739541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5704983) q[0];
sx q[0];
rz(-1.2449188) q[0];
sx q[0];
rz(3.0076497) q[0];
rz(-0.36930034) q[1];
sx q[1];
rz(-1.9422653) q[1];
sx q[1];
rz(1.3901002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62890118) q[0];
sx q[0];
rz(-0.57399625) q[0];
sx q[0];
rz(-1.1074793) q[0];
rz(-1.8758135) q[2];
sx q[2];
rz(-1.8243363) q[2];
sx q[2];
rz(2.5407052) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0005324) q[1];
sx q[1];
rz(-0.62126011) q[1];
sx q[1];
rz(0.48807524) q[1];
x q[2];
rz(1.8769299) q[3];
sx q[3];
rz(-1.2226893) q[3];
sx q[3];
rz(2.8814684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2964581) q[2];
sx q[2];
rz(-0.90201169) q[2];
sx q[2];
rz(-1.585539) q[2];
rz(2.8499917) q[3];
sx q[3];
rz(-0.4762989) q[3];
sx q[3];
rz(0.27879032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6911102) q[0];
sx q[0];
rz(-0.99523681) q[0];
sx q[0];
rz(2.6158748) q[0];
rz(2.763343) q[1];
sx q[1];
rz(-1.5317761) q[1];
sx q[1];
rz(-0.27413109) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12488406) q[0];
sx q[0];
rz(-1.5450461) q[0];
sx q[0];
rz(-2.5514249) q[0];
rz(-pi) q[1];
rz(0.81118213) q[2];
sx q[2];
rz(-0.59181442) q[2];
sx q[2];
rz(2.9331911) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9410206) q[1];
sx q[1];
rz(-1.1686687) q[1];
sx q[1];
rz(0.27718039) q[1];
x q[2];
rz(1.317362) q[3];
sx q[3];
rz(-1.4259286) q[3];
sx q[3];
rz(2.1180017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5963001) q[2];
sx q[2];
rz(-1.9175074) q[2];
sx q[2];
rz(-0.25823414) q[2];
rz(-1.5457414) q[3];
sx q[3];
rz(-1.7786547) q[3];
sx q[3];
rz(-0.405092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.3426568) q[0];
sx q[0];
rz(-0.5181784) q[0];
sx q[0];
rz(0.10547353) q[0];
rz(2.8666829) q[1];
sx q[1];
rz(-1.4242438) q[1];
sx q[1];
rz(-0.28604937) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8559694) q[0];
sx q[0];
rz(-0.21846314) q[0];
sx q[0];
rz(-0.56665786) q[0];
rz(1.8931383) q[2];
sx q[2];
rz(-0.76248705) q[2];
sx q[2];
rz(2.40138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0512684) q[1];
sx q[1];
rz(-2.7503138) q[1];
sx q[1];
rz(1.7404187) q[1];
rz(-pi) q[2];
rz(0.6577095) q[3];
sx q[3];
rz(-2.4793787) q[3];
sx q[3];
rz(-2.7194104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1470571) q[2];
sx q[2];
rz(-2.4592082) q[2];
sx q[2];
rz(2.6915754) q[2];
rz(2.5543645) q[3];
sx q[3];
rz(-1.8000032) q[3];
sx q[3];
rz(-1.2500866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29379544) q[0];
sx q[0];
rz(-1.4601409) q[0];
sx q[0];
rz(-1.9548804) q[0];
rz(2.8222491) q[1];
sx q[1];
rz(-2.1198544) q[1];
sx q[1];
rz(-2.879338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05495435) q[0];
sx q[0];
rz(-2.1146746) q[0];
sx q[0];
rz(-1.8094814) q[0];
rz(1.2683343) q[2];
sx q[2];
rz(-2.710963) q[2];
sx q[2];
rz(2.9521049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9567392) q[1];
sx q[1];
rz(-1.0295086) q[1];
sx q[1];
rz(-0.095603099) q[1];
rz(0.45164177) q[3];
sx q[3];
rz(-1.9596976) q[3];
sx q[3];
rz(-1.9967784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.78476) q[2];
sx q[2];
rz(-0.46669745) q[2];
sx q[2];
rz(-2.812815) q[2];
rz(1.6262936) q[3];
sx q[3];
rz(-1.7833775) q[3];
sx q[3];
rz(2.7330107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3560155) q[0];
sx q[0];
rz(-0.51078904) q[0];
sx q[0];
rz(2.862634) q[0];
rz(-2.6246159) q[1];
sx q[1];
rz(-2.7644283) q[1];
sx q[1];
rz(-1.4252211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12869975) q[0];
sx q[0];
rz(-0.94966799) q[0];
sx q[0];
rz(1.6245882) q[0];
rz(0.80323146) q[2];
sx q[2];
rz(-1.0617501) q[2];
sx q[2];
rz(2.5521297) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9659075) q[1];
sx q[1];
rz(-1.07465) q[1];
sx q[1];
rz(1.3555075) q[1];
rz(-pi) q[2];
rz(-1.9402949) q[3];
sx q[3];
rz(-2.4505121) q[3];
sx q[3];
rz(1.406351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7384537) q[2];
sx q[2];
rz(-2.8454915) q[2];
sx q[2];
rz(2.9816755) q[2];
rz(-2.3219409) q[3];
sx q[3];
rz(-1.22217) q[3];
sx q[3];
rz(-2.405449) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4495471) q[0];
sx q[0];
rz(-1.8501546) q[0];
sx q[0];
rz(2.8458169) q[0];
rz(2.6462818) q[1];
sx q[1];
rz(-2.7073529) q[1];
sx q[1];
rz(1.6326509) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49182941) q[0];
sx q[0];
rz(-1.2163645) q[0];
sx q[0];
rz(1.9948122) q[0];
x q[1];
rz(3.0894439) q[2];
sx q[2];
rz(-1.9226769) q[2];
sx q[2];
rz(3.1089442) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1965643) q[1];
sx q[1];
rz(-1.5032282) q[1];
sx q[1];
rz(0.089971926) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4911777) q[3];
sx q[3];
rz(-0.59594184) q[3];
sx q[3];
rz(-0.015344674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32790023) q[2];
sx q[2];
rz(-1.4699961) q[2];
sx q[2];
rz(0.63391614) q[2];
rz(3.0564803) q[3];
sx q[3];
rz(-1.2462933) q[3];
sx q[3];
rz(-2.3688721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2278628) q[0];
sx q[0];
rz(-1.5416523) q[0];
sx q[0];
rz(-0.1097485) q[0];
rz(-1.9204503) q[1];
sx q[1];
rz(-2.2962062) q[1];
sx q[1];
rz(-2.5795945) q[1];
rz(0.22460266) q[2];
sx q[2];
rz(-1.1418268) q[2];
sx q[2];
rz(1.8662966) q[2];
rz(-1.802465) q[3];
sx q[3];
rz(-1.5944407) q[3];
sx q[3];
rz(-3.1138314) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
