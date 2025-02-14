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
rz(0.8910203) q[0];
sx q[0];
rz(1.855259) q[0];
sx q[0];
rz(8.5387908) q[0];
rz(0.52357829) q[1];
sx q[1];
rz(3.7012586) q[1];
sx q[1];
rz(8.1044365) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1170066) q[0];
sx q[0];
rz(-1.3947233) q[0];
sx q[0];
rz(-1.0278637) q[0];
rz(-pi) q[1];
rz(0.77729887) q[2];
sx q[2];
rz(-1.3242928) q[2];
sx q[2];
rz(2.8706626) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.953254) q[1];
sx q[1];
rz(-1.0085114) q[1];
sx q[1];
rz(3.0694993) q[1];
rz(-pi) q[2];
rz(2.376972) q[3];
sx q[3];
rz(-2.21661) q[3];
sx q[3];
rz(-2.3117206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.77871275) q[2];
sx q[2];
rz(-1.820463) q[2];
sx q[2];
rz(-2.2134181) q[2];
rz(1.5859531) q[3];
sx q[3];
rz(-1.0471683) q[3];
sx q[3];
rz(1.9716523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8856119) q[0];
sx q[0];
rz(-2.8035127) q[0];
sx q[0];
rz(-2.0804491) q[0];
rz(1.2127009) q[1];
sx q[1];
rz(-2.3586528) q[1];
sx q[1];
rz(-1.8276259) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82147017) q[0];
sx q[0];
rz(-0.87813963) q[0];
sx q[0];
rz(-2.0351324) q[0];
rz(-pi) q[1];
rz(1.018386) q[2];
sx q[2];
rz(-0.72478849) q[2];
sx q[2];
rz(-1.2040491) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.81345859) q[1];
sx q[1];
rz(-2.5812006) q[1];
sx q[1];
rz(1.8916393) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84557326) q[3];
sx q[3];
rz(-1.9872287) q[3];
sx q[3];
rz(-1.4786947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5709915) q[2];
sx q[2];
rz(-2.1676895) q[2];
sx q[2];
rz(2.2368597) q[2];
rz(-1.5500801) q[3];
sx q[3];
rz(-1.2856057) q[3];
sx q[3];
rz(1.2929644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078481361) q[0];
sx q[0];
rz(-0.08992973) q[0];
sx q[0];
rz(0.91823804) q[0];
rz(0.69084424) q[1];
sx q[1];
rz(-1.7450688) q[1];
sx q[1];
rz(-2.3450559) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7920216) q[0];
sx q[0];
rz(-0.22788985) q[0];
sx q[0];
rz(1.1538651) q[0];
rz(-1.630322) q[2];
sx q[2];
rz(-1.6135441) q[2];
sx q[2];
rz(1.8190365) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2921089) q[1];
sx q[1];
rz(-2.4359018) q[1];
sx q[1];
rz(1.4035426) q[1];
rz(-pi) q[2];
rz(-0.91751601) q[3];
sx q[3];
rz(-1.2923354) q[3];
sx q[3];
rz(2.3115932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1306661) q[2];
sx q[2];
rz(-0.9504168) q[2];
sx q[2];
rz(-1.2476904) q[2];
rz(2.583336) q[3];
sx q[3];
rz(-1.6853251) q[3];
sx q[3];
rz(-3.1202417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.438544) q[0];
sx q[0];
rz(-3.092364) q[0];
sx q[0];
rz(-2.4141648) q[0];
rz(0.1618596) q[1];
sx q[1];
rz(-1.1715803) q[1];
sx q[1];
rz(-1.1579827) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8901871) q[0];
sx q[0];
rz(-1.6368388) q[0];
sx q[0];
rz(-1.4736045) q[0];
rz(1.7100934) q[2];
sx q[2];
rz(-1.5717531) q[2];
sx q[2];
rz(0.19569163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47745142) q[1];
sx q[1];
rz(-2.913457) q[1];
sx q[1];
rz(-0.99598187) q[1];
x q[2];
rz(2.999008) q[3];
sx q[3];
rz(-0.53053427) q[3];
sx q[3];
rz(-0.50850463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8689279) q[2];
sx q[2];
rz(-1.9395892) q[2];
sx q[2];
rz(-1.391601) q[2];
rz(-2.9113655) q[3];
sx q[3];
rz(-1.9672491) q[3];
sx q[3];
rz(0.19190425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72628438) q[0];
sx q[0];
rz(-3.0400161) q[0];
sx q[0];
rz(-1.104785) q[0];
rz(-3.0305908) q[1];
sx q[1];
rz(-2.0260725) q[1];
sx q[1];
rz(1.312779) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3462145) q[0];
sx q[0];
rz(-1.7722478) q[0];
sx q[0];
rz(-0.18152118) q[0];
x q[1];
rz(2.5778002) q[2];
sx q[2];
rz(-1.6456266) q[2];
sx q[2];
rz(-2.5397391) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5178419) q[1];
sx q[1];
rz(-1.1718796) q[1];
sx q[1];
rz(-0.2174045) q[1];
rz(-pi) q[2];
rz(-2.2897374) q[3];
sx q[3];
rz(-1.2491711) q[3];
sx q[3];
rz(2.3648398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.63867265) q[2];
sx q[2];
rz(-1.0588249) q[2];
sx q[2];
rz(-2.1144833) q[2];
rz(0.98661679) q[3];
sx q[3];
rz(-0.28237453) q[3];
sx q[3];
rz(1.7178887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33243772) q[0];
sx q[0];
rz(-1.9092535) q[0];
sx q[0];
rz(0.80068457) q[0];
rz(2.370131) q[1];
sx q[1];
rz(-1.0779251) q[1];
sx q[1];
rz(1.7652184) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0906935) q[0];
sx q[0];
rz(-0.92779033) q[0];
sx q[0];
rz(-2.0741295) q[0];
rz(-pi) q[1];
rz(3.1065953) q[2];
sx q[2];
rz(-1.1049602) q[2];
sx q[2];
rz(2.310487) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1073947) q[1];
sx q[1];
rz(-1.5575111) q[1];
sx q[1];
rz(-2.4360869) q[1];
rz(-pi) q[2];
rz(0.8326859) q[3];
sx q[3];
rz(-2.247274) q[3];
sx q[3];
rz(3.0665086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8419522) q[2];
sx q[2];
rz(-1.2066634) q[2];
sx q[2];
rz(-0.017814962) q[2];
rz(0.31339112) q[3];
sx q[3];
rz(-2.119795) q[3];
sx q[3];
rz(2.4221086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7844836) q[0];
sx q[0];
rz(-1.5672368) q[0];
sx q[0];
rz(2.9743279) q[0];
rz(-0.74370614) q[1];
sx q[1];
rz(-0.88631648) q[1];
sx q[1];
rz(3.0053265) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24700227) q[0];
sx q[0];
rz(-0.73406363) q[0];
sx q[0];
rz(-1.1511001) q[0];
rz(-pi) q[1];
rz(1.874204) q[2];
sx q[2];
rz(-1.4055411) q[2];
sx q[2];
rz(-0.39656933) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4140461) q[1];
sx q[1];
rz(-0.36484499) q[1];
sx q[1];
rz(-0.18864297) q[1];
rz(-2.9877671) q[3];
sx q[3];
rz(-0.3496146) q[3];
sx q[3];
rz(-1.6817917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3027489) q[2];
sx q[2];
rz(-0.15041298) q[2];
sx q[2];
rz(-1.0813084) q[2];
rz(-0.14351621) q[3];
sx q[3];
rz(-1.84294) q[3];
sx q[3];
rz(1.6941841) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3926587) q[0];
sx q[0];
rz(-1.3329788) q[0];
sx q[0];
rz(0.64312154) q[0];
rz(2.2195623) q[1];
sx q[1];
rz(-2.8755201) q[1];
sx q[1];
rz(-1.5111074) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5095172) q[0];
sx q[0];
rz(-0.38783535) q[0];
sx q[0];
rz(-2.0616421) q[0];
x q[1];
rz(-3.0273075) q[2];
sx q[2];
rz(-1.292789) q[2];
sx q[2];
rz(2.4604527) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0550168) q[1];
sx q[1];
rz(-0.19682238) q[1];
sx q[1];
rz(1.6114525) q[1];
rz(-pi) q[2];
rz(-2.6073205) q[3];
sx q[3];
rz(-2.5639682) q[3];
sx q[3];
rz(-2.3848611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.94613218) q[2];
sx q[2];
rz(-1.5617424) q[2];
sx q[2];
rz(1.3004318) q[2];
rz(-2.2293633) q[3];
sx q[3];
rz(-1.2765086) q[3];
sx q[3];
rz(-2.8279878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4881956) q[0];
sx q[0];
rz(-2.6618239) q[0];
sx q[0];
rz(-2.524014) q[0];
rz(0.081534475) q[1];
sx q[1];
rz(-0.50059861) q[1];
sx q[1];
rz(0.68731442) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0162553) q[0];
sx q[0];
rz(-0.24152405) q[0];
sx q[0];
rz(-1.7707654) q[0];
rz(-pi) q[1];
rz(2.7527037) q[2];
sx q[2];
rz(-0.28770471) q[2];
sx q[2];
rz(-0.012118169) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0159434) q[1];
sx q[1];
rz(-1.9852891) q[1];
sx q[1];
rz(-1.8806367) q[1];
rz(3.0932759) q[3];
sx q[3];
rz(-1.028065) q[3];
sx q[3];
rz(0.84191546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.11757892) q[2];
sx q[2];
rz(-1.1774096) q[2];
sx q[2];
rz(0.071361072) q[2];
rz(1.2355545) q[3];
sx q[3];
rz(-2.8046298) q[3];
sx q[3];
rz(-0.56628847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.5679034) q[0];
sx q[0];
rz(-0.25147831) q[0];
sx q[0];
rz(0.13791826) q[0];
rz(0.57016405) q[1];
sx q[1];
rz(-2.0591996) q[1];
sx q[1];
rz(3.1261442) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041067657) q[0];
sx q[0];
rz(-2.3599632) q[0];
sx q[0];
rz(2.0360721) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49039109) q[2];
sx q[2];
rz(-2.3182959) q[2];
sx q[2];
rz(3.1237912) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8739104) q[1];
sx q[1];
rz(-1.7355669) q[1];
sx q[1];
rz(0.49652892) q[1];
x q[2];
rz(-1.8633719) q[3];
sx q[3];
rz(-0.19147123) q[3];
sx q[3];
rz(-1.1818238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7323759) q[2];
sx q[2];
rz(-2.3541383) q[2];
sx q[2];
rz(0.43238762) q[2];
rz(0.89808291) q[3];
sx q[3];
rz(-0.87528527) q[3];
sx q[3];
rz(1.557365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4846004) q[0];
sx q[0];
rz(-1.1285755) q[0];
sx q[0];
rz(-2.3791671) q[0];
rz(0.55105974) q[1];
sx q[1];
rz(-1.8012128) q[1];
sx q[1];
rz(2.6691379) q[1];
rz(2.77825) q[2];
sx q[2];
rz(-1.6132334) q[2];
sx q[2];
rz(1.2237534) q[2];
rz(-0.95672204) q[3];
sx q[3];
rz(-1.1720368) q[3];
sx q[3];
rz(-1.9839245) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
