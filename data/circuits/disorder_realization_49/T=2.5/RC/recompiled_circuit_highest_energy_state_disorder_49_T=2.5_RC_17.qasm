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
rz(0.024700392) q[0];
sx q[0];
rz(-1.2588809) q[0];
sx q[0];
rz(1.2727241) q[0];
rz(-4.89115) q[1];
sx q[1];
rz(1.7311544) q[1];
sx q[1];
rz(16.601736) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0426038) q[0];
sx q[0];
rz(-1.0013784) q[0];
sx q[0];
rz(0.75890394) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4364901) q[2];
sx q[2];
rz(-2.7573708) q[2];
sx q[2];
rz(2.9085858) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.53613816) q[1];
sx q[1];
rz(-1.0404603) q[1];
sx q[1];
rz(2.082389) q[1];
x q[2];
rz(0.36811604) q[3];
sx q[3];
rz(-1.5541346) q[3];
sx q[3];
rz(-0.9986432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.36665234) q[2];
sx q[2];
rz(-1.6795936) q[2];
sx q[2];
rz(0.53202638) q[2];
rz(2.4074647) q[3];
sx q[3];
rz(-2.8668154) q[3];
sx q[3];
rz(1.1067357) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4183913) q[0];
sx q[0];
rz(-0.98576236) q[0];
sx q[0];
rz(0.080408737) q[0];
rz(0.53781992) q[1];
sx q[1];
rz(-2.0729013) q[1];
sx q[1];
rz(0.72371662) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079637061) q[0];
sx q[0];
rz(-2.300206) q[0];
sx q[0];
rz(-0.60540149) q[0];
x q[1];
rz(2.2494456) q[2];
sx q[2];
rz(-0.79024678) q[2];
sx q[2];
rz(-2.9228766) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5840877) q[1];
sx q[1];
rz(-0.98230108) q[1];
sx q[1];
rz(-3.1397538) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81775093) q[3];
sx q[3];
rz(-1.9649385) q[3];
sx q[3];
rz(-2.6174291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2878652) q[2];
sx q[2];
rz(-1.1327876) q[2];
sx q[2];
rz(0.90616027) q[2];
rz(-2.7791038) q[3];
sx q[3];
rz(-1.5972842) q[3];
sx q[3];
rz(-3.0947065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2694117) q[0];
sx q[0];
rz(-1.325664) q[0];
sx q[0];
rz(1.0915225) q[0];
rz(-0.12256924) q[1];
sx q[1];
rz(-1.7054319) q[1];
sx q[1];
rz(-1.3465808) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1481762) q[0];
sx q[0];
rz(-1.600334) q[0];
sx q[0];
rz(-1.9269153) q[0];
rz(-pi) q[1];
rz(-2.4819751) q[2];
sx q[2];
rz(-1.3085367) q[2];
sx q[2];
rz(-0.66103092) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1344188) q[1];
sx q[1];
rz(-1.5190304) q[1];
sx q[1];
rz(-1.7758796) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4580549) q[3];
sx q[3];
rz(-2.097297) q[3];
sx q[3];
rz(-2.2248588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71318212) q[2];
sx q[2];
rz(-2.0117663) q[2];
sx q[2];
rz(-0.23207363) q[2];
rz(-1.170916) q[3];
sx q[3];
rz(-2.5210095) q[3];
sx q[3];
rz(-0.28461972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.93037) q[0];
sx q[0];
rz(-2.3486597) q[0];
sx q[0];
rz(1.5027745) q[0];
rz(-2.5495095) q[1];
sx q[1];
rz(-1.9860257) q[1];
sx q[1];
rz(-2.2519462) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9735721) q[0];
sx q[0];
rz(-0.73129994) q[0];
sx q[0];
rz(2.2941053) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1355033) q[2];
sx q[2];
rz(-0.5331299) q[2];
sx q[2];
rz(2.1338303) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.020551) q[1];
sx q[1];
rz(-0.36686037) q[1];
sx q[1];
rz(-2.0936784) q[1];
rz(-pi) q[2];
rz(-0.38460807) q[3];
sx q[3];
rz(-1.9536363) q[3];
sx q[3];
rz(2.2205381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.64638102) q[2];
sx q[2];
rz(-1.4633598) q[2];
sx q[2];
rz(-0.67592534) q[2];
rz(1.4365139) q[3];
sx q[3];
rz(-0.33994514) q[3];
sx q[3];
rz(-2.9743312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2109569) q[0];
sx q[0];
rz(-2.8093331) q[0];
sx q[0];
rz(-0.19530547) q[0];
rz(-0.12061128) q[1];
sx q[1];
rz(-0.21574012) q[1];
sx q[1];
rz(2.9511071) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.932175) q[0];
sx q[0];
rz(-0.10917347) q[0];
sx q[0];
rz(0.63430826) q[0];
rz(0.10730524) q[2];
sx q[2];
rz(-2.4621747) q[2];
sx q[2];
rz(-0.36919644) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0234649) q[1];
sx q[1];
rz(-1.6227229) q[1];
sx q[1];
rz(3.0408188) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5267532) q[3];
sx q[3];
rz(-1.9214464) q[3];
sx q[3];
rz(-2.9235554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6327989) q[2];
sx q[2];
rz(-1.4845279) q[2];
sx q[2];
rz(1.18139) q[2];
rz(1.3440291) q[3];
sx q[3];
rz(-0.7414147) q[3];
sx q[3];
rz(2.369829) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33264273) q[0];
sx q[0];
rz(-1.3430261) q[0];
sx q[0];
rz(2.0020265) q[0];
rz(-2.471916) q[1];
sx q[1];
rz(-2.2256336) q[1];
sx q[1];
rz(-1.12961) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.35992) q[0];
sx q[0];
rz(-1.8501213) q[0];
sx q[0];
rz(-1.2646227) q[0];
x q[1];
rz(-0.30252075) q[2];
sx q[2];
rz(-2.4008958) q[2];
sx q[2];
rz(2.4940235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.46907963) q[1];
sx q[1];
rz(-1.3762849) q[1];
sx q[1];
rz(-2.7698293) q[1];
rz(-pi) q[2];
rz(2.9608742) q[3];
sx q[3];
rz(-1.7263432) q[3];
sx q[3];
rz(-2.8441518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.42664042) q[2];
sx q[2];
rz(-1.4098097) q[2];
sx q[2];
rz(-0.40531522) q[2];
rz(-2.3146368) q[3];
sx q[3];
rz(-0.6260286) q[3];
sx q[3];
rz(-0.85957447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7760794) q[0];
sx q[0];
rz(-3.1190393) q[0];
sx q[0];
rz(2.205701) q[0];
rz(0.86839688) q[1];
sx q[1];
rz(-0.52803841) q[1];
sx q[1];
rz(1.7220928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98631937) q[0];
sx q[0];
rz(-2.5012448) q[0];
sx q[0];
rz(-2.5223283) q[0];
rz(-2.0777736) q[2];
sx q[2];
rz(-1.3666743) q[2];
sx q[2];
rz(-0.43694556) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81379997) q[1];
sx q[1];
rz(-0.84974242) q[1];
sx q[1];
rz(-2.2071597) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9299632) q[3];
sx q[3];
rz(-3.1397925) q[3];
sx q[3];
rz(-1.3789919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3203656) q[2];
sx q[2];
rz(-1.3996404) q[2];
sx q[2];
rz(-2.121675) q[2];
rz(-2.6323281) q[3];
sx q[3];
rz(-1.7002707) q[3];
sx q[3];
rz(-0.078484623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7996063) q[0];
sx q[0];
rz(-1.6260363) q[0];
sx q[0];
rz(1.9827783) q[0];
rz(0.47053567) q[1];
sx q[1];
rz(-1.4152941) q[1];
sx q[1];
rz(1.4656969) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62773317) q[0];
sx q[0];
rz(-0.056170551) q[0];
sx q[0];
rz(0.43561952) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.450804) q[2];
sx q[2];
rz(-1.8422519) q[2];
sx q[2];
rz(1.9631621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50087167) q[1];
sx q[1];
rz(-1.4939918) q[1];
sx q[1];
rz(-2.4810664) q[1];
rz(0.64792525) q[3];
sx q[3];
rz(-1.8019391) q[3];
sx q[3];
rz(-2.1794139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1450242) q[2];
sx q[2];
rz(-1.8137167) q[2];
sx q[2];
rz(2.3222951) q[2];
rz(0.79646349) q[3];
sx q[3];
rz(-2.7345246) q[3];
sx q[3];
rz(-1.6527294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9374989) q[0];
sx q[0];
rz(-3.0452073) q[0];
sx q[0];
rz(-0.82164422) q[0];
rz(-2.4687528) q[1];
sx q[1];
rz(-2.5655589) q[1];
sx q[1];
rz(-1.0427262) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4194787) q[0];
sx q[0];
rz(-1.4983777) q[0];
sx q[0];
rz(-0.2054604) q[0];
rz(-0.94093948) q[2];
sx q[2];
rz(-2.2374638) q[2];
sx q[2];
rz(-0.43064865) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.753841) q[1];
sx q[1];
rz(-2.4341704) q[1];
sx q[1];
rz(0.97960569) q[1];
rz(-pi) q[2];
rz(1.4185227) q[3];
sx q[3];
rz(-1.8209753) q[3];
sx q[3];
rz(0.45246779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3599856) q[2];
sx q[2];
rz(-1.8693417) q[2];
sx q[2];
rz(-0.84235111) q[2];
rz(2.602747) q[3];
sx q[3];
rz(-1.4875965) q[3];
sx q[3];
rz(2.6174788) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6738324) q[0];
sx q[0];
rz(-2.9957275) q[0];
sx q[0];
rz(1.2354596) q[0];
rz(-0.94888672) q[1];
sx q[1];
rz(-0.83273879) q[1];
sx q[1];
rz(1.4804776) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2713) q[0];
sx q[0];
rz(-0.90134951) q[0];
sx q[0];
rz(2.7842194) q[0];
x q[1];
rz(0.83809488) q[2];
sx q[2];
rz(-2.3857255) q[2];
sx q[2];
rz(-0.91516337) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.93714156) q[1];
sx q[1];
rz(-1.6455263) q[1];
sx q[1];
rz(0.16507574) q[1];
rz(-0.62678316) q[3];
sx q[3];
rz(-0.39817444) q[3];
sx q[3];
rz(0.45792031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71198717) q[2];
sx q[2];
rz(-1.4236071) q[2];
sx q[2];
rz(2.831366) q[2];
rz(2.6817536) q[3];
sx q[3];
rz(-0.55791563) q[3];
sx q[3];
rz(-1.3742113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9643758) q[0];
sx q[0];
rz(-1.1342659) q[0];
sx q[0];
rz(-2.076617) q[0];
rz(2.2433157) q[1];
sx q[1];
rz(-0.54294642) q[1];
sx q[1];
rz(-0.44566659) q[1];
rz(-2.0736135) q[2];
sx q[2];
rz(-0.72952727) q[2];
sx q[2];
rz(-2.9222957) q[2];
rz(-0.9089009) q[3];
sx q[3];
rz(-1.2596957) q[3];
sx q[3];
rz(-1.1923292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
