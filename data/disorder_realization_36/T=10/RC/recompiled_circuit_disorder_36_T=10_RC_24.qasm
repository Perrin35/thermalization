OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6498123) q[0];
sx q[0];
rz(-0.28591135) q[0];
sx q[0];
rz(0.51529348) q[0];
rz(-1.7973068) q[1];
sx q[1];
rz(-0.15434115) q[1];
sx q[1];
rz(-0.57758346) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56395036) q[0];
sx q[0];
rz(-0.9194153) q[0];
sx q[0];
rz(1.7786068) q[0];
rz(-pi) q[1];
rz(-2.6063927) q[2];
sx q[2];
rz(-2.0173965) q[2];
sx q[2];
rz(1.5314147) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2540993) q[1];
sx q[1];
rz(-1.3026397) q[1];
sx q[1];
rz(1.9894132) q[1];
rz(-pi) q[2];
rz(-0.76831423) q[3];
sx q[3];
rz(-2.3294805) q[3];
sx q[3];
rz(1.0031568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.704533) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(-2.4543767) q[2];
rz(2.1263188) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17094831) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.2600391) q[0];
rz(-1.0062224) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(-0.84567436) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0802404) q[0];
sx q[0];
rz(-0.70525673) q[0];
sx q[0];
rz(-0.7028701) q[0];
x q[1];
rz(3.0078366) q[2];
sx q[2];
rz(-1.6998561) q[2];
sx q[2];
rz(-1.2226906) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0452803) q[1];
sx q[1];
rz(-1.7573866) q[1];
sx q[1];
rz(0.50218302) q[1];
x q[2];
rz(0.91652292) q[3];
sx q[3];
rz(-0.73361165) q[3];
sx q[3];
rz(-0.77378002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3530897) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(0.4804002) q[2];
rz(-1.7885615) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(1.9539179) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95124328) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(-2.3685266) q[0];
rz(0.13126016) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(2.0551596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7074791) q[0];
sx q[0];
rz(-2.0322324) q[0];
sx q[0];
rz(-3.1397318) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.418872) q[2];
sx q[2];
rz(-0.65290367) q[2];
sx q[2];
rz(2.8002847) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.31049) q[1];
sx q[1];
rz(-1.576014) q[1];
sx q[1];
rz(-0.28880854) q[1];
rz(-pi) q[2];
rz(-1.1826035) q[3];
sx q[3];
rz(-2.7079294) q[3];
sx q[3];
rz(0.38503034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1290258) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(1.3712937) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4521769) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(-3.0932328) q[0];
rz(2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(-1.6960467) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500344) q[0];
sx q[0];
rz(-1.7507179) q[0];
sx q[0];
rz(-0.66288373) q[0];
rz(-pi) q[1];
rz(-2.7518026) q[2];
sx q[2];
rz(-0.23332694) q[2];
sx q[2];
rz(1.1167655) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1452892) q[1];
sx q[1];
rz(-0.22876303) q[1];
sx q[1];
rz(-2.8537675) q[1];
rz(-pi) q[2];
rz(2.0479855) q[3];
sx q[3];
rz(-2.5683937) q[3];
sx q[3];
rz(-1.5447865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.066102862) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(2.1172822) q[2];
rz(-1.5284437) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(-0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4399453) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(-1.8540927) q[0];
rz(0.31907407) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(2.2873926) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4860977) q[0];
sx q[0];
rz(-0.71394074) q[0];
sx q[0];
rz(1.5833202) q[0];
x q[1];
rz(2.1335667) q[2];
sx q[2];
rz(-2.3587583) q[2];
sx q[2];
rz(-2.7340739) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.90480587) q[1];
sx q[1];
rz(-0.85610897) q[1];
sx q[1];
rz(-2.3805815) q[1];
rz(-2.3417926) q[3];
sx q[3];
rz(-0.721867) q[3];
sx q[3];
rz(2.0197899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95191082) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(-0.577315) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(2.234941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0034870738) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(0.072120897) q[0];
rz(-1.1068608) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(-3.0153826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5835411) q[0];
sx q[0];
rz(-1.3582894) q[0];
sx q[0];
rz(-1.6030747) q[0];
rz(-pi) q[1];
rz(-0.77504471) q[2];
sx q[2];
rz(-1.8473052) q[2];
sx q[2];
rz(-1.0724049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0969442) q[1];
sx q[1];
rz(-1.3828039) q[1];
sx q[1];
rz(-2.2915927) q[1];
rz(-pi) q[2];
rz(-1.8993127) q[3];
sx q[3];
rz(-1.7528755) q[3];
sx q[3];
rz(-1.3086705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2386027) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(0.77511707) q[2];
rz(-0.827968) q[3];
sx q[3];
rz(-0.29100806) q[3];
sx q[3];
rz(2.0194139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.5489952) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(3.0431842) q[0];
rz(1.9495643) q[1];
sx q[1];
rz(-1.3339309) q[1];
sx q[1];
rz(0.55955204) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0347621) q[0];
sx q[0];
rz(-2.1666399) q[0];
sx q[0];
rz(-1.4164657) q[0];
rz(-pi) q[1];
rz(0.19182972) q[2];
sx q[2];
rz(-2.8755113) q[2];
sx q[2];
rz(1.1508458) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0911078) q[1];
sx q[1];
rz(-2.9851966) q[1];
sx q[1];
rz(-2.1662103) q[1];
rz(-0.63013245) q[3];
sx q[3];
rz(-1.4117985) q[3];
sx q[3];
rz(2.7255448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.45903912) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(-0.34379488) q[2];
rz(2.5750459) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664117) q[0];
sx q[0];
rz(-1.8090929) q[0];
sx q[0];
rz(2.4108316) q[0];
rz(-0.14239755) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(0.87160814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14411892) q[0];
sx q[0];
rz(-1.561957) q[0];
sx q[0];
rz(-3.0661422) q[0];
rz(-1.7224738) q[2];
sx q[2];
rz(-1.9722087) q[2];
sx q[2];
rz(-0.36827189) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3155047) q[1];
sx q[1];
rz(-1.2102038) q[1];
sx q[1];
rz(0.70198595) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96846795) q[3];
sx q[3];
rz(-1.6864711) q[3];
sx q[3];
rz(1.0761716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.72620755) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(-2.0020206) q[2];
rz(-1.4987882) q[3];
sx q[3];
rz(-2.747624) q[3];
sx q[3];
rz(-2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.9386439) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(1.9198445) q[0];
rz(0.16601673) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(-1.6171914) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6451384) q[0];
sx q[0];
rz(-1.5997412) q[0];
sx q[0];
rz(1.4883947) q[0];
rz(-1.379307) q[2];
sx q[2];
rz(-2.0404437) q[2];
sx q[2];
rz(-2.7188403) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.081454885) q[1];
sx q[1];
rz(-2.7739035) q[1];
sx q[1];
rz(2.0245488) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41096656) q[3];
sx q[3];
rz(-0.68115679) q[3];
sx q[3];
rz(-2.4364803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1200072) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(2.7900556) q[2];
rz(2.0848138) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(-2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.4979424) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(-1.836401) q[0];
rz(2.7611043) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(-0.25340733) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8587592) q[0];
sx q[0];
rz(-1.6738322) q[0];
sx q[0];
rz(-1.9673002) q[0];
rz(2.017574) q[2];
sx q[2];
rz(-0.3728711) q[2];
sx q[2];
rz(2.1602221) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0903783) q[1];
sx q[1];
rz(-1.9131294) q[1];
sx q[1];
rz(-1.8263837) q[1];
rz(0.6110544) q[3];
sx q[3];
rz(-1.2130034) q[3];
sx q[3];
rz(-0.65294453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(0.5029451) q[2];
rz(2.2425966) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1702561) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(-2.4304216) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(-0.74264991) q[2];
sx q[2];
rz(-0.37336083) q[2];
sx q[2];
rz(0.20867418) q[2];
rz(-2.2181702) q[3];
sx q[3];
rz(-0.9265201) q[3];
sx q[3];
rz(-2.0410782) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];