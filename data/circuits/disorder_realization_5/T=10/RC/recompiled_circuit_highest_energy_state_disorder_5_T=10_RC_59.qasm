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
rz(-1.9755826) q[0];
sx q[0];
rz(-2.2628885) q[0];
sx q[0];
rz(2.9455844) q[0];
rz(-1.4930383) q[1];
sx q[1];
rz(-0.9382481) q[1];
sx q[1];
rz(0.87747639) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9246793) q[0];
sx q[0];
rz(-1.7857505) q[0];
sx q[0];
rz(-1.1101021) q[0];
rz(-pi) q[1];
rz(-1.2219561) q[2];
sx q[2];
rz(-1.9720805) q[2];
sx q[2];
rz(-0.87043412) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.92447313) q[1];
sx q[1];
rz(-1.8397325) q[1];
sx q[1];
rz(-1.1094138) q[1];
rz(-2.9035197) q[3];
sx q[3];
rz(-2.0364982) q[3];
sx q[3];
rz(0.7326441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8087372) q[2];
sx q[2];
rz(-2.4337807) q[2];
sx q[2];
rz(0.037192496) q[2];
rz(0.91341364) q[3];
sx q[3];
rz(-2.1593058) q[3];
sx q[3];
rz(-1.5090401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7971802) q[0];
sx q[0];
rz(-1.0743112) q[0];
sx q[0];
rz(-2.9194226) q[0];
rz(2.9479058) q[1];
sx q[1];
rz(-1.2682468) q[1];
sx q[1];
rz(0.2598612) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65529275) q[0];
sx q[0];
rz(-1.3278786) q[0];
sx q[0];
rz(1.1134558) q[0];
x q[1];
rz(-1.3919984) q[2];
sx q[2];
rz(-1.6769209) q[2];
sx q[2];
rz(2.6957842) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7227962) q[1];
sx q[1];
rz(-1.3592615) q[1];
sx q[1];
rz(-1.3822894) q[1];
rz(2.6333359) q[3];
sx q[3];
rz(-2.314869) q[3];
sx q[3];
rz(-1.2073769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3746609) q[2];
sx q[2];
rz(-0.58284512) q[2];
sx q[2];
rz(0.15677162) q[2];
rz(-2.3356656) q[3];
sx q[3];
rz(-1.0904049) q[3];
sx q[3];
rz(-2.6167615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2369775) q[0];
sx q[0];
rz(-0.1122864) q[0];
sx q[0];
rz(-1.1861381) q[0];
rz(3.0778432) q[1];
sx q[1];
rz(-1.6155764) q[1];
sx q[1];
rz(-0.41248163) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97226364) q[0];
sx q[0];
rz(-2.0893851) q[0];
sx q[0];
rz(-2.5725468) q[0];
x q[1];
rz(-0.89084741) q[2];
sx q[2];
rz(-1.4483223) q[2];
sx q[2];
rz(0.16799957) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8195678) q[1];
sx q[1];
rz(-1.3680927) q[1];
sx q[1];
rz(2.8611819) q[1];
rz(-0.58150141) q[3];
sx q[3];
rz(-0.99101725) q[3];
sx q[3];
rz(-1.2691154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5731262) q[2];
sx q[2];
rz(-2.7012479) q[2];
sx q[2];
rz(-1.9754515) q[2];
rz(0.79683534) q[3];
sx q[3];
rz(-1.3692057) q[3];
sx q[3];
rz(2.3603175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2260988) q[0];
sx q[0];
rz(-1.508536) q[0];
sx q[0];
rz(1.728212) q[0];
rz(-2.6426897) q[1];
sx q[1];
rz(-1.3830802) q[1];
sx q[1];
rz(-1.2938719) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2436566) q[0];
sx q[0];
rz(-1.5439646) q[0];
sx q[0];
rz(-1.6559589) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0656112) q[2];
sx q[2];
rz(-1.2442028) q[2];
sx q[2];
rz(-3.1103771) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7728987) q[1];
sx q[1];
rz(-0.062090896) q[1];
sx q[1];
rz(0.05120488) q[1];
rz(-0.067469723) q[3];
sx q[3];
rz(-1.6327792) q[3];
sx q[3];
rz(-0.27674473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6125907) q[2];
sx q[2];
rz(-2.4606073) q[2];
sx q[2];
rz(0.16981086) q[2];
rz(-0.42810193) q[3];
sx q[3];
rz(-1.4295108) q[3];
sx q[3];
rz(-2.9812109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.47996461) q[0];
sx q[0];
rz(-2.7492838) q[0];
sx q[0];
rz(-0.83537927) q[0];
rz(-0.63940489) q[1];
sx q[1];
rz(-2.4765922) q[1];
sx q[1];
rz(2.0506052) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24966403) q[0];
sx q[0];
rz(-0.87452379) q[0];
sx q[0];
rz(1.5868756) q[0];
rz(0.29971896) q[2];
sx q[2];
rz(-1.4106361) q[2];
sx q[2];
rz(1.0417787) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.450836) q[1];
sx q[1];
rz(-0.95162205) q[1];
sx q[1];
rz(0.49003933) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8150878) q[3];
sx q[3];
rz(-1.0708263) q[3];
sx q[3];
rz(-2.1799991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7731446) q[2];
sx q[2];
rz(-1.9545363) q[2];
sx q[2];
rz(-0.0040357987) q[2];
rz(-1.076738) q[3];
sx q[3];
rz(-2.4439947) q[3];
sx q[3];
rz(0.18690404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8232329) q[0];
sx q[0];
rz(-1.7802745) q[0];
sx q[0];
rz(2.5079492) q[0];
rz(-0.40114316) q[1];
sx q[1];
rz(-0.70924962) q[1];
sx q[1];
rz(0.69222442) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5790235) q[0];
sx q[0];
rz(-2.2897567) q[0];
sx q[0];
rz(2.1301053) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2622617) q[2];
sx q[2];
rz(-0.82595982) q[2];
sx q[2];
rz(0.334563) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11897381) q[1];
sx q[1];
rz(-1.8198208) q[1];
sx q[1];
rz(-2.157446) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68128392) q[3];
sx q[3];
rz(-2.9046106) q[3];
sx q[3];
rz(-2.9221688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9358518) q[2];
sx q[2];
rz(-0.78448272) q[2];
sx q[2];
rz(1.9007696) q[2];
rz(1.1848909) q[3];
sx q[3];
rz(-1.840206) q[3];
sx q[3];
rz(-1.544781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1320837) q[0];
sx q[0];
rz(-1.3030095) q[0];
sx q[0];
rz(-1.7708154) q[0];
rz(-0.41807434) q[1];
sx q[1];
rz(-1.6590174) q[1];
sx q[1];
rz(-0.48666993) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3330227) q[0];
sx q[0];
rz(-1.6679224) q[0];
sx q[0];
rz(-1.0104695) q[0];
x q[1];
rz(-2.1213221) q[2];
sx q[2];
rz(-1.2472868) q[2];
sx q[2];
rz(-0.24187096) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0907862) q[1];
sx q[1];
rz(-1.3693083) q[1];
sx q[1];
rz(-0.18092107) q[1];
rz(-pi) q[2];
rz(1.6177931) q[3];
sx q[3];
rz(-1.1011657) q[3];
sx q[3];
rz(2.3835973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.439094) q[2];
sx q[2];
rz(-1.0489901) q[2];
sx q[2];
rz(0.19719633) q[2];
rz(2.6321453) q[3];
sx q[3];
rz(-0.26064894) q[3];
sx q[3];
rz(2.313224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8615123) q[0];
sx q[0];
rz(-2.0674288) q[0];
sx q[0];
rz(1.3664838) q[0];
rz(-0.81661433) q[1];
sx q[1];
rz(-2.0520703) q[1];
sx q[1];
rz(-0.28269592) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82474101) q[0];
sx q[0];
rz(-0.62559375) q[0];
sx q[0];
rz(1.2600785) q[0];
x q[1];
rz(-0.28699283) q[2];
sx q[2];
rz(-0.32397917) q[2];
sx q[2];
rz(2.2215077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7273409) q[1];
sx q[1];
rz(-1.7569572) q[1];
sx q[1];
rz(-0.8826137) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62039305) q[3];
sx q[3];
rz(-0.61479688) q[3];
sx q[3];
rz(-2.7063587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1064328) q[2];
sx q[2];
rz(-1.0059493) q[2];
sx q[2];
rz(-2.1584568) q[2];
rz(1.2784917) q[3];
sx q[3];
rz(-2.216414) q[3];
sx q[3];
rz(-2.2127051) q[3];
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
rz(-0.66016692) q[0];
sx q[0];
rz(-1.7919414) q[0];
sx q[0];
rz(-2.8283258) q[0];
rz(0.52472862) q[1];
sx q[1];
rz(-2.8786761) q[1];
sx q[1];
rz(-1.8892586) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1850644) q[0];
sx q[0];
rz(-1.4502429) q[0];
sx q[0];
rz(0.53114364) q[0];
x q[1];
rz(-2.1519442) q[2];
sx q[2];
rz(-2.2454442) q[2];
sx q[2];
rz(2.0591813) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.508751) q[1];
sx q[1];
rz(-1.0392042) q[1];
sx q[1];
rz(0.94137194) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.081324025) q[3];
sx q[3];
rz(-0.81800753) q[3];
sx q[3];
rz(0.52018702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.19920443) q[2];
sx q[2];
rz(-0.62493268) q[2];
sx q[2];
rz(1.3814629) q[2];
rz(2.840461) q[3];
sx q[3];
rz(-0.98265177) q[3];
sx q[3];
rz(0.38100955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23406601) q[0];
sx q[0];
rz(-1.0028239) q[0];
sx q[0];
rz(-1.7940849) q[0];
rz(-0.8849591) q[1];
sx q[1];
rz(-0.62565175) q[1];
sx q[1];
rz(1.7431097) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8475474) q[0];
sx q[0];
rz(-1.9452403) q[0];
sx q[0];
rz(-0.011463077) q[0];
rz(-pi) q[1];
rz(-2.276097) q[2];
sx q[2];
rz(-1.4303606) q[2];
sx q[2];
rz(0.27117929) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3264638) q[1];
sx q[1];
rz(-0.088610709) q[1];
sx q[1];
rz(-0.43472241) q[1];
rz(-pi) q[2];
rz(1.2169525) q[3];
sx q[3];
rz(-2.5501336) q[3];
sx q[3];
rz(-0.22218695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.053190319) q[2];
sx q[2];
rz(-1.3206427) q[2];
sx q[2];
rz(-0.94584805) q[2];
rz(0.69019067) q[3];
sx q[3];
rz(-1.2226356) q[3];
sx q[3];
rz(-1.8405731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90203862) q[0];
sx q[0];
rz(-1.6033462) q[0];
sx q[0];
rz(-0.70815804) q[0];
rz(3.0802849) q[1];
sx q[1];
rz(-2.5262482) q[1];
sx q[1];
rz(1.6940438) q[1];
rz(-1.4881143) q[2];
sx q[2];
rz(-1.8607685) q[2];
sx q[2];
rz(0.49989732) q[2];
rz(0.072515247) q[3];
sx q[3];
rz(-1.7504277) q[3];
sx q[3];
rz(-0.75274368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
