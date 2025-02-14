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
rz(-2.1036086) q[0];
sx q[0];
rz(-1.684364) q[0];
sx q[0];
rz(-1.7888223) q[0];
rz(-5.0985131) q[1];
sx q[1];
rz(2.4689622) q[1];
sx q[1];
rz(11.050635) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7207137) q[0];
sx q[0];
rz(-1.7569245) q[0];
sx q[0];
rz(1.7722919) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47907655) q[2];
sx q[2];
rz(-0.86515364) q[2];
sx q[2];
rz(-2.3912663) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7616854) q[1];
sx q[1];
rz(-2.5406079) q[1];
sx q[1];
rz(2.4883449) q[1];
rz(-pi) q[2];
rz(0.42077371) q[3];
sx q[3];
rz(-1.8340561) q[3];
sx q[3];
rz(0.0076310633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0036014) q[2];
sx q[2];
rz(-0.10819745) q[2];
sx q[2];
rz(0.18599621) q[2];
rz(-3.1229535) q[3];
sx q[3];
rz(-1.2942856) q[3];
sx q[3];
rz(-2.0712461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5828534) q[0];
sx q[0];
rz(-2.5753729) q[0];
sx q[0];
rz(-0.3183611) q[0];
rz(-2.5045577) q[1];
sx q[1];
rz(-2.36167) q[1];
sx q[1];
rz(2.6723518) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1701654) q[0];
sx q[0];
rz(-1.3999456) q[0];
sx q[0];
rz(-0.20194443) q[0];
x q[1];
rz(-1.9465916) q[2];
sx q[2];
rz(-1.6489121) q[2];
sx q[2];
rz(2.5004435) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7982895) q[1];
sx q[1];
rz(-1.251601) q[1];
sx q[1];
rz(-1.3317687) q[1];
rz(0.0057159609) q[3];
sx q[3];
rz(-0.88086956) q[3];
sx q[3];
rz(-1.6315414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8239173) q[2];
sx q[2];
rz(-0.21024148) q[2];
sx q[2];
rz(-2.4483185) q[2];
rz(-1.8215826) q[3];
sx q[3];
rz(-1.3258679) q[3];
sx q[3];
rz(1.0953974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1410809) q[0];
sx q[0];
rz(-2.5040099) q[0];
sx q[0];
rz(-0.47955036) q[0];
rz(-0.50411049) q[1];
sx q[1];
rz(-1.1481552) q[1];
sx q[1];
rz(0.058301059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8970015) q[0];
sx q[0];
rz(-2.069733) q[0];
sx q[0];
rz(1.8656436) q[0];
rz(-2.0928755) q[2];
sx q[2];
rz(-2.1914542) q[2];
sx q[2];
rz(-1.716937) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.95851499) q[1];
sx q[1];
rz(-0.56499186) q[1];
sx q[1];
rz(1.3426258) q[1];
rz(-pi) q[2];
x q[2];
rz(2.635635) q[3];
sx q[3];
rz(-1.803233) q[3];
sx q[3];
rz(-0.94601226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.38640675) q[2];
sx q[2];
rz(-1.9599954) q[2];
sx q[2];
rz(3.1043261) q[2];
rz(-1.5197598) q[3];
sx q[3];
rz(-0.45791364) q[3];
sx q[3];
rz(-2.7601385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.195049) q[0];
sx q[0];
rz(-2.6730972) q[0];
sx q[0];
rz(-3.0750437) q[0];
rz(-1.5244124) q[1];
sx q[1];
rz(-1.8981372) q[1];
sx q[1];
rz(-0.7349416) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3108705) q[0];
sx q[0];
rz(-1.3390216) q[0];
sx q[0];
rz(2.5481497) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3375912) q[2];
sx q[2];
rz(-2.2873023) q[2];
sx q[2];
rz(2.3839714) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1030324) q[1];
sx q[1];
rz(-1.8502697) q[1];
sx q[1];
rz(-0.36906645) q[1];
x q[2];
rz(-2.2913886) q[3];
sx q[3];
rz(-1.6805887) q[3];
sx q[3];
rz(-1.3444855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28896004) q[2];
sx q[2];
rz(-1.5776878) q[2];
sx q[2];
rz(3.0328879) q[2];
rz(2.2147801) q[3];
sx q[3];
rz(-2.1709397) q[3];
sx q[3];
rz(-1.5638117) q[3];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0842593) q[0];
sx q[0];
rz(-2.5045392) q[0];
sx q[0];
rz(-2.0507226) q[0];
rz(-0.57251656) q[1];
sx q[1];
rz(-1.9520091) q[1];
sx q[1];
rz(-2.3172839) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.123468) q[0];
sx q[0];
rz(-1.195692) q[0];
sx q[0];
rz(-0.25183046) q[0];
rz(-pi) q[1];
rz(-0.9277497) q[2];
sx q[2];
rz(-1.130419) q[2];
sx q[2];
rz(-1.7945031) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6047751) q[1];
sx q[1];
rz(-1.0295492) q[1];
sx q[1];
rz(2.6763335) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1371255) q[3];
sx q[3];
rz(-1.1992362) q[3];
sx q[3];
rz(-2.5783758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6606286) q[2];
sx q[2];
rz(-2.169256) q[2];
sx q[2];
rz(0.29120293) q[2];
rz(2.5045942) q[3];
sx q[3];
rz(-1.6773418) q[3];
sx q[3];
rz(-1.8891107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(1.3785192) q[0];
sx q[0];
rz(-2.851649) q[0];
sx q[0];
rz(-3.0105403) q[0];
rz(1.9746926) q[1];
sx q[1];
rz(-2.1578372) q[1];
sx q[1];
rz(0.022620591) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2065434) q[0];
sx q[0];
rz(-0.66842043) q[0];
sx q[0];
rz(0.58983012) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51949595) q[2];
sx q[2];
rz(-0.24458376) q[2];
sx q[2];
rz(1.0102538) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9359259) q[1];
sx q[1];
rz(-2.4956411) q[1];
sx q[1];
rz(-1.0420858) q[1];
rz(-pi) q[2];
rz(0.87225391) q[3];
sx q[3];
rz(-2.8917851) q[3];
sx q[3];
rz(-1.8810617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.268198) q[2];
sx q[2];
rz(-0.70316535) q[2];
sx q[2];
rz(-1.084729) q[2];
rz(1.9966513) q[3];
sx q[3];
rz(-1.8549253) q[3];
sx q[3];
rz(-1.0219319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5893843) q[0];
sx q[0];
rz(-1.5188058) q[0];
sx q[0];
rz(2.6680706) q[0];
rz(1.2208968) q[1];
sx q[1];
rz(-0.41243204) q[1];
sx q[1];
rz(1.268505) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3354123) q[0];
sx q[0];
rz(-1.0632733) q[0];
sx q[0];
rz(0.96854676) q[0];
rz(-pi) q[1];
rz(-2.4051771) q[2];
sx q[2];
rz(-1.0754977) q[2];
sx q[2];
rz(2.8335932) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.699325) q[1];
sx q[1];
rz(-2.0751157) q[1];
sx q[1];
rz(-2.0513351) q[1];
x q[2];
rz(-0.082793391) q[3];
sx q[3];
rz(-1.1361381) q[3];
sx q[3];
rz(2.5744048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1641757) q[2];
sx q[2];
rz(-1.5693393) q[2];
sx q[2];
rz(-2.8951728) q[2];
rz(0.94270802) q[3];
sx q[3];
rz(-2.0556512) q[3];
sx q[3];
rz(-0.76343083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8134269) q[0];
sx q[0];
rz(-1.8712217) q[0];
sx q[0];
rz(-3.0810007) q[0];
rz(-1.5024441) q[1];
sx q[1];
rz(-1.363058) q[1];
sx q[1];
rz(-1.5938119) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13218205) q[0];
sx q[0];
rz(-1.5746467) q[0];
sx q[0];
rz(1.5718223) q[0];
rz(-pi) q[1];
rz(-2.8558989) q[2];
sx q[2];
rz(-2.5473352) q[2];
sx q[2];
rz(2.2069634) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.036435) q[1];
sx q[1];
rz(-2.1686825) q[1];
sx q[1];
rz(2.6041241) q[1];
rz(-pi) q[2];
rz(-0.5129359) q[3];
sx q[3];
rz(-0.78425927) q[3];
sx q[3];
rz(1.605214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26555201) q[2];
sx q[2];
rz(-0.92090845) q[2];
sx q[2];
rz(-1.6843686) q[2];
rz(0.11988104) q[3];
sx q[3];
rz(-0.20457743) q[3];
sx q[3];
rz(2.1129107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2940755) q[0];
sx q[0];
rz(-0.99440044) q[0];
sx q[0];
rz(2.0062398) q[0];
rz(-3.0382233) q[1];
sx q[1];
rz(-1.7366948) q[1];
sx q[1];
rz(-2.1939383) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68610588) q[0];
sx q[0];
rz(-2.4999833) q[0];
sx q[0];
rz(-0.11193402) q[0];
rz(-pi) q[1];
rz(2.728711) q[2];
sx q[2];
rz(-2.2165809) q[2];
sx q[2];
rz(-1.6939163) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3909437) q[1];
sx q[1];
rz(-0.90439561) q[1];
sx q[1];
rz(0.34088582) q[1];
rz(-pi) q[2];
rz(-2.2470948) q[3];
sx q[3];
rz(-0.94493659) q[3];
sx q[3];
rz(1.8602399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.4244298) q[2];
sx q[2];
rz(-0.62817502) q[2];
sx q[2];
rz(-0.90432811) q[2];
rz(-1.2633911) q[3];
sx q[3];
rz(-1.2369913) q[3];
sx q[3];
rz(-0.5164856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.271027) q[0];
sx q[0];
rz(-0.75815433) q[0];
sx q[0];
rz(0.98536056) q[0];
rz(2.789978) q[1];
sx q[1];
rz(-2.1814587) q[1];
sx q[1];
rz(2.7947289) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064929334) q[0];
sx q[0];
rz(-2.4682625) q[0];
sx q[0];
rz(1.9253741) q[0];
x q[1];
rz(-0.11254452) q[2];
sx q[2];
rz(-2.0142609) q[2];
sx q[2];
rz(2.7751768) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4046487) q[1];
sx q[1];
rz(-2.258806) q[1];
sx q[1];
rz(-1.6860123) q[1];
rz(-pi) q[2];
rz(2.4584387) q[3];
sx q[3];
rz(-0.54617631) q[3];
sx q[3];
rz(-0.3953631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2615307) q[2];
sx q[2];
rz(-2.4226649) q[2];
sx q[2];
rz(-3.0066709) q[2];
rz(0.12923446) q[3];
sx q[3];
rz(-2.7415469) q[3];
sx q[3];
rz(1.038704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018910949) q[0];
sx q[0];
rz(-1.8862579) q[0];
sx q[0];
rz(0.12611623) q[0];
rz(2.6799754) q[1];
sx q[1];
rz(-1.7414265) q[1];
sx q[1];
rz(0.19207676) q[1];
rz(1.0151403) q[2];
sx q[2];
rz(-2.3873252) q[2];
sx q[2];
rz(1.5671028) q[2];
rz(-0.25855385) q[3];
sx q[3];
rz(-1.4082485) q[3];
sx q[3];
rz(-0.4859336) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
