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
rz(0.68501002) q[0];
sx q[0];
rz(3.8334414) q[0];
sx q[0];
rz(9.8631996) q[0];
rz(-0.11998478) q[1];
sx q[1];
rz(-0.24155231) q[1];
sx q[1];
rz(-2.144699) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55132574) q[0];
sx q[0];
rz(-2.9661313) q[0];
sx q[0];
rz(-2.0911123) q[0];
rz(2.8710932) q[2];
sx q[2];
rz(-1.5322313) q[2];
sx q[2];
rz(2.6730516) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5146487) q[1];
sx q[1];
rz(-0.46851003) q[1];
sx q[1];
rz(-1.0131939) q[1];
rz(1.3658872) q[3];
sx q[3];
rz(-1.3384322) q[3];
sx q[3];
rz(0.13162498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3991656) q[2];
sx q[2];
rz(-0.70968598) q[2];
sx q[2];
rz(-0.7134552) q[2];
rz(-2.3136638) q[3];
sx q[3];
rz(-1.6498339) q[3];
sx q[3];
rz(2.2152065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4543318) q[0];
sx q[0];
rz(-2.5255272) q[0];
sx q[0];
rz(-1.8808421) q[0];
rz(2.8997391) q[1];
sx q[1];
rz(-1.2956023) q[1];
sx q[1];
rz(2.5557925) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0457092) q[0];
sx q[0];
rz(-1.2564141) q[0];
sx q[0];
rz(0.16736253) q[0];
rz(-pi) q[1];
rz(-0.42795534) q[2];
sx q[2];
rz(-1.578786) q[2];
sx q[2];
rz(-0.34763476) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3812848) q[1];
sx q[1];
rz(-0.80712748) q[1];
sx q[1];
rz(-1.155608) q[1];
rz(-pi) q[2];
rz(1.582566) q[3];
sx q[3];
rz(-1.571072) q[3];
sx q[3];
rz(0.89645146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5213617) q[2];
sx q[2];
rz(-0.62180454) q[2];
sx q[2];
rz(2.8561031) q[2];
rz(-0.96389687) q[3];
sx q[3];
rz(-0.6670835) q[3];
sx q[3];
rz(-2.9662761) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6730839) q[0];
sx q[0];
rz(-0.54247576) q[0];
sx q[0];
rz(-0.28955224) q[0];
rz(-2.8383004) q[1];
sx q[1];
rz(-1.8692317) q[1];
sx q[1];
rz(1.2264651) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68115409) q[0];
sx q[0];
rz(-1.9931721) q[0];
sx q[0];
rz(0.85737164) q[0];
x q[1];
rz(2.5614024) q[2];
sx q[2];
rz(-1.7805432) q[2];
sx q[2];
rz(-2.8856087) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81193594) q[1];
sx q[1];
rz(-1.2878988) q[1];
sx q[1];
rz(0.97269989) q[1];
rz(-pi) q[2];
rz(1.6885593) q[3];
sx q[3];
rz(-2.3533245) q[3];
sx q[3];
rz(0.39358172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.44846416) q[2];
sx q[2];
rz(-2.144564) q[2];
sx q[2];
rz(-2.2779951) q[2];
rz(-3.038285) q[3];
sx q[3];
rz(-1.7938675) q[3];
sx q[3];
rz(-3.0062655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.1358262) q[0];
sx q[0];
rz(-2.5209881) q[0];
sx q[0];
rz(-3.0964858) q[0];
rz(-0.82334423) q[1];
sx q[1];
rz(-2.7594559) q[1];
sx q[1];
rz(-2.8270922) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.452144) q[0];
sx q[0];
rz(-2.2270791) q[0];
sx q[0];
rz(1.6771862) q[0];
rz(-pi) q[1];
rz(2.1067736) q[2];
sx q[2];
rz(-1.2064486) q[2];
sx q[2];
rz(-1.8672158) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38344774) q[1];
sx q[1];
rz(-0.52844344) q[1];
sx q[1];
rz(-2.7695038) q[1];
x q[2];
rz(-2.567611) q[3];
sx q[3];
rz(-0.38061695) q[3];
sx q[3];
rz(1.5215645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0752252) q[2];
sx q[2];
rz(-1.6547357) q[2];
sx q[2];
rz(0.60869795) q[2];
rz(1.4501976) q[3];
sx q[3];
rz(-0.83289731) q[3];
sx q[3];
rz(2.6628185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4920014) q[0];
sx q[0];
rz(-0.86968017) q[0];
sx q[0];
rz(-3.0539404) q[0];
rz(-1.8805257) q[1];
sx q[1];
rz(-1.7365716) q[1];
sx q[1];
rz(2.2671949) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5090655) q[0];
sx q[0];
rz(-1.6035) q[0];
sx q[0];
rz(-2.2421809) q[0];
rz(-pi) q[1];
rz(1.6220785) q[2];
sx q[2];
rz(-0.59133321) q[2];
sx q[2];
rz(1.4695449) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.27694662) q[1];
sx q[1];
rz(-1.0082933) q[1];
sx q[1];
rz(0.87361305) q[1];
rz(0.22190549) q[3];
sx q[3];
rz(-2.1067224) q[3];
sx q[3];
rz(1.3886896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38292357) q[2];
sx q[2];
rz(-2.1047968) q[2];
sx q[2];
rz(-2.5385638) q[2];
rz(2.1488819) q[3];
sx q[3];
rz(-0.38645667) q[3];
sx q[3];
rz(-2.4934798) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25931609) q[0];
sx q[0];
rz(-0.74463212) q[0];
sx q[0];
rz(-0.69389206) q[0];
rz(-2.4248185) q[1];
sx q[1];
rz(-2.0718772) q[1];
sx q[1];
rz(-0.96376354) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19417441) q[0];
sx q[0];
rz(-0.45233417) q[0];
sx q[0];
rz(0.27249713) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18220696) q[2];
sx q[2];
rz(-1.4002396) q[2];
sx q[2];
rz(-0.86917669) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9020667) q[1];
sx q[1];
rz(-2.0888302) q[1];
sx q[1];
rz(2.7076028) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5429872) q[3];
sx q[3];
rz(-0.99861342) q[3];
sx q[3];
rz(-1.5861301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11513772) q[2];
sx q[2];
rz(-2.852735) q[2];
sx q[2];
rz(3.0333983) q[2];
rz(-2.1859956) q[3];
sx q[3];
rz(-3.1028265) q[3];
sx q[3];
rz(-0.55967104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74681246) q[0];
sx q[0];
rz(-2.96947) q[0];
sx q[0];
rz(-0.10699233) q[0];
rz(-0.50210285) q[1];
sx q[1];
rz(-1.5109477) q[1];
sx q[1];
rz(-2.9923901) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7764588) q[0];
sx q[0];
rz(-1.7452876) q[0];
sx q[0];
rz(-0.15326881) q[0];
x q[1];
rz(-0.72609857) q[2];
sx q[2];
rz(-2.4319785) q[2];
sx q[2];
rz(-2.6626056) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8684959) q[1];
sx q[1];
rz(-2.3567216) q[1];
sx q[1];
rz(-2.0071061) q[1];
rz(1.7140237) q[3];
sx q[3];
rz(-1.9960072) q[3];
sx q[3];
rz(0.79637209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8788098) q[2];
sx q[2];
rz(-2.171319) q[2];
sx q[2];
rz(0.61857569) q[2];
rz(0.68459073) q[3];
sx q[3];
rz(-0.22817831) q[3];
sx q[3];
rz(-0.98606199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743643) q[0];
sx q[0];
rz(-1.3716797) q[0];
sx q[0];
rz(-0.57269639) q[0];
rz(1.0193846) q[1];
sx q[1];
rz(-2.9044594) q[1];
sx q[1];
rz(-3.0651029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6394326) q[0];
sx q[0];
rz(-1.0442088) q[0];
sx q[0];
rz(2.9798085) q[0];
x q[1];
rz(-1.7911081) q[2];
sx q[2];
rz(-1.2800084) q[2];
sx q[2];
rz(-2.5875768) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1566211) q[1];
sx q[1];
rz(-2.467359) q[1];
sx q[1];
rz(-1.7471501) q[1];
x q[2];
rz(-0.88249607) q[3];
sx q[3];
rz(-1.232649) q[3];
sx q[3];
rz(-2.8201617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1028334) q[2];
sx q[2];
rz(-0.88492727) q[2];
sx q[2];
rz(0.49212512) q[2];
rz(-2.7627908) q[3];
sx q[3];
rz(-2.7087961) q[3];
sx q[3];
rz(-2.2545599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4019796) q[0];
sx q[0];
rz(-2.6519863) q[0];
sx q[0];
rz(-0.42994764) q[0];
rz(-2.6509189) q[1];
sx q[1];
rz(-2.6538167) q[1];
sx q[1];
rz(-1.0027764) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96129629) q[0];
sx q[0];
rz(-1.8334098) q[0];
sx q[0];
rz(1.0680593) q[0];
x q[1];
rz(-2.241469) q[2];
sx q[2];
rz(-1.6311833) q[2];
sx q[2];
rz(-0.65298572) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7007737) q[1];
sx q[1];
rz(-2.6731468) q[1];
sx q[1];
rz(-3.0354795) q[1];
rz(-pi) q[2];
rz(-2.0315038) q[3];
sx q[3];
rz(-0.92383251) q[3];
sx q[3];
rz(1.6443192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7542725) q[2];
sx q[2];
rz(-2.5405799) q[2];
sx q[2];
rz(0.69166541) q[2];
rz(0.12868853) q[3];
sx q[3];
rz(-1.5494989) q[3];
sx q[3];
rz(-2.5531829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.9733031) q[0];
sx q[0];
rz(-3.0423218) q[0];
sx q[0];
rz(-0.6231935) q[0];
rz(0.19459952) q[1];
sx q[1];
rz(-2.01229) q[1];
sx q[1];
rz(2.2653939) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7565544) q[0];
sx q[0];
rz(-1.4557739) q[0];
sx q[0];
rz(0.043941078) q[0];
rz(-pi) q[1];
rz(-1.7224465) q[2];
sx q[2];
rz(-0.83497161) q[2];
sx q[2];
rz(0.73260546) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9651803) q[1];
sx q[1];
rz(-1.1438055) q[1];
sx q[1];
rz(1.7199442) q[1];
rz(0.51214062) q[3];
sx q[3];
rz(-2.509553) q[3];
sx q[3];
rz(-1.8919945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1759922) q[2];
sx q[2];
rz(-2.8922562) q[2];
sx q[2];
rz(-2.5052137) q[2];
rz(-1.5198358) q[3];
sx q[3];
rz(-2.2607925) q[3];
sx q[3];
rz(-2.9145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31430055) q[0];
sx q[0];
rz(-1.5416332) q[0];
sx q[0];
rz(2.736349) q[0];
rz(1.4018519) q[1];
sx q[1];
rz(-1.4973462) q[1];
sx q[1];
rz(-3.0037465) q[1];
rz(-2.9707303) q[2];
sx q[2];
rz(-2.5124585) q[2];
sx q[2];
rz(-0.72825904) q[2];
rz(0.059611353) q[3];
sx q[3];
rz(-2.6047299) q[3];
sx q[3];
rz(-2.4951618) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
