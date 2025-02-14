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
rz(1.4577515) q[0];
sx q[0];
rz(-3.0384851) q[0];
sx q[0];
rz(-0.74080324) q[0];
rz(-0.15751547) q[1];
sx q[1];
rz(-2.4060251) q[1];
sx q[1];
rz(-0.29677376) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65771253) q[0];
sx q[0];
rz(-1.1917416) q[0];
sx q[0];
rz(-2.0188278) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27539092) q[2];
sx q[2];
rz(-2.4045334) q[2];
sx q[2];
rz(-0.136497) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9903914) q[1];
sx q[1];
rz(-1.5784653) q[1];
sx q[1];
rz(2.5911314) q[1];
rz(0.17916174) q[3];
sx q[3];
rz(-1.1350313) q[3];
sx q[3];
rz(0.023438862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.890581) q[2];
sx q[2];
rz(-3.0457532) q[2];
sx q[2];
rz(1.7354234) q[2];
rz(0.74875325) q[3];
sx q[3];
rz(-0.65776062) q[3];
sx q[3];
rz(2.9296618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1521848) q[0];
sx q[0];
rz(-2.2693372) q[0];
sx q[0];
rz(2.8023791) q[0];
rz(0.57505125) q[1];
sx q[1];
rz(-1.9564068) q[1];
sx q[1];
rz(2.7599879) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1763827) q[0];
sx q[0];
rz(-1.8032224) q[0];
sx q[0];
rz(1.1430955) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8782042) q[2];
sx q[2];
rz(-1.8991422) q[2];
sx q[2];
rz(1.245861) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2927865) q[1];
sx q[1];
rz(-1.5611575) q[1];
sx q[1];
rz(1.9697492) q[1];
rz(-pi) q[2];
rz(-0.53839243) q[3];
sx q[3];
rz(-0.88148553) q[3];
sx q[3];
rz(1.6834843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7316458) q[2];
sx q[2];
rz(-2.4531328) q[2];
sx q[2];
rz(0.48639578) q[2];
rz(2.5977123) q[3];
sx q[3];
rz(-2.0982274) q[3];
sx q[3];
rz(-0.16415183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.630702) q[0];
sx q[0];
rz(-0.4158026) q[0];
sx q[0];
rz(-2.1422332) q[0];
rz(-2.4103145) q[1];
sx q[1];
rz(-2.6731698) q[1];
sx q[1];
rz(0.17459248) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0492741) q[0];
sx q[0];
rz(-1.1241372) q[0];
sx q[0];
rz(-3.1389278) q[0];
rz(-pi) q[1];
rz(-1.393494) q[2];
sx q[2];
rz(-1.4127222) q[2];
sx q[2];
rz(-0.47249913) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0583971) q[1];
sx q[1];
rz(-0.10279142) q[1];
sx q[1];
rz(1.0344857) q[1];
rz(0.2184775) q[3];
sx q[3];
rz(-1.1214773) q[3];
sx q[3];
rz(1.9868509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.752855) q[2];
sx q[2];
rz(-0.85192215) q[2];
sx q[2];
rz(1.277415) q[2];
rz(-2.3169005) q[3];
sx q[3];
rz(-0.84452355) q[3];
sx q[3];
rz(-0.1674913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.326062) q[0];
sx q[0];
rz(-2.675246) q[0];
sx q[0];
rz(-1.9933568) q[0];
rz(-2.8094021) q[1];
sx q[1];
rz(-0.59993184) q[1];
sx q[1];
rz(-1.0848328) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058131071) q[0];
sx q[0];
rz(-1.8634184) q[0];
sx q[0];
rz(2.9954442) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81697322) q[2];
sx q[2];
rz(-1.3187871) q[2];
sx q[2];
rz(-1.1622045) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12382896) q[1];
sx q[1];
rz(-0.74018407) q[1];
sx q[1];
rz(0.077880903) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6779019) q[3];
sx q[3];
rz(-1.1775908) q[3];
sx q[3];
rz(2.5785407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.79160249) q[2];
sx q[2];
rz(-0.76452667) q[2];
sx q[2];
rz(0.0066268607) q[2];
rz(2.918112) q[3];
sx q[3];
rz(-1.0379182) q[3];
sx q[3];
rz(-2.3124783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70581907) q[0];
sx q[0];
rz(-0.75278246) q[0];
sx q[0];
rz(1.1821795) q[0];
rz(-1.301282) q[1];
sx q[1];
rz(-0.92295206) q[1];
sx q[1];
rz(2.9822947) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1263016) q[0];
sx q[0];
rz(-1.8582004) q[0];
sx q[0];
rz(-1.6096344) q[0];
x q[1];
rz(-1.0403522) q[2];
sx q[2];
rz(-0.36200501) q[2];
sx q[2];
rz(1.1623882) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0634126) q[1];
sx q[1];
rz(-1.0920807) q[1];
sx q[1];
rz(-0.53005752) q[1];
rz(-0.90386541) q[3];
sx q[3];
rz(-1.9124036) q[3];
sx q[3];
rz(1.2056392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84676877) q[2];
sx q[2];
rz(-2.939665) q[2];
sx q[2];
rz(1.798604) q[2];
rz(-0.35074562) q[3];
sx q[3];
rz(-2.0028159) q[3];
sx q[3];
rz(-0.32575592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89966929) q[0];
sx q[0];
rz(-2.3133008) q[0];
sx q[0];
rz(-2.9747466) q[0];
rz(2.9150561) q[1];
sx q[1];
rz(-1.3275361) q[1];
sx q[1];
rz(-2.66364) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8628381) q[0];
sx q[0];
rz(-2.4669381) q[0];
sx q[0];
rz(-2.1341896) q[0];
x q[1];
rz(0.40845334) q[2];
sx q[2];
rz(-2.3801) q[2];
sx q[2];
rz(-3.048427) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8357827) q[1];
sx q[1];
rz(-2.1205582) q[1];
sx q[1];
rz(2.1533986) q[1];
rz(-pi) q[2];
rz(-2.0178363) q[3];
sx q[3];
rz(-1.9634878) q[3];
sx q[3];
rz(-2.6156619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6106674) q[2];
sx q[2];
rz(-1.764856) q[2];
sx q[2];
rz(-1.2923856) q[2];
rz(2.2409706) q[3];
sx q[3];
rz(-0.77018046) q[3];
sx q[3];
rz(-1.5952516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47386277) q[0];
sx q[0];
rz(-2.168499) q[0];
sx q[0];
rz(1.6170172) q[0];
rz(-1.698311) q[1];
sx q[1];
rz(-2.2459005) q[1];
sx q[1];
rz(0.5009833) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9568142) q[0];
sx q[0];
rz(-1.2010472) q[0];
sx q[0];
rz(2.2201204) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47246859) q[2];
sx q[2];
rz(-2.5295265) q[2];
sx q[2];
rz(-1.2304359) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0308548) q[1];
sx q[1];
rz(-0.99910611) q[1];
sx q[1];
rz(-0.60380377) q[1];
rz(-pi) q[2];
rz(-2.8328308) q[3];
sx q[3];
rz(-1.4523066) q[3];
sx q[3];
rz(-2.091696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3524126) q[2];
sx q[2];
rz(-0.12426201) q[2];
sx q[2];
rz(-0.46335709) q[2];
rz(0.067666791) q[3];
sx q[3];
rz(-1.8384408) q[3];
sx q[3];
rz(0.1629924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25245923) q[0];
sx q[0];
rz(-1.2754138) q[0];
sx q[0];
rz(0.5109936) q[0];
rz(2.4976318) q[1];
sx q[1];
rz(-1.124758) q[1];
sx q[1];
rz(2.8535829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12342425) q[0];
sx q[0];
rz(-0.19561681) q[0];
sx q[0];
rz(0.19851144) q[0];
rz(-pi) q[1];
rz(-2.2637469) q[2];
sx q[2];
rz(-1.0509509) q[2];
sx q[2];
rz(2.9016888) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0916413) q[1];
sx q[1];
rz(-2.5041951) q[1];
sx q[1];
rz(-2.6693488) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37841017) q[3];
sx q[3];
rz(-1.6749951) q[3];
sx q[3];
rz(-0.7507594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94706941) q[2];
sx q[2];
rz(-1.1627407) q[2];
sx q[2];
rz(0.21491773) q[2];
rz(-1.8042709) q[3];
sx q[3];
rz(-0.53842068) q[3];
sx q[3];
rz(2.9609093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8682078) q[0];
sx q[0];
rz(-2.2053563) q[0];
sx q[0];
rz(-2.0599763) q[0];
rz(0.99010211) q[1];
sx q[1];
rz(-1.5847881) q[1];
sx q[1];
rz(-0.33499151) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43083999) q[0];
sx q[0];
rz(-0.3246322) q[0];
sx q[0];
rz(1.6765094) q[0];
rz(-0.0063517687) q[2];
sx q[2];
rz(-1.1489023) q[2];
sx q[2];
rz(-0.61074257) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44030467) q[1];
sx q[1];
rz(-3.0106643) q[1];
sx q[1];
rz(-2.4314315) q[1];
x q[2];
rz(-1.5509299) q[3];
sx q[3];
rz(-0.40106138) q[3];
sx q[3];
rz(3.0326518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0830903) q[2];
sx q[2];
rz(-0.52512705) q[2];
sx q[2];
rz(0.38880175) q[2];
rz(0.36924103) q[3];
sx q[3];
rz(-2.8856314) q[3];
sx q[3];
rz(0.84454876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9487172) q[0];
sx q[0];
rz(-0.95272869) q[0];
sx q[0];
rz(2.4396851) q[0];
rz(-1.6096055) q[1];
sx q[1];
rz(-2.46789) q[1];
sx q[1];
rz(-2.7598377) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3385388) q[0];
sx q[0];
rz(-3.038254) q[0];
sx q[0];
rz(-2.4843012) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4456621) q[2];
sx q[2];
rz(-1.9788673) q[2];
sx q[2];
rz(-3.0681075) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.020479105) q[1];
sx q[1];
rz(-1.319863) q[1];
sx q[1];
rz(-2.8931055) q[1];
rz(-pi) q[2];
rz(2.3424861) q[3];
sx q[3];
rz(-1.0630084) q[3];
sx q[3];
rz(0.72572177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4664885) q[2];
sx q[2];
rz(-2.1110822) q[2];
sx q[2];
rz(2.4934736) q[2];
rz(-1.0068007) q[3];
sx q[3];
rz(-1.9961793) q[3];
sx q[3];
rz(-3.0692611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61160144) q[0];
sx q[0];
rz(-1.4016822) q[0];
sx q[0];
rz(0.66521426) q[0];
rz(-0.29162677) q[1];
sx q[1];
rz(-1.7886152) q[1];
sx q[1];
rz(2.0033966) q[1];
rz(-2.6219528) q[2];
sx q[2];
rz(-1.8815132) q[2];
sx q[2];
rz(2.5757679) q[2];
rz(2.9838647) q[3];
sx q[3];
rz(-1.986151) q[3];
sx q[3];
rz(2.3421866) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
