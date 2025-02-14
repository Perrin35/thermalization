OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0187674) q[0];
sx q[0];
rz(-3.1376165) q[0];
sx q[0];
rz(3.0206326) q[0];
rz(-2.6236985) q[1];
sx q[1];
rz(-2.8105812) q[1];
sx q[1];
rz(-1.9537227) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8796381) q[0];
sx q[0];
rz(-1.2722135) q[0];
sx q[0];
rz(1.6613597) q[0];
rz(-0.14622525) q[2];
sx q[2];
rz(-1.0492965) q[2];
sx q[2];
rz(2.5086049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29562274) q[1];
sx q[1];
rz(-0.69087183) q[1];
sx q[1];
rz(-0.9102896) q[1];
rz(-pi) q[2];
rz(-0.066066001) q[3];
sx q[3];
rz(-2.5089536) q[3];
sx q[3];
rz(-0.047601117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0594222) q[2];
sx q[2];
rz(-0.94258451) q[2];
sx q[2];
rz(1.4284632) q[2];
rz(-1.7510471) q[3];
sx q[3];
rz(-2.3640552) q[3];
sx q[3];
rz(-2.9318504) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9222337) q[0];
sx q[0];
rz(-1.5606422) q[0];
sx q[0];
rz(-1.4571762) q[0];
rz(-0.68151418) q[1];
sx q[1];
rz(-0.69029713) q[1];
sx q[1];
rz(-2.1843279) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10892222) q[0];
sx q[0];
rz(-1.3650948) q[0];
sx q[0];
rz(0.47423307) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8419015) q[2];
sx q[2];
rz(-0.77163358) q[2];
sx q[2];
rz(0.51811213) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73814002) q[1];
sx q[1];
rz(-2.7305718) q[1];
sx q[1];
rz(1.6938126) q[1];
rz(-1.5019341) q[3];
sx q[3];
rz(-1.5810229) q[3];
sx q[3];
rz(-0.49730147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7822781) q[2];
sx q[2];
rz(-1.0474397) q[2];
sx q[2];
rz(2.2347343) q[2];
rz(-2.4685229) q[3];
sx q[3];
rz(-0.38452092) q[3];
sx q[3];
rz(-1.2955906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4562255) q[0];
sx q[0];
rz(-2.4816368) q[0];
sx q[0];
rz(-2.5947156) q[0];
rz(1.7710255) q[1];
sx q[1];
rz(-1.9799045) q[1];
sx q[1];
rz(0.10017698) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7179095) q[0];
sx q[0];
rz(-0.29742213) q[0];
sx q[0];
rz(2.0398159) q[0];
x q[1];
rz(1.9630394) q[2];
sx q[2];
rz(-1.0970976) q[2];
sx q[2];
rz(-0.77609962) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.55771577) q[1];
sx q[1];
rz(-0.85266948) q[1];
sx q[1];
rz(-1.8750983) q[1];
rz(2.4002714) q[3];
sx q[3];
rz(-2.5466515) q[3];
sx q[3];
rz(-0.5812318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.748041) q[2];
sx q[2];
rz(-1.2057722) q[2];
sx q[2];
rz(2.098341) q[2];
rz(2.4539963) q[3];
sx q[3];
rz(-0.50191534) q[3];
sx q[3];
rz(2.8858394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8461175) q[0];
sx q[0];
rz(-2.6961374) q[0];
sx q[0];
rz(-1.0365781) q[0];
rz(-1.8702033) q[1];
sx q[1];
rz(-0.82415736) q[1];
sx q[1];
rz(1.2870671) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4117779) q[0];
sx q[0];
rz(-2.3049816) q[0];
sx q[0];
rz(-0.69363026) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9997226) q[2];
sx q[2];
rz(-2.1867036) q[2];
sx q[2];
rz(-2.3162637) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27582622) q[1];
sx q[1];
rz(-2.4600907) q[1];
sx q[1];
rz(1.1033415) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1274476) q[3];
sx q[3];
rz(-1.8487159) q[3];
sx q[3];
rz(1.2025361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.53671545) q[2];
sx q[2];
rz(-2.776919) q[2];
sx q[2];
rz(-3.1216915) q[2];
rz(-0.59602916) q[3];
sx q[3];
rz(-0.28087956) q[3];
sx q[3];
rz(1.9635487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6562281) q[0];
sx q[0];
rz(-1.2475659) q[0];
sx q[0];
rz(-1.9670991) q[0];
rz(-0.30612048) q[1];
sx q[1];
rz(-1.0447634) q[1];
sx q[1];
rz(-1.7721446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4547538) q[0];
sx q[0];
rz(-0.0072221998) q[0];
sx q[0];
rz(-1.104284) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8218064) q[2];
sx q[2];
rz(-2.5282871) q[2];
sx q[2];
rz(-1.5410739) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7983129) q[1];
sx q[1];
rz(-1.7213744) q[1];
sx q[1];
rz(-1.9184789) q[1];
rz(2.6275614) q[3];
sx q[3];
rz(-1.2054878) q[3];
sx q[3];
rz(1.1082054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8451763) q[2];
sx q[2];
rz(-0.74411074) q[2];
sx q[2];
rz(-0.80294341) q[2];
rz(2.0261649) q[3];
sx q[3];
rz(-0.69474703) q[3];
sx q[3];
rz(-1.7723134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95336103) q[0];
sx q[0];
rz(-0.53934923) q[0];
sx q[0];
rz(-0.11235919) q[0];
rz(2.864783) q[1];
sx q[1];
rz(-1.9748297) q[1];
sx q[1];
rz(0.64754957) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43418068) q[0];
sx q[0];
rz(-2.4828504) q[0];
sx q[0];
rz(2.1408268) q[0];
rz(1.1109839) q[2];
sx q[2];
rz(-0.95470482) q[2];
sx q[2];
rz(-1.3715708) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.12880691) q[1];
sx q[1];
rz(-1.4527078) q[1];
sx q[1];
rz(-0.4507555) q[1];
rz(-0.93490113) q[3];
sx q[3];
rz(-0.29014698) q[3];
sx q[3];
rz(-1.4842509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4623146) q[2];
sx q[2];
rz(-0.15028149) q[2];
sx q[2];
rz(1.798299) q[2];
rz(2.6278833) q[3];
sx q[3];
rz(-1.4880344) q[3];
sx q[3];
rz(-2.546379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7996456) q[0];
sx q[0];
rz(-1.8193614) q[0];
sx q[0];
rz(-2.7213726) q[0];
rz(-1.3601903) q[1];
sx q[1];
rz(-1.1667292) q[1];
sx q[1];
rz(0.80418783) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73690945) q[0];
sx q[0];
rz(-1.0361639) q[0];
sx q[0];
rz(0.44628365) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8560195) q[2];
sx q[2];
rz(-2.336506) q[2];
sx q[2];
rz(2.7325163) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.11629352) q[1];
sx q[1];
rz(-2.5466047) q[1];
sx q[1];
rz(-1.8354675) q[1];
x q[2];
rz(-0.98919373) q[3];
sx q[3];
rz(-2.9448754) q[3];
sx q[3];
rz(0.9087874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5899902) q[2];
sx q[2];
rz(-1.0162105) q[2];
sx q[2];
rz(-1.1691947) q[2];
rz(2.9679838) q[3];
sx q[3];
rz(-2.2326525) q[3];
sx q[3];
rz(1.199523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46840295) q[0];
sx q[0];
rz(-1.0262187) q[0];
sx q[0];
rz(-1.2029368) q[0];
rz(-0.24083336) q[1];
sx q[1];
rz(-2.2151561) q[1];
sx q[1];
rz(0.9333207) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4543264) q[0];
sx q[0];
rz(-0.51652241) q[0];
sx q[0];
rz(-0.83965001) q[0];
rz(-pi) q[1];
rz(-0.9983409) q[2];
sx q[2];
rz(-1.3476552) q[2];
sx q[2];
rz(2.9024692) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62040239) q[1];
sx q[1];
rz(-0.4333843) q[1];
sx q[1];
rz(0.49232011) q[1];
rz(-pi) q[2];
rz(-1.0068137) q[3];
sx q[3];
rz(-1.4419334) q[3];
sx q[3];
rz(-3.0430678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11747083) q[2];
sx q[2];
rz(-1.6880219) q[2];
sx q[2];
rz(-0.13548279) q[2];
rz(2.1425715) q[3];
sx q[3];
rz(-0.24805598) q[3];
sx q[3];
rz(2.5854056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67021543) q[0];
sx q[0];
rz(-2.0129634) q[0];
sx q[0];
rz(1.3742597) q[0];
rz(-2.0595835) q[1];
sx q[1];
rz(-0.78649414) q[1];
sx q[1];
rz(1.5793922) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6830272) q[0];
sx q[0];
rz(-2.6011254) q[0];
sx q[0];
rz(-1.5066654) q[0];
rz(-pi) q[1];
rz(-2.5815046) q[2];
sx q[2];
rz(-1.240452) q[2];
sx q[2];
rz(-0.3696839) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1154707) q[1];
sx q[1];
rz(-1.9023832) q[1];
sx q[1];
rz(-2.0381171) q[1];
x q[2];
rz(2.7855603) q[3];
sx q[3];
rz(-2.1904328) q[3];
sx q[3];
rz(2.5621264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9595327) q[2];
sx q[2];
rz(-1.5740732) q[2];
sx q[2];
rz(0.62425295) q[2];
rz(0.67638451) q[3];
sx q[3];
rz(-1.9703777) q[3];
sx q[3];
rz(0.83522767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4312209) q[0];
sx q[0];
rz(-1.2483163) q[0];
sx q[0];
rz(1.7927908) q[0];
rz(0.30385941) q[1];
sx q[1];
rz(-1.6108578) q[1];
sx q[1];
rz(2.2644728) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22280414) q[0];
sx q[0];
rz(-2.3979146) q[0];
sx q[0];
rz(2.9271462) q[0];
rz(-pi) q[1];
rz(-3.0325012) q[2];
sx q[2];
rz(-1.2477836) q[2];
sx q[2];
rz(0.10650466) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7925229) q[1];
sx q[1];
rz(-1.9510837) q[1];
sx q[1];
rz(-2.040634) q[1];
rz(-pi) q[2];
rz(0.56352625) q[3];
sx q[3];
rz(-1.0358264) q[3];
sx q[3];
rz(0.68171147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0861686) q[2];
sx q[2];
rz(-2.7358027) q[2];
sx q[2];
rz(-1.0322734) q[2];
rz(1.0885193) q[3];
sx q[3];
rz(-1.6531331) q[3];
sx q[3];
rz(-0.77435875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1879723) q[0];
sx q[0];
rz(-0.80360501) q[0];
sx q[0];
rz(-3.0927717) q[0];
rz(1.8232518) q[1];
sx q[1];
rz(-1.7346458) q[1];
sx q[1];
rz(0.55191747) q[1];
rz(-1.8404519) q[2];
sx q[2];
rz(-1.1674623) q[2];
sx q[2];
rz(1.8415378) q[2];
rz(2.0069176) q[3];
sx q[3];
rz(-1.5584617) q[3];
sx q[3];
rz(-0.8212318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
