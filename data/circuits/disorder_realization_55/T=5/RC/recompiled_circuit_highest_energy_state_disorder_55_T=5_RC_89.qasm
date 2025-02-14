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
rz(0.060184181) q[0];
sx q[0];
rz(-2.0826075) q[0];
sx q[0];
rz(-1.0101779) q[0];
rz(2.3330359) q[1];
sx q[1];
rz(-2.8454236) q[1];
sx q[1];
rz(-0.30997601) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5910019) q[0];
sx q[0];
rz(-0.43568107) q[0];
sx q[0];
rz(-2.8753619) q[0];
x q[1];
rz(-1.5484137) q[2];
sx q[2];
rz(-1.3655919) q[2];
sx q[2];
rz(-2.4930645) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7594362) q[1];
sx q[1];
rz(-1.2431989) q[1];
sx q[1];
rz(2.9345678) q[1];
rz(-pi) q[2];
rz(1.3499851) q[3];
sx q[3];
rz(-0.75149262) q[3];
sx q[3];
rz(1.8210851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6628722) q[2];
sx q[2];
rz(-2.8933849) q[2];
sx q[2];
rz(-3.0459246) q[2];
rz(3.0293448) q[3];
sx q[3];
rz(-2.2662558) q[3];
sx q[3];
rz(1.4314502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9005168) q[0];
sx q[0];
rz(-0.88923419) q[0];
sx q[0];
rz(0.61479968) q[0];
rz(-2.2531033) q[1];
sx q[1];
rz(-1.5526086) q[1];
sx q[1];
rz(2.6420171) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0083778) q[0];
sx q[0];
rz(-1.3691804) q[0];
sx q[0];
rz(1.9033543) q[0];
x q[1];
rz(-0.52421661) q[2];
sx q[2];
rz(-1.4337863) q[2];
sx q[2];
rz(1.3772587) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0608637) q[1];
sx q[1];
rz(-1.398096) q[1];
sx q[1];
rz(-2.2279146) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8039475) q[3];
sx q[3];
rz(-2.1269375) q[3];
sx q[3];
rz(0.82586702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2483612) q[2];
sx q[2];
rz(-1.2373368) q[2];
sx q[2];
rz(3.1207747) q[2];
rz(1.9013532) q[3];
sx q[3];
rz(-1.4695243) q[3];
sx q[3];
rz(-1.1851236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6388539) q[0];
sx q[0];
rz(-0.26832142) q[0];
sx q[0];
rz(2.8079206) q[0];
rz(-2.4036713) q[1];
sx q[1];
rz(-1.2260022) q[1];
sx q[1];
rz(-0.12942448) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11331081) q[0];
sx q[0];
rz(-0.89963642) q[0];
sx q[0];
rz(0.20361237) q[0];
x q[1];
rz(-0.67141165) q[2];
sx q[2];
rz(-2.0575626) q[2];
sx q[2];
rz(1.1306131) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8379535) q[1];
sx q[1];
rz(-0.15175125) q[1];
sx q[1];
rz(0.93317731) q[1];
x q[2];
rz(-0.78227104) q[3];
sx q[3];
rz(-1.2201628) q[3];
sx q[3];
rz(2.3334437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1239329) q[2];
sx q[2];
rz(-1.5256226) q[2];
sx q[2];
rz(-1.370149) q[2];
rz(0.74639368) q[3];
sx q[3];
rz(-1.3179702) q[3];
sx q[3];
rz(-0.082775041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.133404) q[0];
sx q[0];
rz(-1.2097825) q[0];
sx q[0];
rz(0.56513894) q[0];
rz(-0.42916974) q[1];
sx q[1];
rz(-2.4950835) q[1];
sx q[1];
rz(0.82829222) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6713326) q[0];
sx q[0];
rz(-1.4909706) q[0];
sx q[0];
rz(0.74175394) q[0];
x q[1];
rz(2.056869) q[2];
sx q[2];
rz(-1.3157433) q[2];
sx q[2];
rz(1.6549095) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8569717) q[1];
sx q[1];
rz(-1.0880252) q[1];
sx q[1];
rz(-1.8716145) q[1];
x q[2];
rz(-1.9902112) q[3];
sx q[3];
rz(-0.55636251) q[3];
sx q[3];
rz(2.0030947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4565178) q[2];
sx q[2];
rz(-1.061331) q[2];
sx q[2];
rz(-1.7100517) q[2];
rz(-0.93959129) q[3];
sx q[3];
rz(-0.51274931) q[3];
sx q[3];
rz(-0.00069869839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.13570304) q[0];
sx q[0];
rz(-2.8764184) q[0];
sx q[0];
rz(1.3294719) q[0];
rz(1.1520518) q[1];
sx q[1];
rz(-2.0538797) q[1];
sx q[1];
rz(0.00024814127) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084293289) q[0];
sx q[0];
rz(-1.7547742) q[0];
sx q[0];
rz(1.8198245) q[0];
rz(-pi) q[1];
rz(-1.7176487) q[2];
sx q[2];
rz(-0.53895437) q[2];
sx q[2];
rz(-2.1951064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2622288) q[1];
sx q[1];
rz(-1.2254224) q[1];
sx q[1];
rz(2.097258) q[1];
rz(-0.080623589) q[3];
sx q[3];
rz(-2.4084457) q[3];
sx q[3];
rz(-2.6997363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1382711) q[2];
sx q[2];
rz(-1.891581) q[2];
sx q[2];
rz(-0.0058343466) q[2];
rz(0.88998574) q[3];
sx q[3];
rz(-0.68166387) q[3];
sx q[3];
rz(-1.1267004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24790813) q[0];
sx q[0];
rz(-1.4434781) q[0];
sx q[0];
rz(1.6538612) q[0];
rz(0.030390175) q[1];
sx q[1];
rz(-2.0149714) q[1];
sx q[1];
rz(-2.1991275) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68719351) q[0];
sx q[0];
rz(-1.8158578) q[0];
sx q[0];
rz(-1.3769727) q[0];
rz(-2.4320388) q[2];
sx q[2];
rz(-2.7559469) q[2];
sx q[2];
rz(-2.8288159) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2561431) q[1];
sx q[1];
rz(-1.1032618) q[1];
sx q[1];
rz(-0.3216775) q[1];
rz(-pi) q[2];
rz(-1.0932572) q[3];
sx q[3];
rz(-0.52553383) q[3];
sx q[3];
rz(-2.9346187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5421062) q[2];
sx q[2];
rz(-2.3607871) q[2];
sx q[2];
rz(1.3055118) q[2];
rz(0.54715884) q[3];
sx q[3];
rz(-1.9360417) q[3];
sx q[3];
rz(-2.3356596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18043537) q[0];
sx q[0];
rz(-1.0338217) q[0];
sx q[0];
rz(-0.10051522) q[0];
rz(-0.48249498) q[1];
sx q[1];
rz(-1.1588187) q[1];
sx q[1];
rz(2.2844792) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094497546) q[0];
sx q[0];
rz(-0.73323868) q[0];
sx q[0];
rz(0.49219699) q[0];
x q[1];
rz(-1.521342) q[2];
sx q[2];
rz(-2.2170545) q[2];
sx q[2];
rz(0.44909278) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7776612) q[1];
sx q[1];
rz(-0.45501935) q[1];
sx q[1];
rz(2.6520686) q[1];
rz(-pi) q[2];
rz(2.8842912) q[3];
sx q[3];
rz(-2.1273779) q[3];
sx q[3];
rz(-1.0986932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40019217) q[2];
sx q[2];
rz(-0.97101784) q[2];
sx q[2];
rz(2.8391489) q[2];
rz(2.1619469) q[3];
sx q[3];
rz(-1.0023578) q[3];
sx q[3];
rz(-1.8625331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.361146) q[0];
sx q[0];
rz(-0.13872153) q[0];
sx q[0];
rz(-3.1106023) q[0];
rz(0.57394761) q[1];
sx q[1];
rz(-1.4731044) q[1];
sx q[1];
rz(-1.2773638) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.325186) q[0];
sx q[0];
rz(-1.2379097) q[0];
sx q[0];
rz(-1.9247967) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1436815) q[2];
sx q[2];
rz(-2.6027205) q[2];
sx q[2];
rz(0.32921916) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7614224) q[1];
sx q[1];
rz(-1.4974277) q[1];
sx q[1];
rz(2.0580893) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3627618) q[3];
sx q[3];
rz(-0.29104189) q[3];
sx q[3];
rz(2.7179949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0370827) q[2];
sx q[2];
rz(-0.13585486) q[2];
sx q[2];
rz(2.6541397) q[2];
rz(0.69665748) q[3];
sx q[3];
rz(-2.2127377) q[3];
sx q[3];
rz(-0.063974403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0697407) q[0];
sx q[0];
rz(-2.6672279) q[0];
sx q[0];
rz(-2.790614) q[0];
rz(2.9674496) q[1];
sx q[1];
rz(-1.6330279) q[1];
sx q[1];
rz(1.4016271) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.642386) q[0];
sx q[0];
rz(-2.8642352) q[0];
sx q[0];
rz(-0.46210285) q[0];
rz(-pi) q[1];
rz(-2.3910037) q[2];
sx q[2];
rz(-1.2341414) q[2];
sx q[2];
rz(-2.2845993) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3199215) q[1];
sx q[1];
rz(-0.18583365) q[1];
sx q[1];
rz(-1.7286848) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97934874) q[3];
sx q[3];
rz(-2.2205847) q[3];
sx q[3];
rz(-1.5929008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1104687) q[2];
sx q[2];
rz(-0.68507552) q[2];
sx q[2];
rz(-2.5049211) q[2];
rz(1.1104256) q[3];
sx q[3];
rz(-1.597155) q[3];
sx q[3];
rz(-2.6231664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7745895) q[0];
sx q[0];
rz(-2.84802) q[0];
sx q[0];
rz(2.0830182) q[0];
rz(0.038657945) q[1];
sx q[1];
rz(-1.5267173) q[1];
sx q[1];
rz(-1.0640594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.80262) q[0];
sx q[0];
rz(-1.5741328) q[0];
sx q[0];
rz(0.0044857684) q[0];
x q[1];
rz(-1.8594785) q[2];
sx q[2];
rz(-1.836986) q[2];
sx q[2];
rz(2.4634354) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.56434332) q[1];
sx q[1];
rz(-1.5303622) q[1];
sx q[1];
rz(3.0811429) q[1];
rz(-1.6056772) q[3];
sx q[3];
rz(-1.1115171) q[3];
sx q[3];
rz(1.1920795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4260063) q[2];
sx q[2];
rz(-2.6435489) q[2];
sx q[2];
rz(3.0675724) q[2];
rz(2.7600539) q[3];
sx q[3];
rz(-1.8266725) q[3];
sx q[3];
rz(-2.7874302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.030180177) q[0];
sx q[0];
rz(-1.8288061) q[0];
sx q[0];
rz(-0.70963138) q[0];
rz(2.5297655) q[1];
sx q[1];
rz(-0.71129967) q[1];
sx q[1];
rz(1.5536972) q[1];
rz(-2.6371018) q[2];
sx q[2];
rz(-1.2332543) q[2];
sx q[2];
rz(2.6674657) q[2];
rz(-2.0430123) q[3];
sx q[3];
rz(-0.87971148) q[3];
sx q[3];
rz(1.9426027) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
