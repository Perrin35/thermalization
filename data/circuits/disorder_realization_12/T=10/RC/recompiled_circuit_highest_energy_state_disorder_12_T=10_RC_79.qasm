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
rz(2.501261) q[0];
sx q[0];
rz(-2.2202272) q[0];
sx q[0];
rz(1.5322354) q[0];
rz(0.20707239) q[1];
sx q[1];
rz(4.3011811) q[1];
sx q[1];
rz(11.094697) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4764413) q[0];
sx q[0];
rz(-1.8855321) q[0];
sx q[0];
rz(-1.2405618) q[0];
rz(0.81469131) q[2];
sx q[2];
rz(-2.1213795) q[2];
sx q[2];
rz(-1.4514635) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0589367) q[1];
sx q[1];
rz(-1.3322222) q[1];
sx q[1];
rz(0.10175565) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26443414) q[3];
sx q[3];
rz(-0.61226058) q[3];
sx q[3];
rz(1.0972638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6472935) q[2];
sx q[2];
rz(-0.99227253) q[2];
sx q[2];
rz(2.5956019) q[2];
rz(-2.2030988) q[3];
sx q[3];
rz(-0.5135082) q[3];
sx q[3];
rz(-2.688496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0168125) q[0];
sx q[0];
rz(-2.0240968) q[0];
sx q[0];
rz(2.5002531) q[0];
rz(-1.9524139) q[1];
sx q[1];
rz(-0.42354217) q[1];
sx q[1];
rz(2.0372527) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3743312) q[0];
sx q[0];
rz(-0.64909726) q[0];
sx q[0];
rz(-2.8968206) q[0];
rz(-0.32764999) q[2];
sx q[2];
rz(-2.0830685) q[2];
sx q[2];
rz(-1.9180026) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1321261) q[1];
sx q[1];
rz(-0.75003549) q[1];
sx q[1];
rz(2.2723531) q[1];
rz(-pi) q[2];
rz(-0.087314815) q[3];
sx q[3];
rz(-1.5305291) q[3];
sx q[3];
rz(-0.95424679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7426593) q[2];
sx q[2];
rz(-2.0519966) q[2];
sx q[2];
rz(-2.147414) q[2];
rz(1.5587156) q[3];
sx q[3];
rz(-0.67928687) q[3];
sx q[3];
rz(-2.5905632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0068552103) q[0];
sx q[0];
rz(-1.0900494) q[0];
sx q[0];
rz(-2.5110974) q[0];
rz(-0.14671239) q[1];
sx q[1];
rz(-1.2598597) q[1];
sx q[1];
rz(0.86348081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3849626) q[0];
sx q[0];
rz(-0.73548404) q[0];
sx q[0];
rz(1.1138659) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21196694) q[2];
sx q[2];
rz(-1.2951998) q[2];
sx q[2];
rz(0.26108643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8161728) q[1];
sx q[1];
rz(-2.6920941) q[1];
sx q[1];
rz(-1.8155431) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8276958) q[3];
sx q[3];
rz(-1.6180929) q[3];
sx q[3];
rz(1.4236049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8063987) q[2];
sx q[2];
rz(-1.6890182) q[2];
sx q[2];
rz(-2.2550968) q[2];
rz(2.7116306) q[3];
sx q[3];
rz(-2.6468247) q[3];
sx q[3];
rz(2.2494242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6524413) q[0];
sx q[0];
rz(-2.8471071) q[0];
sx q[0];
rz(0.44892204) q[0];
rz(2.4849675) q[1];
sx q[1];
rz(-1.4337599) q[1];
sx q[1];
rz(2.8112559) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2097004) q[0];
sx q[0];
rz(-2.7468514) q[0];
sx q[0];
rz(2.1861211) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22808157) q[2];
sx q[2];
rz(-1.9794399) q[2];
sx q[2];
rz(3.0774866) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8790414) q[1];
sx q[1];
rz(-1.8617668) q[1];
sx q[1];
rz(2.0249184) q[1];
rz(-pi) q[2];
rz(1.4497595) q[3];
sx q[3];
rz(-1.5618298) q[3];
sx q[3];
rz(-3.1154273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0535447) q[2];
sx q[2];
rz(-1.6036754) q[2];
sx q[2];
rz(-0.21928445) q[2];
rz(-0.80026475) q[3];
sx q[3];
rz(-2.181668) q[3];
sx q[3];
rz(0.76103359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4735755) q[0];
sx q[0];
rz(-2.0125084) q[0];
sx q[0];
rz(1.4833204) q[0];
rz(0.56770101) q[1];
sx q[1];
rz(-1.2115819) q[1];
sx q[1];
rz(1.6426881) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9358112) q[0];
sx q[0];
rz(-2.09778) q[0];
sx q[0];
rz(-1.4006268) q[0];
rz(-0.11168555) q[2];
sx q[2];
rz(-3.0431261) q[2];
sx q[2];
rz(-2.6278327) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3659542) q[1];
sx q[1];
rz(-1.3729551) q[1];
sx q[1];
rz(-0.61590804) q[1];
x q[2];
rz(1.7622225) q[3];
sx q[3];
rz(-1.6903121) q[3];
sx q[3];
rz(-1.1690948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6296926) q[2];
sx q[2];
rz(-2.3834159) q[2];
sx q[2];
rz(-2.8153815) q[2];
rz(-0.1263667) q[3];
sx q[3];
rz(-0.56505239) q[3];
sx q[3];
rz(0.27157426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7997953) q[0];
sx q[0];
rz(-1.837715) q[0];
sx q[0];
rz(2.6190992) q[0];
rz(0.74784589) q[1];
sx q[1];
rz(-1.6035085) q[1];
sx q[1];
rz(-1.1572256) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8337473) q[0];
sx q[0];
rz(-0.54824775) q[0];
sx q[0];
rz(-2.1566216) q[0];
x q[1];
rz(2.0339478) q[2];
sx q[2];
rz(-1.9310631) q[2];
sx q[2];
rz(2.2778794) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.044837601) q[1];
sx q[1];
rz(-1.4787349) q[1];
sx q[1];
rz(1.7141059) q[1];
rz(-pi) q[2];
rz(0.22208235) q[3];
sx q[3];
rz(-0.33470585) q[3];
sx q[3];
rz(0.10192733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.010765643) q[2];
sx q[2];
rz(-2.6289434) q[2];
sx q[2];
rz(-1.958468) q[2];
rz(-0.042923953) q[3];
sx q[3];
rz(-1.639651) q[3];
sx q[3];
rz(2.8973798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6965028) q[0];
sx q[0];
rz(-2.255991) q[0];
sx q[0];
rz(2.4687299) q[0];
rz(-0.54975763) q[1];
sx q[1];
rz(-1.8472698) q[1];
sx q[1];
rz(-1.6514282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2047574) q[0];
sx q[0];
rz(-2.1278739) q[0];
sx q[0];
rz(1.7812935) q[0];
rz(-0.76810683) q[2];
sx q[2];
rz(-0.66181493) q[2];
sx q[2];
rz(-2.2827374) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0615599) q[1];
sx q[1];
rz(-1.4747689) q[1];
sx q[1];
rz(-0.43635861) q[1];
rz(-pi) q[2];
rz(0.96236046) q[3];
sx q[3];
rz(-0.058789805) q[3];
sx q[3];
rz(0.99942452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.043639) q[2];
sx q[2];
rz(-2.592228) q[2];
sx q[2];
rz(-1.9926386) q[2];
rz(0.39135459) q[3];
sx q[3];
rz(-1.2201744) q[3];
sx q[3];
rz(-0.84827387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5963762) q[0];
sx q[0];
rz(-1.9356118) q[0];
sx q[0];
rz(2.4429831) q[0];
rz(2.2907603) q[1];
sx q[1];
rz(-0.60092503) q[1];
sx q[1];
rz(0.94815475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9382309) q[0];
sx q[0];
rz(-1.6530418) q[0];
sx q[0];
rz(-2.4191678) q[0];
rz(1.9004563) q[2];
sx q[2];
rz(-1.621581) q[2];
sx q[2];
rz(-2.9542758) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5531456) q[1];
sx q[1];
rz(-1.9737893) q[1];
sx q[1];
rz(-2.1827667) q[1];
rz(0.73581477) q[3];
sx q[3];
rz(-2.1664005) q[3];
sx q[3];
rz(-0.99887139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66905388) q[2];
sx q[2];
rz(-1.4488181) q[2];
sx q[2];
rz(1.6834458) q[2];
rz(-1.091188) q[3];
sx q[3];
rz(-2.2444221) q[3];
sx q[3];
rz(0.18479656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75337306) q[0];
sx q[0];
rz(-2.9599157) q[0];
sx q[0];
rz(1.8411807) q[0];
rz(0.45937195) q[1];
sx q[1];
rz(-1.6875024) q[1];
sx q[1];
rz(0.5836817) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3193286) q[0];
sx q[0];
rz(-0.92389014) q[0];
sx q[0];
rz(0.95563332) q[0];
rz(1.5351717) q[2];
sx q[2];
rz(-1.6093264) q[2];
sx q[2];
rz(-0.57254475) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4286268) q[1];
sx q[1];
rz(-1.9216864) q[1];
sx q[1];
rz(1.0246197) q[1];
rz(-pi) q[2];
rz(-0.2512989) q[3];
sx q[3];
rz(-1.150395) q[3];
sx q[3];
rz(-0.29743089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.61648458) q[2];
sx q[2];
rz(-2.8664092) q[2];
sx q[2];
rz(-0.52213651) q[2];
rz(-0.76256847) q[3];
sx q[3];
rz(-1.4985761) q[3];
sx q[3];
rz(-0.34408072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46227118) q[0];
sx q[0];
rz(-1.2079879) q[0];
sx q[0];
rz(0.59360498) q[0];
rz(0.69315928) q[1];
sx q[1];
rz(-2.3489372) q[1];
sx q[1];
rz(-0.75137442) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63735139) q[0];
sx q[0];
rz(-0.48287409) q[0];
sx q[0];
rz(0.58321799) q[0];
rz(-0.85752731) q[2];
sx q[2];
rz(-0.81999841) q[2];
sx q[2];
rz(1.6423026) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.33022949) q[1];
sx q[1];
rz(-0.6682446) q[1];
sx q[1];
rz(3.038637) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9727207) q[3];
sx q[3];
rz(-1.4948625) q[3];
sx q[3];
rz(2.8266265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1410602) q[2];
sx q[2];
rz(-2.4324721) q[2];
sx q[2];
rz(0.023690311) q[2];
rz(1.116811) q[3];
sx q[3];
rz(-1.6938554) q[3];
sx q[3];
rz(1.1344596) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58917191) q[0];
sx q[0];
rz(-1.9235274) q[0];
sx q[0];
rz(-2.0148475) q[0];
rz(1.3202271) q[1];
sx q[1];
rz(-1.8658493) q[1];
sx q[1];
rz(1.0125926) q[1];
rz(2.7350551) q[2];
sx q[2];
rz(-0.58264049) q[2];
sx q[2];
rz(-1.1140136) q[2];
rz(2.9145225) q[3];
sx q[3];
rz(-1.7961226) q[3];
sx q[3];
rz(0.097573438) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
