OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55750027) q[0];
sx q[0];
rz(-3.1182365) q[0];
sx q[0];
rz(-0.93552247) q[0];
rz(-1.6159396) q[1];
sx q[1];
rz(-1.5204117) q[1];
sx q[1];
rz(2.867155) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4532758) q[0];
sx q[0];
rz(-0.098628086) q[0];
sx q[0];
rz(-2.3345956) q[0];
x q[1];
rz(-2.6018167) q[2];
sx q[2];
rz(-1.1640884) q[2];
sx q[2];
rz(-1.1193898) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1573616) q[1];
sx q[1];
rz(-1.5604291) q[1];
sx q[1];
rz(1.5342516) q[1];
x q[2];
rz(-3.0679061) q[3];
sx q[3];
rz(-1.7409803) q[3];
sx q[3];
rz(-0.88995349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9258257) q[2];
sx q[2];
rz(-0.01130686) q[2];
sx q[2];
rz(-1.0781778) q[2];
rz(-0.82873851) q[3];
sx q[3];
rz(-1.5107892) q[3];
sx q[3];
rz(0.80882788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0593798) q[0];
sx q[0];
rz(-1.8635211) q[0];
sx q[0];
rz(1.3905806) q[0];
rz(-1.7104205) q[1];
sx q[1];
rz(-0.0043967604) q[1];
sx q[1];
rz(1.4336525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9427458) q[0];
sx q[0];
rz(-2.2042243) q[0];
sx q[0];
rz(3.01157) q[0];
x q[1];
rz(2.0145871) q[2];
sx q[2];
rz(-1.5559042) q[2];
sx q[2];
rz(-3.1010951) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2428618) q[1];
sx q[1];
rz(-3.1325097) q[1];
sx q[1];
rz(-1.5255949) q[1];
x q[2];
rz(-1.6363899) q[3];
sx q[3];
rz(-0.7842614) q[3];
sx q[3];
rz(2.0900871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37890515) q[2];
sx q[2];
rz(-1.5964369) q[2];
sx q[2];
rz(1.572466) q[2];
rz(0.76111859) q[3];
sx q[3];
rz(-3.0877536) q[3];
sx q[3];
rz(-1.9332473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.222027) q[0];
sx q[0];
rz(-0.61028218) q[0];
sx q[0];
rz(-0.56184226) q[0];
rz(-1.5737083) q[1];
sx q[1];
rz(-1.578873) q[1];
sx q[1];
rz(3.124253) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2527232) q[0];
sx q[0];
rz(-2.4742352) q[0];
sx q[0];
rz(0.86020893) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2566913) q[2];
sx q[2];
rz(-1.1807311) q[2];
sx q[2];
rz(-1.6194413) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8120809) q[1];
sx q[1];
rz(-1.9332262) q[1];
sx q[1];
rz(-0.017048841) q[1];
rz(1.5937599) q[3];
sx q[3];
rz(-1.0606137) q[3];
sx q[3];
rz(1.9397473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5138381) q[2];
sx q[2];
rz(-1.7226115) q[2];
sx q[2];
rz(2.5886152) q[2];
rz(-1.9364457) q[3];
sx q[3];
rz(-1.5744934) q[3];
sx q[3];
rz(1.5945826) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39128458) q[0];
sx q[0];
rz(-1.0306083) q[0];
sx q[0];
rz(-1.3543825) q[0];
rz(-1.3722108) q[1];
sx q[1];
rz(-0.0028227614) q[1];
sx q[1];
rz(1.7599958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8244265) q[0];
sx q[0];
rz(-1.1689416) q[0];
sx q[0];
rz(2.0749932) q[0];
x q[1];
rz(3.1388248) q[2];
sx q[2];
rz(-1.5684897) q[2];
sx q[2];
rz(-0.053611343) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.957037) q[1];
sx q[1];
rz(-2.364675) q[1];
sx q[1];
rz(1.0724111) q[1];
rz(-2.8546811) q[3];
sx q[3];
rz(-2.3182456) q[3];
sx q[3];
rz(0.2940184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0733205) q[2];
sx q[2];
rz(-3.1219411) q[2];
sx q[2];
rz(-1.9015296) q[2];
rz(2.9114919) q[3];
sx q[3];
rz(-0.0041882526) q[3];
sx q[3];
rz(2.7219462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0536026) q[0];
sx q[0];
rz(-0.78730655) q[0];
sx q[0];
rz(1.8863652) q[0];
rz(-3.1322196) q[1];
sx q[1];
rz(-1.7724937) q[1];
sx q[1];
rz(3.110041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7566236) q[0];
sx q[0];
rz(-1.0393123) q[0];
sx q[0];
rz(-0.062180659) q[0];
rz(1.8589694) q[2];
sx q[2];
rz(-2.8409578) q[2];
sx q[2];
rz(-0.86821454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42571354) q[1];
sx q[1];
rz(-1.7919994) q[1];
sx q[1];
rz(-2.1137266) q[1];
rz(-pi) q[2];
rz(2.191014) q[3];
sx q[3];
rz(-1.0744541) q[3];
sx q[3];
rz(0.84211189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8072529) q[2];
sx q[2];
rz(-3.1355317) q[2];
sx q[2];
rz(0.80017153) q[2];
rz(2.3470894) q[3];
sx q[3];
rz(-3.1092293) q[3];
sx q[3];
rz(-2.2728424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.054258) q[0];
sx q[0];
rz(-0.15905173) q[0];
sx q[0];
rz(1.6416838) q[0];
rz(-0.17290393) q[1];
sx q[1];
rz(-3.0992295) q[1];
sx q[1];
rz(0.049887966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2931516) q[0];
sx q[0];
rz(-1.3209136) q[0];
sx q[0];
rz(-0.94012733) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4535657) q[2];
sx q[2];
rz(-2.3701326) q[2];
sx q[2];
rz(-1.7457123) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1357291) q[1];
sx q[1];
rz(-1.8195489) q[1];
sx q[1];
rz(0.84524024) q[1];
rz(-pi) q[2];
rz(-0.53277758) q[3];
sx q[3];
rz(-2.2190337) q[3];
sx q[3];
rz(-1.6870354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.738203) q[2];
sx q[2];
rz(-0.047813606) q[2];
sx q[2];
rz(-1.8955463) q[2];
rz(1.3561148) q[3];
sx q[3];
rz(-3.1062283) q[3];
sx q[3];
rz(-1.7252007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7127011) q[0];
sx q[0];
rz(-2.3172947) q[0];
sx q[0];
rz(-1.4225381) q[0];
rz(1.2369583) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(0.2027771) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13064676) q[0];
sx q[0];
rz(-1.8955064) q[0];
sx q[0];
rz(-2.3492011) q[0];
rz(-2.4325763) q[2];
sx q[2];
rz(-1.0142027) q[2];
sx q[2];
rz(-0.32278827) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7084658) q[1];
sx q[1];
rz(-0.34516343) q[1];
sx q[1];
rz(0.16394798) q[1];
rz(0.15722991) q[3];
sx q[3];
rz(-1.0649324) q[3];
sx q[3];
rz(1.0613522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.612959) q[2];
sx q[2];
rz(-0.10047675) q[2];
sx q[2];
rz(2.4670777) q[2];
rz(1.7965192) q[3];
sx q[3];
rz(-0.14480545) q[3];
sx q[3];
rz(1.323918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26789185) q[0];
sx q[0];
rz(-2.38509) q[0];
sx q[0];
rz(-0.85195136) q[0];
rz(0.1918699) q[1];
sx q[1];
rz(-3.1286897) q[1];
sx q[1];
rz(0.26564863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.083355) q[0];
sx q[0];
rz(-1.7449656) q[0];
sx q[0];
rz(2.7252498) q[0];
rz(-pi) q[1];
rz(0.61518367) q[2];
sx q[2];
rz(-1.2968204) q[2];
sx q[2];
rz(-3.0469303) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3959256) q[1];
sx q[1];
rz(-3.0394181) q[1];
sx q[1];
rz(2.5273782) q[1];
x q[2];
rz(-1.9453796) q[3];
sx q[3];
rz(-1.0455971) q[3];
sx q[3];
rz(1.7708667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77136451) q[2];
sx q[2];
rz(-0.092265487) q[2];
sx q[2];
rz(-0.51221687) q[2];
rz(0.15906119) q[3];
sx q[3];
rz(-0.03511196) q[3];
sx q[3];
rz(-1.7062645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8856186) q[0];
sx q[0];
rz(-1.6640478) q[0];
sx q[0];
rz(1.0783827) q[0];
rz(-1.501561) q[1];
sx q[1];
rz(-2.9402132) q[1];
sx q[1];
rz(-1.5837502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5756861) q[0];
sx q[0];
rz(-2.4543336) q[0];
sx q[0];
rz(-3.1172196) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0556765) q[2];
sx q[2];
rz(-2.4011302) q[2];
sx q[2];
rz(-1.6892576) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2557775) q[1];
sx q[1];
rz(-3.1193135) q[1];
sx q[1];
rz(2.9905969) q[1];
rz(1.1224062) q[3];
sx q[3];
rz(-0.70647722) q[3];
sx q[3];
rz(2.024533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.25643361) q[2];
sx q[2];
rz(-3.1270471) q[2];
sx q[2];
rz(0.069615901) q[2];
rz(3.0749248) q[3];
sx q[3];
rz(-2.127141) q[3];
sx q[3];
rz(0.67924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2081864) q[0];
sx q[0];
rz(-1.2766159) q[0];
sx q[0];
rz(2.8896914) q[0];
rz(-1.4725641) q[1];
sx q[1];
rz(-2.9258969) q[1];
sx q[1];
rz(-0.073898166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65739261) q[0];
sx q[0];
rz(-1.9354685) q[0];
sx q[0];
rz(0.075693746) q[0];
rz(-pi) q[1];
rz(-0.017357512) q[2];
sx q[2];
rz(-1.5394754) q[2];
sx q[2];
rz(-0.90085122) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.042797814) q[1];
sx q[1];
rz(-0.86314161) q[1];
sx q[1];
rz(-1.8902197) q[1];
x q[2];
rz(-0.68040397) q[3];
sx q[3];
rz(-0.26016737) q[3];
sx q[3];
rz(-2.5870067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.68022388) q[2];
sx q[2];
rz(-3.1344487) q[2];
sx q[2];
rz(-0.76772493) q[2];
rz(1.7420306) q[3];
sx q[3];
rz(-3.1407686) q[3];
sx q[3];
rz(2.6373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6298228) q[0];
sx q[0];
rz(-0.98291021) q[0];
sx q[0];
rz(1.7194189) q[0];
rz(0.024367532) q[1];
sx q[1];
rz(-0.15932803) q[1];
sx q[1];
rz(-2.9111964) q[1];
rz(1.1047614) q[2];
sx q[2];
rz(-0.73095041) q[2];
sx q[2];
rz(-1.6292844) q[2];
rz(2.7244115) q[3];
sx q[3];
rz(-1.4818896) q[3];
sx q[3];
rz(2.8703574) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
