OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5840924) q[0];
sx q[0];
rz(-0.023356181) q[0];
sx q[0];
rz(-2.2060702) q[0];
rz(-1.6159396) q[1];
sx q[1];
rz(-1.5204117) q[1];
sx q[1];
rz(2.867155) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68831681) q[0];
sx q[0];
rz(-3.0429646) q[0];
sx q[0];
rz(-2.3345956) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53977592) q[2];
sx q[2];
rz(-1.1640884) q[2];
sx q[2];
rz(-2.0222029) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86286592) q[1];
sx q[1];
rz(-0.037986156) q[1];
sx q[1];
rz(-1.2943062) q[1];
rz(-pi) q[2];
rz(-1.9755483) q[3];
sx q[3];
rz(-0.18530986) q[3];
sx q[3];
rz(-2.6626793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9258257) q[2];
sx q[2];
rz(-0.01130686) q[2];
sx q[2];
rz(-2.0634148) q[2];
rz(0.82873851) q[3];
sx q[3];
rz(-1.5107892) q[3];
sx q[3];
rz(2.3327648) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082212903) q[0];
sx q[0];
rz(-1.8635211) q[0];
sx q[0];
rz(-1.7510121) q[0];
rz(-1.7104205) q[1];
sx q[1];
rz(-0.0043967604) q[1];
sx q[1];
rz(1.4336525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19884685) q[0];
sx q[0];
rz(-2.2042243) q[0];
sx q[0];
rz(3.01157) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1251034) q[2];
sx q[2];
rz(-1.1270583) q[2];
sx q[2];
rz(1.6042142) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8987309) q[1];
sx q[1];
rz(-0.0090829385) q[1];
sx q[1];
rz(-1.6159978) q[1];
x q[2];
rz(0.78761132) q[3];
sx q[3];
rz(-1.6171086) q[3];
sx q[3];
rz(2.5758343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7626875) q[2];
sx q[2];
rz(-1.5964369) q[2];
sx q[2];
rz(1.5691266) q[2];
rz(-0.76111859) q[3];
sx q[3];
rz(-0.053839024) q[3];
sx q[3];
rz(-1.9332473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91956562) q[0];
sx q[0];
rz(-0.61028218) q[0];
sx q[0];
rz(-0.56184226) q[0];
rz(-1.5678844) q[1];
sx q[1];
rz(-1.578873) q[1];
sx q[1];
rz(-3.124253) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0837096) q[0];
sx q[0];
rz(-2.0590933) q[0];
sx q[0];
rz(-0.47476193) q[0];
x q[1];
rz(-2.2566913) q[2];
sx q[2];
rz(-1.1807311) q[2];
sx q[2];
rz(1.6194413) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8120809) q[1];
sx q[1];
rz(-1.9332262) q[1];
sx q[1];
rz(3.1245438) q[1];
rz(-pi) q[2];
x q[2];
rz(0.041009237) q[3];
sx q[3];
rz(-2.6309391) q[3];
sx q[3];
rz(1.8927495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6277546) q[2];
sx q[2];
rz(-1.4189812) q[2];
sx q[2];
rz(-2.5886152) q[2];
rz(-1.9364457) q[3];
sx q[3];
rz(-1.5744934) q[3];
sx q[3];
rz(-1.5470101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.39128458) q[0];
sx q[0];
rz(-1.0306083) q[0];
sx q[0];
rz(-1.3543825) q[0];
rz(-1.7693819) q[1];
sx q[1];
rz(-0.0028227614) q[1];
sx q[1];
rz(1.3815968) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3171661) q[0];
sx q[0];
rz(-1.1689416) q[0];
sx q[0];
rz(-2.0749932) q[0];
rz(-pi) q[1];
rz(-2.4468378) q[2];
sx q[2];
rz(-0.0036029795) q[2];
sx q[2];
rz(-0.82243332) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.957037) q[1];
sx q[1];
rz(-2.364675) q[1];
sx q[1];
rz(2.0691815) q[1];
rz(-1.2744585) q[3];
sx q[3];
rz(-2.3510286) q[3];
sx q[3];
rz(2.4380655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0733205) q[2];
sx q[2];
rz(-0.019651532) q[2];
sx q[2];
rz(1.9015296) q[2];
rz(0.23010075) q[3];
sx q[3];
rz(-0.0041882526) q[3];
sx q[3];
rz(0.41964644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0536026) q[0];
sx q[0];
rz(-2.3542861) q[0];
sx q[0];
rz(1.8863652) q[0];
rz(-3.1322196) q[1];
sx q[1];
rz(-1.3690989) q[1];
sx q[1];
rz(-3.110041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2627336) q[0];
sx q[0];
rz(-0.53476214) q[0];
sx q[0];
rz(-1.6761024) q[0];
x q[1];
rz(1.8589694) q[2];
sx q[2];
rz(-2.8409578) q[2];
sx q[2];
rz(2.2733781) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0134558) q[1];
sx q[1];
rz(-2.0990879) q[1];
sx q[1];
rz(-2.8847413) q[1];
x q[2];
rz(2.3208951) q[3];
sx q[3];
rz(-0.77337556) q[3];
sx q[3];
rz(-0.14062961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.3343398) q[2];
sx q[2];
rz(-0.0060609239) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0873347) q[0];
sx q[0];
rz(-2.9825409) q[0];
sx q[0];
rz(-1.4999088) q[0];
rz(0.17290393) q[1];
sx q[1];
rz(-3.0992295) q[1];
sx q[1];
rz(-0.049887966) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.04258) q[0];
sx q[0];
rz(-2.178971) q[0];
sx q[0];
rz(2.8355168) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11325963) q[2];
sx q[2];
rz(-0.80598968) q[2];
sx q[2];
rz(1.2330556) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7916471) q[1];
sx q[1];
rz(-0.87222067) q[1];
sx q[1];
rz(-2.8142745) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53277758) q[3];
sx q[3];
rz(-0.92255892) q[3];
sx q[3];
rz(1.6870354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40338966) q[2];
sx q[2];
rz(-3.093779) q[2];
sx q[2];
rz(1.8955463) q[2];
rz(1.7854779) q[3];
sx q[3];
rz(-3.1062283) q[3];
sx q[3];
rz(1.7252007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7127011) q[0];
sx q[0];
rz(-2.3172947) q[0];
sx q[0];
rz(-1.7190546) q[0];
rz(1.9046344) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(2.9388156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7449489) q[0];
sx q[0];
rz(-2.2988964) q[0];
sx q[0];
rz(0.44162314) q[0];
rz(0.70901633) q[2];
sx q[2];
rz(-2.1273899) q[2];
sx q[2];
rz(0.32278827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2591035) q[1];
sx q[1];
rz(-1.2304502) q[1];
sx q[1];
rz(1.5121786) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2953128) q[3];
sx q[3];
rz(-2.6138895) q[3];
sx q[3];
rz(0.74515158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.612959) q[2];
sx q[2];
rz(-3.0411159) q[2];
sx q[2];
rz(0.67451492) q[2];
rz(-1.3450735) q[3];
sx q[3];
rz(-0.14480545) q[3];
sx q[3];
rz(-1.8176746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26789185) q[0];
sx q[0];
rz(-0.75650263) q[0];
sx q[0];
rz(-2.2896413) q[0];
rz(-2.9497228) q[1];
sx q[1];
rz(-3.1286897) q[1];
sx q[1];
rz(0.26564863) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56391956) q[0];
sx q[0];
rz(-1.980459) q[0];
sx q[0];
rz(-1.760861) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2393565) q[2];
sx q[2];
rz(-0.98173117) q[2];
sx q[2];
rz(1.6650496) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.354968) q[1];
sx q[1];
rz(-1.6296128) q[1];
sx q[1];
rz(3.0579964) q[1];
rz(-pi) q[2];
rz(0.55691584) q[3];
sx q[3];
rz(-1.2486826) q[3];
sx q[3];
rz(-2.7469001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8856186) q[0];
sx q[0];
rz(-1.6640478) q[0];
sx q[0];
rz(2.0632099) q[0];
rz(1.501561) q[1];
sx q[1];
rz(-2.9402132) q[1];
sx q[1];
rz(1.5837502) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5659065) q[0];
sx q[0];
rz(-0.68725902) q[0];
sx q[0];
rz(-0.024373011) q[0];
x q[1];
rz(-0.89084004) q[2];
sx q[2];
rz(-1.8906697) q[2];
sx q[2];
rz(0.25242463) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4068102) q[1];
sx q[1];
rz(-1.592822) q[1];
sx q[1];
rz(1.5674445) q[1];
rz(-pi) q[2];
rz(2.0191865) q[3];
sx q[3];
rz(-2.4351154) q[3];
sx q[3];
rz(-1.1170596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.885159) q[2];
sx q[2];
rz(-3.1270471) q[2];
sx q[2];
rz(-3.0719768) q[2];
rz(3.0749248) q[3];
sx q[3];
rz(-1.0144517) q[3];
sx q[3];
rz(-0.67924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93340623) q[0];
sx q[0];
rz(-1.8649768) q[0];
sx q[0];
rz(0.25190121) q[0];
rz(1.4725641) q[1];
sx q[1];
rz(-0.21569574) q[1];
sx q[1];
rz(3.0676945) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2011482) q[0];
sx q[0];
rz(-1.5000888) q[0];
sx q[0];
rz(-1.9364249) q[0];
rz(-2.0766856) q[2];
sx q[2];
rz(-3.1057851) q[2];
sx q[2];
rz(-0.39469013) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.51346362) q[1];
sx q[1];
rz(-2.3766915) q[1];
sx q[1];
rz(2.7897863) q[1];
rz(-2.9375495) q[3];
sx q[3];
rz(-1.7333441) q[3];
sx q[3];
rz(-0.35248392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4613688) q[2];
sx q[2];
rz(-3.1344487) q[2];
sx q[2];
rz(-2.3738677) q[2];
rz(-1.399562) q[3];
sx q[3];
rz(-3.1407686) q[3];
sx q[3];
rz(-0.50423938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
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
rz(-2.2461535) q[2];
sx q[2];
rz(-1.2661305) q[2];
sx q[2];
rz(0.2998395) q[2];
rz(1.4735994) q[3];
sx q[3];
rz(-1.9862277) q[3];
sx q[3];
rz(1.338892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
