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
rz(2.6858202) q[0];
sx q[0];
rz(-2.0210285) q[0];
sx q[0];
rz(-1.4135345) q[0];
rz(-2.7780374) q[1];
sx q[1];
rz(-1.8283365) q[1];
sx q[1];
rz(0.29120905) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094452642) q[0];
sx q[0];
rz(-1.5507076) q[0];
sx q[0];
rz(1.8909341) q[0];
x q[1];
rz(-1.9092567) q[2];
sx q[2];
rz(-2.7833653) q[2];
sx q[2];
rz(-0.097964935) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6068961) q[1];
sx q[1];
rz(-1.7217858) q[1];
sx q[1];
rz(0.57735898) q[1];
rz(-pi) q[2];
rz(-2.324027) q[3];
sx q[3];
rz(-2.052784) q[3];
sx q[3];
rz(0.40460247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1353961) q[2];
sx q[2];
rz(-1.5781032) q[2];
sx q[2];
rz(-3.1175933) q[2];
rz(2.7855347) q[3];
sx q[3];
rz(-2.4672716) q[3];
sx q[3];
rz(-3.1168028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085670797) q[0];
sx q[0];
rz(-2.4503158) q[0];
sx q[0];
rz(2.5066277) q[0];
rz(2.3786646) q[1];
sx q[1];
rz(-1.678391) q[1];
sx q[1];
rz(-0.91711226) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0046154) q[0];
sx q[0];
rz(-0.03923035) q[0];
sx q[0];
rz(-2.0580388) q[0];
rz(-pi) q[1];
rz(2.6659662) q[2];
sx q[2];
rz(-2.1948819) q[2];
sx q[2];
rz(-1.2211354) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.13027993) q[1];
sx q[1];
rz(-2.6940799) q[1];
sx q[1];
rz(-1.7496197) q[1];
x q[2];
rz(0.21288721) q[3];
sx q[3];
rz(-2.1799984) q[3];
sx q[3];
rz(1.4573801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.94747535) q[2];
sx q[2];
rz(-2.6185991) q[2];
sx q[2];
rz(-2.4083162) q[2];
rz(-2.0745847) q[3];
sx q[3];
rz(-1.4209483) q[3];
sx q[3];
rz(-2.9616621) q[3];
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
rz(0.5697524) q[0];
sx q[0];
rz(-1.952992) q[0];
sx q[0];
rz(0.52823129) q[0];
rz(-2.275548) q[1];
sx q[1];
rz(-1.8128315) q[1];
sx q[1];
rz(-2.6285062) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7932803) q[0];
sx q[0];
rz(-1.5282807) q[0];
sx q[0];
rz(-2.0094064) q[0];
rz(-pi) q[1];
rz(1.971286) q[2];
sx q[2];
rz(-0.68953625) q[2];
sx q[2];
rz(0.51338085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9218623) q[1];
sx q[1];
rz(-0.37223909) q[1];
sx q[1];
rz(-1.9277186) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5993027) q[3];
sx q[3];
rz(-2.0703201) q[3];
sx q[3];
rz(-2.0716425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.036006007) q[2];
sx q[2];
rz(-1.0893818) q[2];
sx q[2];
rz(3.1028683) q[2];
rz(0.45768467) q[3];
sx q[3];
rz(-0.80515146) q[3];
sx q[3];
rz(1.3409486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4337346) q[0];
sx q[0];
rz(-2.7271294) q[0];
sx q[0];
rz(1.1210972) q[0];
rz(-2.3838249) q[1];
sx q[1];
rz(-1.517375) q[1];
sx q[1];
rz(-2.4574492) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4829041) q[0];
sx q[0];
rz(-1.7724121) q[0];
sx q[0];
rz(-2.3134507) q[0];
rz(-pi) q[1];
rz(2.5252417) q[2];
sx q[2];
rz(-0.01238298) q[2];
sx q[2];
rz(-1.826926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5110179) q[1];
sx q[1];
rz(-1.8749494) q[1];
sx q[1];
rz(1.8777147) q[1];
rz(-pi) q[2];
rz(0.87881256) q[3];
sx q[3];
rz(-2.4185816) q[3];
sx q[3];
rz(-2.3634152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6571558) q[2];
sx q[2];
rz(-1.5269205) q[2];
sx q[2];
rz(0.26051513) q[2];
rz(-0.83812964) q[3];
sx q[3];
rz(-1.1476436) q[3];
sx q[3];
rz(3.1269791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92358661) q[0];
sx q[0];
rz(-1.3777233) q[0];
sx q[0];
rz(-2.6309784) q[0];
rz(-1.249373) q[1];
sx q[1];
rz(-1.1621954) q[1];
sx q[1];
rz(1.4543264) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.175486) q[0];
sx q[0];
rz(-1.5990055) q[0];
sx q[0];
rz(-2.9278446) q[0];
rz(-2.8122105) q[2];
sx q[2];
rz(-0.14757809) q[2];
sx q[2];
rz(-1.7640698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1087423) q[1];
sx q[1];
rz(-1.8287484) q[1];
sx q[1];
rz(0.56405073) q[1];
rz(-0.14438676) q[3];
sx q[3];
rz(-1.4098415) q[3];
sx q[3];
rz(-0.99100366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.55383596) q[2];
sx q[2];
rz(-0.38148701) q[2];
sx q[2];
rz(-0.39056632) q[2];
rz(0.25911123) q[3];
sx q[3];
rz(-1.5026389) q[3];
sx q[3];
rz(2.794877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3542341) q[0];
sx q[0];
rz(-2.7095095) q[0];
sx q[0];
rz(-2.4969192) q[0];
rz(1.8395754) q[1];
sx q[1];
rz(-1.6141067) q[1];
sx q[1];
rz(0.34861809) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2520599) q[0];
sx q[0];
rz(-2.4675998) q[0];
sx q[0];
rz(0.78788449) q[0];
rz(-pi) q[1];
rz(2.5041144) q[2];
sx q[2];
rz(-0.63796746) q[2];
sx q[2];
rz(-1.5359985) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4233089) q[1];
sx q[1];
rz(-1.8656785) q[1];
sx q[1];
rz(-2.0096094) q[1];
x q[2];
rz(1.5362605) q[3];
sx q[3];
rz(-1.6876965) q[3];
sx q[3];
rz(1.7021321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0384486) q[2];
sx q[2];
rz(-2.1861173) q[2];
sx q[2];
rz(3.1032584) q[2];
rz(2.1740055) q[3];
sx q[3];
rz(-1.7318232) q[3];
sx q[3];
rz(-2.7819395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86647812) q[0];
sx q[0];
rz(-2.622128) q[0];
sx q[0];
rz(-0.073609784) q[0];
rz(1.4069125) q[1];
sx q[1];
rz(-2.584447) q[1];
sx q[1];
rz(0.73307347) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778962) q[0];
sx q[0];
rz(-1.0390345) q[0];
sx q[0];
rz(-0.97490208) q[0];
x q[1];
rz(0.91821155) q[2];
sx q[2];
rz(-0.84647734) q[2];
sx q[2];
rz(1.0670964) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0907369) q[1];
sx q[1];
rz(-2.6466718) q[1];
sx q[1];
rz(-0.75116091) q[1];
rz(0.81871512) q[3];
sx q[3];
rz(-2.8041841) q[3];
sx q[3];
rz(-0.81160566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1072032) q[2];
sx q[2];
rz(-1.723039) q[2];
sx q[2];
rz(0.011073152) q[2];
rz(-2.8360046) q[3];
sx q[3];
rz(-1.1253076) q[3];
sx q[3];
rz(-1.8886458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8082064) q[0];
sx q[0];
rz(-0.60788637) q[0];
sx q[0];
rz(2.407684) q[0];
rz(-0.58087307) q[1];
sx q[1];
rz(-1.0092694) q[1];
sx q[1];
rz(-0.098800585) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0694612) q[0];
sx q[0];
rz(-1.0476351) q[0];
sx q[0];
rz(-0.69055542) q[0];
rz(-0.67811857) q[2];
sx q[2];
rz(-1.0762012) q[2];
sx q[2];
rz(1.2505609) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96896763) q[1];
sx q[1];
rz(-1.8207112) q[1];
sx q[1];
rz(-2.6843798) q[1];
x q[2];
rz(-1.4905761) q[3];
sx q[3];
rz(-1.7003683) q[3];
sx q[3];
rz(0.097565325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8431479) q[2];
sx q[2];
rz(-1.2848022) q[2];
sx q[2];
rz(-1.0996381) q[2];
rz(-3.1067276) q[3];
sx q[3];
rz(-2.3979082) q[3];
sx q[3];
rz(-2.5449424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1140887) q[0];
sx q[0];
rz(-1.1570258) q[0];
sx q[0];
rz(-1.9644568) q[0];
rz(1.4741395) q[1];
sx q[1];
rz(-0.93282229) q[1];
sx q[1];
rz(1.4449545) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.225017) q[0];
sx q[0];
rz(-1.1780043) q[0];
sx q[0];
rz(-1.265822) q[0];
rz(-pi) q[1];
rz(2.1423526) q[2];
sx q[2];
rz(-0.88724595) q[2];
sx q[2];
rz(2.7627166) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18981537) q[1];
sx q[1];
rz(-0.51945247) q[1];
sx q[1];
rz(-1.7096667) q[1];
rz(0.8952462) q[3];
sx q[3];
rz(-1.352868) q[3];
sx q[3];
rz(-2.982614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.16704796) q[2];
sx q[2];
rz(-1.5844774) q[2];
sx q[2];
rz(2.9696999) q[2];
rz(0.78957549) q[3];
sx q[3];
rz(-2.004345) q[3];
sx q[3];
rz(-1.8062228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0563141) q[0];
sx q[0];
rz(-2.9264937) q[0];
sx q[0];
rz(2.6000182) q[0];
rz(-1.9821292) q[1];
sx q[1];
rz(-1.3075202) q[1];
sx q[1];
rz(-2.3084739) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.706382) q[0];
sx q[0];
rz(-2.4738389) q[0];
sx q[0];
rz(0.45489648) q[0];
rz(-pi) q[1];
rz(-1.2255185) q[2];
sx q[2];
rz(-2.08689) q[2];
sx q[2];
rz(2.6654625) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.009506735) q[1];
sx q[1];
rz(-0.888538) q[1];
sx q[1];
rz(-0.29341523) q[1];
rz(-pi) q[2];
rz(-1.572979) q[3];
sx q[3];
rz(-2.4139521) q[3];
sx q[3];
rz(1.997987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1200166) q[2];
sx q[2];
rz(-2.3638201) q[2];
sx q[2];
rz(2.1982101) q[2];
rz(1.0776445) q[3];
sx q[3];
rz(-1.5871983) q[3];
sx q[3];
rz(0.32850346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42644603) q[0];
sx q[0];
rz(-2.1052512) q[0];
sx q[0];
rz(-0.25682009) q[0];
rz(2.9019451) q[1];
sx q[1];
rz(-0.9050723) q[1];
sx q[1];
rz(-1.1368652) q[1];
rz(-2.4113833) q[2];
sx q[2];
rz(-2.3013023) q[2];
sx q[2];
rz(-0.58462044) q[2];
rz(-0.67962525) q[3];
sx q[3];
rz(-0.39120318) q[3];
sx q[3];
rz(-0.92728566) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
