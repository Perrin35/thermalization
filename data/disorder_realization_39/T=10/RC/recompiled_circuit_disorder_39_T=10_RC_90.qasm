OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(2.6775223) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(1.8571412) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61304898) q[0];
sx q[0];
rz(-1.1819981) q[0];
sx q[0];
rz(0.31624985) q[0];
rz(-pi) q[1];
rz(-2.651865) q[2];
sx q[2];
rz(-1.5416607) q[2];
sx q[2];
rz(-2.4907128) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1661108) q[1];
sx q[1];
rz(-1.2088641) q[1];
sx q[1];
rz(2.6544177) q[1];
x q[2];
rz(-3.0331217) q[3];
sx q[3];
rz(-2.9900108) q[3];
sx q[3];
rz(0.82419318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39711943) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(2.822067) q[2];
rz(-2.5630991) q[3];
sx q[3];
rz(-0.47839034) q[3];
sx q[3];
rz(2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7330866) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(2.615036) q[0];
rz(-0.56354848) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(2.3449576) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54290402) q[0];
sx q[0];
rz(-1.2447378) q[0];
sx q[0];
rz(2.23404) q[0];
rz(0.30226207) q[2];
sx q[2];
rz(-2.5455591) q[2];
sx q[2];
rz(-0.1421393) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4097152) q[1];
sx q[1];
rz(-0.60740031) q[1];
sx q[1];
rz(-3.1190447) q[1];
rz(-pi) q[2];
rz(-0.49780952) q[3];
sx q[3];
rz(-2.1216765) q[3];
sx q[3];
rz(2.2450972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9600296) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(-0.78655085) q[2];
rz(2.6484047) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86984533) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(1.440381) q[0];
rz(-2.4213743) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(2.4386141) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1709135) q[0];
sx q[0];
rz(-1.2909856) q[0];
sx q[0];
rz(0.12165102) q[0];
rz(-pi) q[1];
rz(-1.9269283) q[2];
sx q[2];
rz(-1.0991569) q[2];
sx q[2];
rz(-2.4871662) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.30333334) q[1];
sx q[1];
rz(-1.6175744) q[1];
sx q[1];
rz(-0.10951885) q[1];
rz(-2.968077) q[3];
sx q[3];
rz(-2.1870038) q[3];
sx q[3];
rz(-2.4254352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8748223) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(-2.2276145) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5660969) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(-0.77600586) q[0];
rz(1.874118) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-2.5783096) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8877836) q[0];
sx q[0];
rz(-0.94622181) q[0];
sx q[0];
rz(2.8028691) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4513676) q[2];
sx q[2];
rz(-1.8340655) q[2];
sx q[2];
rz(-3.1379679) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5717585) q[1];
sx q[1];
rz(-2.8697439) q[1];
sx q[1];
rz(-3.1175201) q[1];
x q[2];
rz(2.4241583) q[3];
sx q[3];
rz(-1.7955901) q[3];
sx q[3];
rz(-0.7625398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6528066) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(2.9768067) q[2];
rz(2.9131043) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27424681) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(-0.85246032) q[0];
rz(-0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(-1.16211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1253818) q[0];
sx q[0];
rz(-2.1124766) q[0];
sx q[0];
rz(1.1075695) q[0];
rz(-2.7439793) q[2];
sx q[2];
rz(-2.0587454) q[2];
sx q[2];
rz(-2.5401126) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0218378) q[1];
sx q[1];
rz(-0.86958414) q[1];
sx q[1];
rz(-0.55418684) q[1];
rz(-pi) q[2];
x q[2];
rz(0.01708548) q[3];
sx q[3];
rz(-0.72165976) q[3];
sx q[3];
rz(-2.7929896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.44624415) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(1.2472786) q[2];
rz(-0.013109664) q[3];
sx q[3];
rz(-1.40991) q[3];
sx q[3];
rz(2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(-2.2728249) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(1.3060588) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36201492) q[0];
sx q[0];
rz(-1.0581731) q[0];
sx q[0];
rz(-0.4431475) q[0];
x q[1];
rz(-1.3681709) q[2];
sx q[2];
rz(-1.1154004) q[2];
sx q[2];
rz(2.1855598) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8514511) q[1];
sx q[1];
rz(-2.4212004) q[1];
sx q[1];
rz(0.011522567) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4218876) q[3];
sx q[3];
rz(-2.1544666) q[3];
sx q[3];
rz(-2.7910809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7541472) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(-0.60069096) q[2];
rz(2.1389652) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7464741) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(1.612161) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(-0.41710645) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0102651) q[0];
sx q[0];
rz(-1.7908887) q[0];
sx q[0];
rz(0.025028153) q[0];
rz(-pi) q[1];
rz(0.30714005) q[2];
sx q[2];
rz(-1.8516314) q[2];
sx q[2];
rz(-1.3020696) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.68361359) q[1];
sx q[1];
rz(-0.90404592) q[1];
sx q[1];
rz(-0.97138202) q[1];
x q[2];
rz(-0.38325558) q[3];
sx q[3];
rz(-2.4517422) q[3];
sx q[3];
rz(-2.3243429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0059011857) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(-2.588429) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(-2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(-1.1897855) q[0];
rz(-1.7143543) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(-3.0292125) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0880786) q[0];
sx q[0];
rz(-1.7435762) q[0];
sx q[0];
rz(2.0515576) q[0];
rz(-pi) q[1];
rz(0.33371146) q[2];
sx q[2];
rz(-1.6534272) q[2];
sx q[2];
rz(-2.4090648) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97265128) q[1];
sx q[1];
rz(-2.2691138) q[1];
sx q[1];
rz(-2.6886743) q[1];
rz(1.8295248) q[3];
sx q[3];
rz(-1.4530164) q[3];
sx q[3];
rz(-0.044101276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(1.3486264) q[2];
rz(-1.9366692) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(3.1066185) q[0];
rz(0.84683013) q[1];
sx q[1];
rz(-1.433082) q[1];
sx q[1];
rz(-0.91167489) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.002279) q[0];
sx q[0];
rz(-0.3542491) q[0];
sx q[0];
rz(2.2646963) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20810017) q[2];
sx q[2];
rz(-0.55317438) q[2];
sx q[2];
rz(-2.2262239) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94757838) q[1];
sx q[1];
rz(-1.7742426) q[1];
sx q[1];
rz(-0.46781637) q[1];
rz(-0.73260143) q[3];
sx q[3];
rz(-2.4676975) q[3];
sx q[3];
rz(-3.1042838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70790616) q[2];
sx q[2];
rz(-0.5138548) q[2];
sx q[2];
rz(-0.24924499) q[2];
rz(0.76672673) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3578167) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(-2.7695079) q[0];
rz(-0.58139873) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(1.4153597) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5370731) q[0];
sx q[0];
rz(-1.8329289) q[0];
sx q[0];
rz(1.7781236) q[0];
rz(-pi) q[1];
x q[1];
rz(2.084311) q[2];
sx q[2];
rz(-1.0143177) q[2];
sx q[2];
rz(1.2531812) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.806554) q[1];
sx q[1];
rz(-1.5015258) q[1];
sx q[1];
rz(0.075242234) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7578027) q[3];
sx q[3];
rz(-2.0818315) q[3];
sx q[3];
rz(-0.25110652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.32594484) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(-2.9369205) q[2];
rz(1.4137911) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(2.0457101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50080147) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(1.5851371) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(2.0251705) q[2];
sx q[2];
rz(-1.4618256) q[2];
sx q[2];
rz(1.7088919) q[2];
rz(-3.0588991) q[3];
sx q[3];
rz(-2.1413998) q[3];
sx q[3];
rz(-1.0257046) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
