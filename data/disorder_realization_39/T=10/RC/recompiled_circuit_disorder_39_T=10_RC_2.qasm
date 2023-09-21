OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0948148) q[0];
sx q[0];
rz(4.2098213) q[0];
sx q[0];
rz(9.8888483) q[0];
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
rz(0.099079236) q[0];
sx q[0];
rz(-0.49603841) q[0];
sx q[0];
rz(2.2201559) q[0];
rz(0.48972763) q[2];
sx q[2];
rz(-1.5999319) q[2];
sx q[2];
rz(-0.65087989) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40988906) q[1];
sx q[1];
rz(-1.1176425) q[1];
sx q[1];
rz(-1.975592) q[1];
rz(-0.10847096) q[3];
sx q[3];
rz(-0.15158187) q[3];
sx q[3];
rz(-2.3173995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7444732) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(0.31952566) q[2];
rz(0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(0.67392504) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7330866) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(-2.615036) q[0];
rz(-0.56354848) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(2.3449576) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63882534) q[0];
sx q[0];
rz(-0.728038) q[0];
sx q[0];
rz(1.0685705) q[0];
x q[1];
rz(1.37155) q[2];
sx q[2];
rz(-1.0052048) q[2];
sx q[2];
rz(2.9233962) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7318774) q[1];
sx q[1];
rz(-2.5341923) q[1];
sx q[1];
rz(-0.022547988) q[1];
rz(0.96062406) q[3];
sx q[3];
rz(-1.1517797) q[3];
sx q[3];
rz(0.39715365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9600296) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(0.78655085) q[2];
rz(-0.49318796) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(-2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86984533) q[0];
sx q[0];
rz(-0.87690502) q[0];
sx q[0];
rz(-1.440381) q[0];
rz(-2.4213743) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(-0.70297855) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7754606) q[0];
sx q[0];
rz(-1.6876939) q[0];
sx q[0];
rz(1.2890105) q[0];
rz(0.59962745) q[2];
sx q[2];
rz(-2.5587974) q[2];
sx q[2];
rz(-1.3404913) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30333334) q[1];
sx q[1];
rz(-1.5240182) q[1];
sx q[1];
rz(3.0320738) q[1];
rz(-pi) q[2];
rz(-0.17351563) q[3];
sx q[3];
rz(-0.95458889) q[3];
sx q[3];
rz(0.71615744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8748223) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(0.91153574) q[2];
rz(-2.2276145) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(-0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754958) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(0.77600586) q[0];
rz(1.874118) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(0.56328303) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8877836) q[0];
sx q[0];
rz(-0.94622181) q[0];
sx q[0];
rz(-2.8028691) q[0];
x q[1];
rz(2.4513676) q[2];
sx q[2];
rz(-1.8340655) q[2];
sx q[2];
rz(-3.1379679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5717585) q[1];
sx q[1];
rz(-0.27184871) q[1];
sx q[1];
rz(-0.024072577) q[1];
rz(-pi) q[2];
rz(0.71743439) q[3];
sx q[3];
rz(-1.7955901) q[3];
sx q[3];
rz(0.7625398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6528066) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(-0.164786) q[2];
rz(0.22848836) q[3];
sx q[3];
rz(-0.27992862) q[3];
sx q[3];
rz(2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27424681) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(-2.2891323) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(1.16211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4441372) q[0];
sx q[0];
rz(-1.9636969) q[0];
sx q[0];
rz(0.59209728) q[0];
rz(-1.048462) q[2];
sx q[2];
rz(-1.221721) q[2];
sx q[2];
rz(2.366684) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.64486849) q[1];
sx q[1];
rz(-2.2779896) q[1];
sx q[1];
rz(2.1281388) q[1];
rz(-0.72158738) q[3];
sx q[3];
rz(-1.5820832) q[3];
sx q[3];
rz(1.9065726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6953485) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(-1.8943141) q[2];
rz(3.128483) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(-2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86876774) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(-1.8355339) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7346142) q[0];
sx q[0];
rz(-2.4771871) q[0];
sx q[0];
rz(0.91974308) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38988955) q[2];
sx q[2];
rz(-2.6460558) q[2];
sx q[2];
rz(1.3930266) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8522779) q[1];
sx q[1];
rz(-1.5783974) q[1];
sx q[1];
rz(-0.72035933) q[1];
x q[2];
rz(1.7197051) q[3];
sx q[3];
rz(-2.1544666) q[3];
sx q[3];
rz(-2.7910809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7541472) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(-0.60069096) q[2];
rz(2.1389652) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(-2.7887662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7464741) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(-2.0948998) q[0];
rz(-1.612161) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(-2.7244862) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0171623) q[0];
sx q[0];
rz(-2.9201047) q[0];
sx q[0];
rz(1.4593967) q[0];
x q[1];
rz(2.3796758) q[2];
sx q[2];
rz(-0.41315213) q[2];
sx q[2];
rz(2.6921536) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.48737803) q[1];
sx q[1];
rz(-1.1114792) q[1];
sx q[1];
rz(-0.76141255) q[1];
rz(-pi) q[2];
rz(0.38325558) q[3];
sx q[3];
rz(-2.4517422) q[3];
sx q[3];
rz(-0.81724973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1356915) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(-0.55316365) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(-2.6161391) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.816514) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(1.9518071) q[0];
rz(1.4272383) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(0.11238012) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0880786) q[0];
sx q[0];
rz(-1.7435762) q[0];
sx q[0];
rz(-1.0900351) q[0];
x q[1];
rz(2.8939395) q[2];
sx q[2];
rz(-2.7981749) q[2];
sx q[2];
rz(-2.0695956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97265128) q[1];
sx q[1];
rz(-2.2691138) q[1];
sx q[1];
rz(0.45291839) q[1];
rz(-pi) q[2];
rz(-0.12179575) q[3];
sx q[3];
rz(-1.8276916) q[3];
sx q[3];
rz(1.6459873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0743951) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(-1.3486264) q[2];
rz(1.9366692) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(-0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(-3.1066185) q[0];
rz(-2.2947625) q[1];
sx q[1];
rz(-1.433082) q[1];
sx q[1];
rz(2.2299178) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23110403) q[0];
sx q[0];
rz(-1.7945053) q[0];
sx q[0];
rz(1.847812) q[0];
rz(-pi) q[1];
rz(-0.20810017) q[2];
sx q[2];
rz(-2.5884183) q[2];
sx q[2];
rz(-0.91536872) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8988077) q[1];
sx q[1];
rz(-0.50711942) q[1];
sx q[1];
rz(0.42906638) q[1];
rz(-pi) q[2];
rz(1.0802286) q[3];
sx q[3];
rz(-1.0883696) q[3];
sx q[3];
rz(2.3232943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(2.8923477) q[2];
rz(-2.3748659) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3578167) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(2.7695079) q[0];
rz(-0.58139873) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(-1.4153597) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5370731) q[0];
sx q[0];
rz(-1.8329289) q[0];
sx q[0];
rz(1.7781236) q[0];
rz(-2.5214228) q[2];
sx q[2];
rz(-2.0010741) q[2];
sx q[2];
rz(-0.60713965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6485939) q[1];
sx q[1];
rz(-0.10222888) q[1];
sx q[1];
rz(0.74536721) q[1];
rz(-pi) q[2];
rz(1.0269208) q[3];
sx q[3];
rz(-1.9034991) q[3];
sx q[3];
rz(-1.1247016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32594484) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(0.20467219) q[2];
rz(-1.4137911) q[3];
sx q[3];
rz(-1.1433583) q[3];
sx q[3];
rz(2.0457101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50080147) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(-1.5564556) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(3.0204308) q[2];
sx q[2];
rz(-1.1193174) q[2];
sx q[2];
rz(0.08502273) q[2];
rz(-0.08269357) q[3];
sx q[3];
rz(-1.0001928) q[3];
sx q[3];
rz(2.115888) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
