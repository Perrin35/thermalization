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
rz(-1.2844515) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5285437) q[0];
sx q[0];
rz(-1.1819981) q[0];
sx q[0];
rz(-0.31624985) q[0];
rz(-2.651865) q[2];
sx q[2];
rz(-1.5999319) q[2];
sx q[2];
rz(2.4907128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1661108) q[1];
sx q[1];
rz(-1.9327285) q[1];
sx q[1];
rz(-0.48717498) q[1];
rz(-pi) q[2];
rz(-2.9908882) q[3];
sx q[3];
rz(-1.5871443) q[3];
sx q[3];
rz(2.5022262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7444732) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(-0.31952566) q[2];
rz(0.5784936) q[3];
sx q[3];
rz(-0.47839034) q[3];
sx q[3];
rz(2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7330866) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(-2.615036) q[0];
rz(0.56354848) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(2.3449576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63882534) q[0];
sx q[0];
rz(-2.4135547) q[0];
sx q[0];
rz(-1.0685705) q[0];
rz(-pi) q[1];
rz(-2.566922) q[2];
sx q[2];
rz(-1.4029014) q[2];
sx q[2];
rz(-1.6811973) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4097152) q[1];
sx q[1];
rz(-0.60740031) q[1];
sx q[1];
rz(-0.022547988) q[1];
rz(-2.1809686) q[3];
sx q[3];
rz(-1.1517797) q[3];
sx q[3];
rz(-2.744439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.18156302) q[2];
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
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86984533) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(-1.7012117) q[0];
rz(-2.4213743) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(-0.70297855) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5876578) q[0];
sx q[0];
rz(-0.30447391) q[0];
sx q[0];
rz(1.1712043) q[0];
rz(-0.59962745) q[2];
sx q[2];
rz(-2.5587974) q[2];
sx q[2];
rz(1.3404913) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2726047) q[1];
sx q[1];
rz(-1.4613978) q[1];
sx q[1];
rz(-1.617856) q[1];
rz(-pi) q[2];
rz(-1.8099144) q[3];
sx q[3];
rz(-2.5044887) q[3];
sx q[3];
rz(1.0106196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8748223) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(-2.2300569) q[2];
rz(0.91397816) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5660969) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(-0.77600586) q[0];
rz(1.2674747) q[1];
sx q[1];
rz(-2.0327366) q[1];
sx q[1];
rz(0.56328303) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3455428) q[0];
sx q[0];
rz(-2.4420218) q[0];
sx q[0];
rz(-2.0027341) q[0];
rz(-0.40043719) q[2];
sx q[2];
rz(-2.4106328) q[2];
sx q[2];
rz(1.8797344) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54484425) q[1];
sx q[1];
rz(-1.8425643) q[1];
sx q[1];
rz(1.5775058) q[1];
x q[2];
rz(-0.33470811) q[3];
sx q[3];
rz(-2.3957806) q[3];
sx q[3];
rz(-2.0832182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6528066) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(2.9768067) q[2];
rz(2.9131043) q[3];
sx q[3];
rz(-0.27992862) q[3];
sx q[3];
rz(0.94648186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27424681) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(0.85246032) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(1.16211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64415414) q[0];
sx q[0];
rz(-0.69735202) q[0];
sx q[0];
rz(2.5028412) q[0];
rz(-2.2010872) q[2];
sx q[2];
rz(-0.61912196) q[2];
sx q[2];
rz(-1.3319912) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4967242) q[1];
sx q[1];
rz(-0.86360303) q[1];
sx q[1];
rz(2.1281388) q[1];
rz(-pi) q[2];
rz(1.5557628) q[3];
sx q[3];
rz(-2.2923277) q[3];
sx q[3];
rz(2.815747) q[3];
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
rz(0.013109664) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(-0.75063467) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(1.3060588) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7040946) q[0];
sx q[0];
rz(-1.9537582) q[0];
sx q[0];
rz(-2.1279446) q[0];
x q[1];
rz(-1.3681709) q[2];
sx q[2];
rz(-2.0261923) q[2];
sx q[2];
rz(0.95603285) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8522779) q[1];
sx q[1];
rz(-1.5631952) q[1];
sx q[1];
rz(-2.4212333) q[1];
rz(2.9206198) q[3];
sx q[3];
rz(-0.60022012) q[3];
sx q[3];
rz(2.5252987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3874454) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(-2.5409017) q[2];
rz(1.0026275) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(-2.7887662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-1.7464741) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(-1.612161) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(2.7244862) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0171623) q[0];
sx q[0];
rz(-2.9201047) q[0];
sx q[0];
rz(1.4593967) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76191683) q[2];
sx q[2];
rz(-0.41315213) q[2];
sx q[2];
rz(-0.44943902) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.48737803) q[1];
sx q[1];
rz(-1.1114792) q[1];
sx q[1];
rz(0.76141255) q[1];
rz(-pi) q[2];
rz(2.7583371) q[3];
sx q[3];
rz(-2.4517422) q[3];
sx q[3];
rz(0.81724973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0059011857) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(-1.7741514) q[2];
rz(-2.588429) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(-2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32507867) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(1.1897855) q[0];
rz(-1.7143543) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(3.0292125) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6067137) q[0];
sx q[0];
rz(-1.0977854) q[0];
sx q[0];
rz(0.19434778) q[0];
rz(0.33371146) q[2];
sx q[2];
rz(-1.4881655) q[2];
sx q[2];
rz(-0.73252788) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1689414) q[1];
sx q[1];
rz(-2.2691138) q[1];
sx q[1];
rz(0.45291839) q[1];
x q[2];
rz(1.3120679) q[3];
sx q[3];
rz(-1.4530164) q[3];
sx q[3];
rz(0.044101276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(-1.7929662) q[2];
rz(-1.2049234) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(-3.1066185) q[0];
rz(-0.84683013) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(-0.91167489) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9104886) q[0];
sx q[0];
rz(-1.3470874) q[0];
sx q[0];
rz(-1.847812) q[0];
x q[1];
rz(0.20810017) q[2];
sx q[2];
rz(-2.5884183) q[2];
sx q[2];
rz(-2.2262239) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.94757838) q[1];
sx q[1];
rz(-1.7742426) q[1];
sx q[1];
rz(2.6737763) q[1];
rz(-pi) q[2];
rz(0.73260143) q[3];
sx q[3];
rz(-2.4676975) q[3];
sx q[3];
rz(-0.037308824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.70790616) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-0.24924499) q[2];
rz(2.3748659) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(0.39961091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578167) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(2.7695079) q[0];
rz(-0.58139873) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(-1.7262329) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6045195) q[0];
sx q[0];
rz(-1.3086638) q[0];
sx q[0];
rz(-1.3634691) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62016983) q[2];
sx q[2];
rz(-2.0010741) q[2];
sx q[2];
rz(-2.534453) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.806554) q[1];
sx q[1];
rz(-1.6400669) q[1];
sx q[1];
rz(-3.0663504) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0269208) q[3];
sx q[3];
rz(-1.2380935) q[3];
sx q[3];
rz(1.1247016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8156478) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(-2.9369205) q[2];
rz(1.4137911) q[3];
sx q[3];
rz(-1.1433583) q[3];
sx q[3];
rz(-2.0457101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407912) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(-1.5564556) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(-0.12116184) q[2];
sx q[2];
rz(-1.1193174) q[2];
sx q[2];
rz(0.08502273) q[2];
rz(1.4428044) q[3];
sx q[3];
rz(-2.5656869) q[3];
sx q[3];
rz(2.2681469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];