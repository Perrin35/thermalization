OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6157827) q[0];
sx q[0];
rz(-1.4178185) q[0];
sx q[0];
rz(-2.5807227) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(-1.2150432) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1808704) q[0];
sx q[0];
rz(-1.7579494) q[0];
sx q[0];
rz(0.10589177) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66081337) q[2];
sx q[2];
rz(-1.3408957) q[2];
sx q[2];
rz(-1.7125318) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.58303761) q[1];
sx q[1];
rz(-1.2416632) q[1];
sx q[1];
rz(-1.0298883) q[1];
x q[2];
rz(-1.0662765) q[3];
sx q[3];
rz(-2.8404232) q[3];
sx q[3];
rz(-2.1198213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3540196) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(2.9585178) q[2];
rz(-0.37781528) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(0.29418501) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29782444) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(-3.0644754) q[0];
rz(0.33879694) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(1.5391301) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2167643) q[0];
sx q[0];
rz(-2.1568858) q[0];
sx q[0];
rz(3.0490962) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73866663) q[2];
sx q[2];
rz(-1.7324442) q[2];
sx q[2];
rz(2.6612298) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51255608) q[1];
sx q[1];
rz(-1.7450383) q[1];
sx q[1];
rz(-2.0577355) q[1];
x q[2];
rz(0.48861309) q[3];
sx q[3];
rz(-2.6693137) q[3];
sx q[3];
rz(-2.3924535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(-0.65845931) q[2];
rz(2.9902839) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(2.4466799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.113134) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(0.43310305) q[0];
rz(1.1921047) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(2.5862397) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0387602) q[0];
sx q[0];
rz(-2.0534678) q[0];
sx q[0];
rz(-0.81911266) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8279282) q[2];
sx q[2];
rz(-2.1137538) q[2];
sx q[2];
rz(-2.0344337) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8086116) q[1];
sx q[1];
rz(-0.65578991) q[1];
sx q[1];
rz(-1.8379184) q[1];
rz(-1.025612) q[3];
sx q[3];
rz(-1.7435929) q[3];
sx q[3];
rz(1.9219414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4042523) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(1.8910485) q[2];
rz(-0.2441497) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(1.4499433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(0.81480169) q[0];
rz(1.3793777) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(-2.8864158) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6107451) q[0];
sx q[0];
rz(-1.8662211) q[0];
sx q[0];
rz(2.7644964) q[0];
x q[1];
rz(-1.2404664) q[2];
sx q[2];
rz(-0.9579881) q[2];
sx q[2];
rz(-1.4532879) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42923388) q[1];
sx q[1];
rz(-0.36839596) q[1];
sx q[1];
rz(2.1894987) q[1];
rz(-pi) q[2];
rz(-1.6595483) q[3];
sx q[3];
rz(-0.52270652) q[3];
sx q[3];
rz(-0.50597092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2531551) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(-0.17318428) q[2];
rz(-0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(0.1023275) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2816876) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(-1.7657071) q[0];
rz(1.2777404) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(3.0854991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027095196) q[0];
sx q[0];
rz(-1.2074911) q[0];
sx q[0];
rz(-2.5278805) q[0];
rz(0.9390097) q[2];
sx q[2];
rz(-0.54982215) q[2];
sx q[2];
rz(0.75013559) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5406815) q[1];
sx q[1];
rz(-1.2680506) q[1];
sx q[1];
rz(-1.5555698) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7124743) q[3];
sx q[3];
rz(-0.56483993) q[3];
sx q[3];
rz(-2.7431938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1405979) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(-1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(-3.0284001) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(1.2619069) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-0.9544968) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6765103) q[0];
sx q[0];
rz(-1.0445147) q[0];
sx q[0];
rz(2.9617873) q[0];
rz(-pi) q[1];
rz(-3.1068222) q[2];
sx q[2];
rz(-2.2010942) q[2];
sx q[2];
rz(2.6882753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8050025) q[1];
sx q[1];
rz(-1.0953566) q[1];
sx q[1];
rz(2.5617983) q[1];
rz(-1.0682085) q[3];
sx q[3];
rz(-2.8121901) q[3];
sx q[3];
rz(-2.9941032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.548617) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(-2.0992289) q[2];
rz(-0.43867612) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(-1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2838659) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(1.5076393) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(-0.61002237) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3726495) q[0];
sx q[0];
rz(-0.56108755) q[0];
sx q[0];
rz(-2.6555496) q[0];
x q[1];
rz(2.2248613) q[2];
sx q[2];
rz(-1.4787294) q[2];
sx q[2];
rz(-2.9031861) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.89027379) q[1];
sx q[1];
rz(-2.2379025) q[1];
sx q[1];
rz(3.1097079) q[1];
rz(-2.0537297) q[3];
sx q[3];
rz(-0.42290877) q[3];
sx q[3];
rz(-1.203323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(-1.5504799) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.780705) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(1.5135182) q[0];
rz(2.6121415) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(-0.73658529) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7685331) q[0];
sx q[0];
rz(-3.1096418) q[0];
sx q[0];
rz(-1.2241227) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2484776) q[2];
sx q[2];
rz(-1.467448) q[2];
sx q[2];
rz(1.581574) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8229586) q[1];
sx q[1];
rz(-1.4689323) q[1];
sx q[1];
rz(2.676079) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1804817) q[3];
sx q[3];
rz(-0.52934066) q[3];
sx q[3];
rz(1.5101658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0344051) q[2];
sx q[2];
rz(-1.9341058) q[2];
sx q[2];
rz(-0.68391189) q[2];
rz(-1.2290139) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(-1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.1972315) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(2.9558682) q[0];
rz(-2.1445403) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(2.396778) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1574402) q[0];
sx q[0];
rz(-2.3072349) q[0];
sx q[0];
rz(-2.0122583) q[0];
x q[1];
rz(2.7460329) q[2];
sx q[2];
rz(-2.308508) q[2];
sx q[2];
rz(-1.8849444) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8150755) q[1];
sx q[1];
rz(-0.41985598) q[1];
sx q[1];
rz(-2.144787) q[1];
rz(-pi) q[2];
rz(0.78299384) q[3];
sx q[3];
rz(-1.2851614) q[3];
sx q[3];
rz(-2.1180958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0361438) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(1.194681) q[2];
rz(-0.99669325) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(-0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3982518) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(0.40400305) q[0];
rz(0.031127302) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(1.9706479) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4963213) q[0];
sx q[0];
rz(-1.5491345) q[0];
sx q[0];
rz(0.42692703) q[0];
rz(-pi) q[1];
rz(-0.77565907) q[2];
sx q[2];
rz(-1.3752898) q[2];
sx q[2];
rz(2.069371) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78466258) q[1];
sx q[1];
rz(-1.5934048) q[1];
sx q[1];
rz(-0.008634062) q[1];
rz(-pi) q[2];
rz(-0.42697866) q[3];
sx q[3];
rz(-1.3037762) q[3];
sx q[3];
rz(0.094735183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4460454) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(0.5919624) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(-1.5238354) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3175209) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(3.042165) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(-2.4026985) q[2];
sx q[2];
rz(-1.0514435) q[2];
sx q[2];
rz(-2.4098868) q[2];
rz(1.6173784) q[3];
sx q[3];
rz(-0.33380476) q[3];
sx q[3];
rz(-0.87915626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
