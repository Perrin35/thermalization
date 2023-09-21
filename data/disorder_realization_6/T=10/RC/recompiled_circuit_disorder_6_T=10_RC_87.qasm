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
rz(0.56086993) q[0];
rz(4.2545118) q[1];
sx q[1];
rz(1.7634044) q[1];
sx q[1];
rz(7.4982285) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44170609) q[0];
sx q[0];
rz(-2.9268648) q[0];
sx q[0];
rz(2.0798652) q[0];
x q[1];
rz(-2.4807793) q[2];
sx q[2];
rz(-1.3408957) q[2];
sx q[2];
rz(1.7125318) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.58303761) q[1];
sx q[1];
rz(-1.8999294) q[1];
sx q[1];
rz(1.0298883) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3052985) q[3];
sx q[3];
rz(-1.7146829) q[3];
sx q[3];
rz(-3.0778411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.78757301) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(-0.18307486) q[2];
rz(0.37781528) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(-0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29782444) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(-3.0644754) q[0];
rz(-0.33879694) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(-1.6024626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9248283) q[0];
sx q[0];
rz(-0.98470682) q[0];
sx q[0];
rz(-3.0490962) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23761959) q[2];
sx q[2];
rz(-2.3887206) q[2];
sx q[2];
rz(-1.876229) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96670818) q[1];
sx q[1];
rz(-1.0918573) q[1];
sx q[1];
rz(2.944988) q[1];
x q[2];
rz(1.8061403) q[3];
sx q[3];
rz(-1.9841521) q[3];
sx q[3];
rz(-1.2873161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(2.4831333) q[2];
rz(0.15130875) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(2.4466799) q[3];
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
rz(3.113134) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(-2.7084896) q[0];
rz(1.1921047) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(-2.5862397) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99291891) q[0];
sx q[0];
rz(-2.2745471) q[0];
sx q[0];
rz(0.91627319) q[0];
rz(-2.136134) q[2];
sx q[2];
rz(-1.3034504) q[2];
sx q[2];
rz(0.29758673) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1173646) q[1];
sx q[1];
rz(-1.4091361) q[1];
sx q[1];
rz(2.2092186) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8954574) q[3];
sx q[3];
rz(-0.56926308) q[3];
sx q[3];
rz(3.0666921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4042523) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(-1.8910485) q[2];
rz(0.2441497) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(-1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8811532) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(-2.326791) q[0];
rz(1.3793777) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(2.8864158) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4677306) q[0];
sx q[0];
rz(-0.47463372) q[0];
sx q[0];
rz(2.4509096) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7093676) q[2];
sx q[2];
rz(-0.6859633) q[2];
sx q[2];
rz(-1.9908817) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.081101) q[1];
sx q[1];
rz(-1.8685409) q[1];
sx q[1];
rz(2.9213419) q[1];
x q[2];
rz(-1.4820443) q[3];
sx q[3];
rz(-2.6188861) q[3];
sx q[3];
rz(2.6356217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8884376) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(-0.17318428) q[2];
rz(2.611768) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.859905) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(1.3758855) q[0];
rz(-1.8638523) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(3.0854991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.010904) q[0];
sx q[0];
rz(-2.4405257) q[0];
sx q[0];
rz(2.5581193) q[0];
rz(-1.1115083) q[2];
sx q[2];
rz(-1.2570724) q[2];
sx q[2];
rz(-2.8788061) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5406815) q[1];
sx q[1];
rz(-1.2680506) q[1];
sx q[1];
rz(1.5555698) q[1];
rz(-1.4291184) q[3];
sx q[3];
rz(-0.56483993) q[3];
sx q[3];
rz(-2.7431938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.00099480199) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(0.43219217) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11319259) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(-1.2619069) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-2.1870959) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4650824) q[0];
sx q[0];
rz(-1.0445147) q[0];
sx q[0];
rz(0.17980534) q[0];
rz(1.6184094) q[2];
sx q[2];
rz(-0.63112586) q[2];
sx q[2];
rz(2.7472251) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8442321) q[1];
sx q[1];
rz(-2.4095222) q[1];
sx q[1];
rz(-0.75433235) q[1];
rz(-0.16320634) q[3];
sx q[3];
rz(-1.8582134) q[3];
sx q[3];
rz(-2.4678469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(-1.0423638) q[2];
rz(-0.43867612) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(1.3180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2838659) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(0.74321157) q[0];
rz(1.6339533) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(2.5315703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22247032) q[0];
sx q[0];
rz(-1.3195992) q[0];
sx q[0];
rz(2.6343976) q[0];
x q[1];
rz(2.2248613) q[2];
sx q[2];
rz(-1.4787294) q[2];
sx q[2];
rz(-2.9031861) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66079084) q[1];
sx q[1];
rz(-1.5457488) q[1];
sx q[1];
rz(0.90344306) q[1];
rz(-pi) q[2];
rz(-0.20603541) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(0.68148617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(-0.91840333) q[2];
rz(-1.5911128) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.780705) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(-1.6280744) q[0];
rz(0.52945119) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(0.73658529) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37305957) q[0];
sx q[0];
rz(-0.031950843) q[0];
sx q[0];
rz(-1.2241227) q[0];
rz(-pi) q[1];
rz(-0.10891624) q[2];
sx q[2];
rz(-1.8913336) q[2];
sx q[2];
rz(0.045217302) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.96111091) q[3];
sx q[3];
rz(-2.612252) q[3];
sx q[3];
rz(-1.5101658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0344051) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(-0.68391189) q[2];
rz(-1.9125787) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(-1.3945403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.1972315) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(0.18572447) q[0];
rz(0.99705237) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(0.7448147) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72070044) q[0];
sx q[0];
rz(-1.2487131) q[0];
sx q[0];
rz(2.3548404) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1698193) q[2];
sx q[2];
rz(-0.81911659) q[2];
sx q[2];
rz(0.70105201) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77764952) q[1];
sx q[1];
rz(-1.793982) q[1];
sx q[1];
rz(1.2121483) q[1];
x q[2];
rz(2.3585988) q[3];
sx q[3];
rz(-1.8564312) q[3];
sx q[3];
rz(-2.1180958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0361438) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.194681) q[2];
rz(-0.99669325) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74334082) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(-0.031127302) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(-1.1709447) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4963213) q[0];
sx q[0];
rz(-1.5491345) q[0];
sx q[0];
rz(-2.7146656) q[0];
rz(-pi) q[1];
rz(0.77565907) q[2];
sx q[2];
rz(-1.3752898) q[2];
sx q[2];
rz(-2.069371) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78466258) q[1];
sx q[1];
rz(-1.5481879) q[1];
sx q[1];
rz(-0.008634062) q[1];
rz(-1.8627432) q[3];
sx q[3];
rz(-1.9817096) q[3];
sx q[3];
rz(1.7850072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4460454) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(0.5919624) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(0.73889417) q[2];
sx q[2];
rz(-1.0514435) q[2];
sx q[2];
rz(-2.4098868) q[2];
rz(0.016146544) q[3];
sx q[3];
rz(-1.9042249) q[3];
sx q[3];
rz(-0.92845542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
