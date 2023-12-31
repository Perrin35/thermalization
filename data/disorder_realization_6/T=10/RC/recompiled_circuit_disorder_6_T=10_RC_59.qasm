OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(-1.7237741) q[0];
sx q[0];
rz(-0.56086993) q[0];
rz(4.2545118) q[1];
sx q[1];
rz(1.7634044) q[1];
sx q[1];
rz(7.4982285) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44170609) q[0];
sx q[0];
rz(-0.21472782) q[0];
sx q[0];
rz(-1.0617274) q[0];
rz(0.66081337) q[2];
sx q[2];
rz(-1.3408957) q[2];
sx q[2];
rz(1.7125318) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.558555) q[1];
sx q[1];
rz(-1.8999294) q[1];
sx q[1];
rz(-1.0298883) q[1];
x q[2];
rz(1.8362942) q[3];
sx q[3];
rz(-1.7146829) q[3];
sx q[3];
rz(3.0778411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.78757301) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(-2.9585178) q[2];
rz(-0.37781528) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-0.29418501) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29782444) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(-3.0644754) q[0];
rz(-2.8027957) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(1.6024626) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30277006) q[0];
sx q[0];
rz(-1.6478224) q[0];
sx q[0];
rz(-2.1588615) q[0];
rz(-pi) q[1];
rz(-1.3537172) q[2];
sx q[2];
rz(-0.84393822) q[2];
sx q[2];
rz(-0.94490563) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7669249) q[1];
sx q[1];
rz(-0.51480773) q[1];
sx q[1];
rz(-1.9306081) q[1];
rz(-1.3354524) q[3];
sx q[3];
rz(-1.9841521) q[3];
sx q[3];
rz(-1.2873161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.845528) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(2.4831333) q[2];
rz(0.15130875) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(0.69491274) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.113134) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(2.7084896) q[0];
rz(1.1921047) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(2.5862397) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0387602) q[0];
sx q[0];
rz(-2.0534678) q[0];
sx q[0];
rz(-2.32248) q[0];
rz(1.0981512) q[2];
sx q[2];
rz(-0.61908365) q[2];
sx q[2];
rz(1.4738136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8086116) q[1];
sx q[1];
rz(-0.65578991) q[1];
sx q[1];
rz(-1.8379184) q[1];
rz(1.025612) q[3];
sx q[3];
rz(-1.3979997) q[3];
sx q[3];
rz(1.9219414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.73734036) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(1.8910485) q[2];
rz(-0.2441497) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(-1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(2.326791) q[0];
rz(1.3793777) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(2.8864158) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92514738) q[0];
sx q[0];
rz(-1.9307923) q[0];
sx q[0];
rz(-1.2544592) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43222506) q[2];
sx q[2];
rz(-0.6859633) q[2];
sx q[2];
rz(1.150711) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7123588) q[1];
sx q[1];
rz(-2.7731967) q[1];
sx q[1];
rz(-0.952094) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0917986) q[3];
sx q[3];
rz(-1.61506) q[3];
sx q[3];
rz(1.1417768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8884376) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(2.9684084) q[2];
rz(-0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(-3.0392652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(-1.7657071) q[0];
rz(-1.8638523) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(0.056093562) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1306886) q[0];
sx q[0];
rz(-2.4405257) q[0];
sx q[0];
rz(-2.5581193) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1115083) q[2];
sx q[2];
rz(-1.8845203) q[2];
sx q[2];
rz(2.8788061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4896506) q[1];
sx q[1];
rz(-0.30311668) q[1];
sx q[1];
rz(-0.048708212) q[1];
x q[2];
rz(-3.0523473) q[3];
sx q[3];
rz(-2.1292994) q[3];
sx q[3];
rz(-2.9104779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1405979) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(-0.43219217) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11319259) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(0.64754852) q[0];
rz(-1.8796857) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-0.9544968) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19676767) q[0];
sx q[0];
rz(-1.4155354) q[0];
sx q[0];
rz(1.0374271) q[0];
rz(-1.6184094) q[2];
sx q[2];
rz(-0.63112586) q[2];
sx q[2];
rz(0.39436755) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8442321) q[1];
sx q[1];
rz(-2.4095222) q[1];
sx q[1];
rz(-0.75433235) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0682085) q[3];
sx q[3];
rz(-2.8121901) q[3];
sx q[3];
rz(2.9941032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(-2.0992289) q[2];
rz(0.43867612) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(-1.3180805) q[3];
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
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2838659) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(-2.3983811) q[0];
rz(-1.6339533) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(2.5315703) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2110721) q[0];
sx q[0];
rz(-2.0606344) q[0];
sx q[0];
rz(1.2852438) q[0];
x q[1];
rz(-0.91673135) q[2];
sx q[2];
rz(-1.4787294) q[2];
sx q[2];
rz(-2.9031861) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66079084) q[1];
sx q[1];
rz(-1.5958438) q[1];
sx q[1];
rz(-0.90344306) q[1];
rz(-pi) q[2];
rz(-2.0537297) q[3];
sx q[3];
rz(-0.42290877) q[3];
sx q[3];
rz(-1.203323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(-1.5504799) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(-0.38890719) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36088762) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(-1.5135182) q[0];
rz(-2.6121415) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(2.4050074) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37305957) q[0];
sx q[0];
rz(-3.1096418) q[0];
sx q[0];
rz(-1.91747) q[0];
x q[1];
rz(3.0326764) q[2];
sx q[2];
rz(-1.250259) q[2];
sx q[2];
rz(3.0963754) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8229586) q[1];
sx q[1];
rz(-1.4689323) q[1];
sx q[1];
rz(-2.676079) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96111091) q[3];
sx q[3];
rz(-2.612252) q[3];
sx q[3];
rz(1.6314268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1071876) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(1.2290139) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(-1.3945403) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443611) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(2.9558682) q[0];
rz(-2.1445403) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(-0.7448147) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4208922) q[0];
sx q[0];
rz(-1.8928796) q[0];
sx q[0];
rz(-2.3548404) q[0];
x q[1];
rz(2.7460329) q[2];
sx q[2];
rz(-2.308508) q[2];
sx q[2];
rz(-1.8849444) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4312268) q[1];
sx q[1];
rz(-1.9201628) q[1];
sx q[1];
rz(2.9037895) q[1];
x q[2];
rz(1.9635779) q[3];
sx q[3];
rz(-0.82735705) q[3];
sx q[3];
rz(-2.3208997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1054489) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(1.194681) q[2];
rz(0.99669325) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74334082) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(-3.1104654) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(-1.9706479) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02689657) q[0];
sx q[0];
rz(-0.42744246) q[0];
sx q[0];
rz(-0.052274152) q[0];
rz(1.8413713) q[2];
sx q[2];
rz(-2.3279394) q[2];
sx q[2];
rz(2.4547581) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3569301) q[1];
sx q[1];
rz(-1.5934048) q[1];
sx q[1];
rz(3.1329586) q[1];
rz(-2.714614) q[3];
sx q[3];
rz(-1.8378165) q[3];
sx q[3];
rz(-3.0468575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4460454) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(-2.5496303) q[2];
rz(0.56636089) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(-1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3175209) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(3.042165) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(0.912491) q[2];
sx q[2];
rz(-0.9463263) q[2];
sx q[2];
rz(1.8778388) q[2];
rz(1.5242143) q[3];
sx q[3];
rz(-2.8077879) q[3];
sx q[3];
rz(2.2624364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
