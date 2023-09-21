OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3368971) q[0];
sx q[0];
rz(-2.1043632) q[0];
sx q[0];
rz(-0.35559911) q[0];
rz(2.8388677) q[1];
sx q[1];
rz(-1.0441138) q[1];
sx q[1];
rz(1.27966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0227538) q[0];
sx q[0];
rz(-1.4316598) q[0];
sx q[0];
rz(-1.3214146) q[0];
rz(-pi) q[1];
rz(-1.3470634) q[2];
sx q[2];
rz(-1.1587843) q[2];
sx q[2];
rz(-1.7681233) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.31561786) q[1];
sx q[1];
rz(-0.74233913) q[1];
sx q[1];
rz(-2.5151392) q[1];
x q[2];
rz(-1.3650595) q[3];
sx q[3];
rz(-1.7997026) q[3];
sx q[3];
rz(-0.3068026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0063643) q[2];
sx q[2];
rz(-0.93868119) q[2];
sx q[2];
rz(-0.65650666) q[2];
rz(-0.73900977) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083369) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(2.546229) q[0];
rz(-3.0796675) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(2.6541236) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66747626) q[0];
sx q[0];
rz(-1.6352909) q[0];
sx q[0];
rz(-1.6518031) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1335282) q[2];
sx q[2];
rz(-1.9568866) q[2];
sx q[2];
rz(-1.776945) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6007874) q[1];
sx q[1];
rz(-1.4122737) q[1];
sx q[1];
rz(1.3807138) q[1];
x q[2];
rz(-1.8879714) q[3];
sx q[3];
rz(-1.653169) q[3];
sx q[3];
rz(1.9272643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9618824) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(-0.3953735) q[2];
rz(2.0987434) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9449126) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(2.9512067) q[0];
rz(-3.0186675) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(-2.9188459) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5771772) q[0];
sx q[0];
rz(-0.555493) q[0];
sx q[0];
rz(0.76340686) q[0];
x q[1];
rz(-1.0388971) q[2];
sx q[2];
rz(-2.7173018) q[2];
sx q[2];
rz(1.588856) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.682444) q[1];
sx q[1];
rz(-0.61668452) q[1];
sx q[1];
rz(-0.1174121) q[1];
x q[2];
rz(-1.907903) q[3];
sx q[3];
rz(-1.3385696) q[3];
sx q[3];
rz(-3.0984578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8824076) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(1.1941236) q[2];
rz(-1.0549226) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(-0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38917437) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(1.4032723) q[0];
rz(-1.6482884) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(-0.9202252) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46566761) q[0];
sx q[0];
rz(-1.9786069) q[0];
sx q[0];
rz(-0.84665914) q[0];
x q[1];
rz(-2.9158981) q[2];
sx q[2];
rz(-2.7719471) q[2];
sx q[2];
rz(-2.8452578) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8060382) q[1];
sx q[1];
rz(-1.1127377) q[1];
sx q[1];
rz(0.6946509) q[1];
rz(-pi) q[2];
x q[2];
rz(2.782015) q[3];
sx q[3];
rz(-1.3782856) q[3];
sx q[3];
rz(1.0958835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(1.2801923) q[2];
rz(-0.81930339) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6827877) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(1.4047594) q[0];
rz(2.3732896) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(-0.4531025) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4015409) q[0];
sx q[0];
rz(-1.6443559) q[0];
sx q[0];
rz(3.1042276) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7136953) q[2];
sx q[2];
rz(-1.5362527) q[2];
sx q[2];
rz(0.75682109) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0616152) q[1];
sx q[1];
rz(-1.4887759) q[1];
sx q[1];
rz(-1.3498989) q[1];
rz(-0.94282486) q[3];
sx q[3];
rz(-2.4403265) q[3];
sx q[3];
rz(-3.0758465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0129464) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033427514) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(-0.23813716) q[0];
rz(-0.37995964) q[1];
sx q[1];
rz(-1.0521051) q[1];
sx q[1];
rz(-1.628081) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0673163) q[0];
sx q[0];
rz(-2.5551676) q[0];
sx q[0];
rz(-2.5308454) q[0];
rz(0.061821826) q[2];
sx q[2];
rz(-1.1218058) q[2];
sx q[2];
rz(0.87473727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.39735079) q[1];
sx q[1];
rz(-2.2480818) q[1];
sx q[1];
rz(0.97126295) q[1];
rz(-0.38576491) q[3];
sx q[3];
rz(-1.1499975) q[3];
sx q[3];
rz(2.5667131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(-1.1435821) q[2];
rz(0.13051662) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1259574) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(-0.39757279) q[0];
rz(-1.3145087) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(-1.1869173) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2884739) q[0];
sx q[0];
rz(-2.0679727) q[0];
sx q[0];
rz(-0.32062809) q[0];
x q[1];
rz(1.5415807) q[2];
sx q[2];
rz(-1.511682) q[2];
sx q[2];
rz(-2.7616449) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3992607) q[1];
sx q[1];
rz(-0.84335828) q[1];
sx q[1];
rz(-2.9620693) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5116634) q[3];
sx q[3];
rz(-1.3327206) q[3];
sx q[3];
rz(-1.5743953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26178965) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(-1.7857893) q[2];
rz(1.5073744) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81925201) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(0.85574714) q[0];
rz(-0.021727173) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(2.1112679) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5949769) q[0];
sx q[0];
rz(-0.82406509) q[0];
sx q[0];
rz(-2.0440408) q[0];
rz(-pi) q[1];
rz(1.0439992) q[2];
sx q[2];
rz(-1.5287405) q[2];
sx q[2];
rz(2.2823357) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1460261) q[1];
sx q[1];
rz(-1.4194173) q[1];
sx q[1];
rz(1.2703018) q[1];
rz(2.6183073) q[3];
sx q[3];
rz(-2.0998294) q[3];
sx q[3];
rz(2.7494591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8490303) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(-0.070177468) q[2];
rz(-2.3146546) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(-0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6475911) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-2.8253187) q[0];
rz(-1.0519741) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(-2.0064328) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2410256) q[0];
sx q[0];
rz(-1.3824438) q[0];
sx q[0];
rz(-0.39649773) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4831545) q[2];
sx q[2];
rz(-0.32155514) q[2];
sx q[2];
rz(-1.9784387) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70442048) q[1];
sx q[1];
rz(-1.6295027) q[1];
sx q[1];
rz(-1.0874332) q[1];
rz(1.6039861) q[3];
sx q[3];
rz(-2.5182708) q[3];
sx q[3];
rz(2.8896796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(1.5967782) q[2];
rz(-0.67772135) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.233376) q[0];
sx q[0];
rz(-2.4514618) q[0];
sx q[0];
rz(0.35183516) q[0];
rz(2.8219163) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(0.19616729) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10712121) q[0];
sx q[0];
rz(-2.2959318) q[0];
sx q[0];
rz(-0.79191533) q[0];
rz(1.2162131) q[2];
sx q[2];
rz(-0.14419975) q[2];
sx q[2];
rz(-1.2186288) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6876467) q[1];
sx q[1];
rz(-1.5554264) q[1];
sx q[1];
rz(0.83666283) q[1];
rz(2.2833061) q[3];
sx q[3];
rz(-0.86601102) q[3];
sx q[3];
rz(-1.5376877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7508042) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(0.081136726) q[2];
rz(-1.9598512) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(1.0749764) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022973013) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(-1.3148057) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(2.9076004) q[2];
sx q[2];
rz(-1.5237612) q[2];
sx q[2];
rz(0.16618726) q[2];
rz(-0.046394596) q[3];
sx q[3];
rz(-1.3413324) q[3];
sx q[3];
rz(-3.0975292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];