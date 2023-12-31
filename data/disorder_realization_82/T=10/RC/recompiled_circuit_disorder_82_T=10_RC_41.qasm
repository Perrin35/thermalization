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
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(-1.27966) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0227538) q[0];
sx q[0];
rz(-1.7099329) q[0];
sx q[0];
rz(-1.3214146) q[0];
rz(1.3470634) q[2];
sx q[2];
rz(-1.1587843) q[2];
sx q[2];
rz(-1.3734693) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8259748) q[1];
sx q[1];
rz(-0.74233913) q[1];
sx q[1];
rz(-2.5151392) q[1];
rz(-pi) q[2];
rz(-1.3650595) q[3];
sx q[3];
rz(-1.7997026) q[3];
sx q[3];
rz(2.8347901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0063643) q[2];
sx q[2];
rz(-0.93868119) q[2];
sx q[2];
rz(-0.65650666) q[2];
rz(-0.73900977) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0083369) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(0.59536368) q[0];
rz(0.061925109) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(2.6541236) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90855234) q[0];
sx q[0];
rz(-1.4899583) q[0];
sx q[0];
rz(-3.0768865) q[0];
rz(-pi) q[1];
rz(1.1335282) q[2];
sx q[2];
rz(-1.184706) q[2];
sx q[2];
rz(-1.776945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.060354787) q[1];
sx q[1];
rz(-1.3831257) q[1];
sx q[1];
rz(-0.16138046) q[1];
x q[2];
rz(-1.8879714) q[3];
sx q[3];
rz(-1.4884236) q[3];
sx q[3];
rz(1.2143283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9618824) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(0.3953735) q[2];
rz(-2.0987434) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(2.9512067) q[0];
rz(-3.0186675) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(-0.22274676) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73244625) q[0];
sx q[0];
rz(-1.961686) q[0];
sx q[0];
rz(-1.976165) q[0];
x q[1];
rz(2.1026956) q[2];
sx q[2];
rz(-0.42429081) q[2];
sx q[2];
rz(-1.588856) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5388515) q[1];
sx q[1];
rz(-2.1826084) q[1];
sx q[1];
rz(1.4879423) q[1];
rz(-2.1915216) q[3];
sx q[3];
rz(-0.40682236) q[3];
sx q[3];
rz(-2.1086958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8824076) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(-1.9474691) q[2];
rz(2.0866701) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7524183) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(-1.7383204) q[0];
rz(-1.6482884) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(2.2213675) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3738149) q[0];
sx q[0];
rz(-2.2245363) q[0];
sx q[0];
rz(-0.52315229) q[0];
rz(-pi) q[1];
rz(-0.22569457) q[2];
sx q[2];
rz(-0.36964551) q[2];
sx q[2];
rz(-2.8452578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0238266) q[1];
sx q[1];
rz(-0.95925602) q[1];
sx q[1];
rz(2.141342) q[1];
x q[2];
rz(1.3654877) q[3];
sx q[3];
rz(-1.9234386) q[3];
sx q[3];
rz(2.7384788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.197864) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(1.8614004) q[2];
rz(2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(-1.357648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.458805) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(-1.7368332) q[0];
rz(2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(0.4531025) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8720854) q[0];
sx q[0];
rz(-0.082490248) q[0];
sx q[0];
rz(2.0399658) q[0];
rz(1.4278973) q[2];
sx q[2];
rz(-1.5362527) q[2];
sx q[2];
rz(2.3847716) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0799775) q[1];
sx q[1];
rz(-1.4887759) q[1];
sx q[1];
rz(-1.7916937) q[1];
rz(-pi) q[2];
rz(-2.1987678) q[3];
sx q[3];
rz(-2.4403265) q[3];
sx q[3];
rz(3.0758465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(0.72367469) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(-0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081651) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(0.23813716) q[0];
rz(-0.37995964) q[1];
sx q[1];
rz(-1.0521051) q[1];
sx q[1];
rz(1.5135117) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0673163) q[0];
sx q[0];
rz(-2.5551676) q[0];
sx q[0];
rz(-0.61074722) q[0];
rz(1.4432625) q[2];
sx q[2];
rz(-0.45293929) q[2];
sx q[2];
rz(0.73308257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7102393) q[1];
sx q[1];
rz(-2.2696886) q[1];
sx q[1];
rz(0.61183521) q[1];
rz(-pi) q[2];
rz(2.2699039) q[3];
sx q[3];
rz(-2.5786434) q[3];
sx q[3];
rz(2.9339919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0825519) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(1.9980105) q[2];
rz(3.011076) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(-0.17091621) q[3];
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
rz(2.7440199) q[0];
rz(1.3145087) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(1.9546753) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2884739) q[0];
sx q[0];
rz(-1.0736199) q[0];
sx q[0];
rz(0.32062809) q[0];
rz(0.45849623) q[2];
sx q[2];
rz(-0.065932238) q[2];
sx q[2];
rz(0.83930783) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0901511) q[1];
sx q[1];
rz(-1.4370343) q[1];
sx q[1];
rz(2.306288) q[1];
x q[2];
rz(-2.5116634) q[3];
sx q[3];
rz(-1.3327206) q[3];
sx q[3];
rz(-1.5743953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223406) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(0.85574714) q[0];
rz(3.1198655) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(-1.0303248) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4521342) q[0];
sx q[0];
rz(-1.2297213) q[0];
sx q[0];
rz(0.80490168) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.048642283) q[2];
sx q[2];
rz(-2.0970793) q[2];
sx q[2];
rz(0.7359879) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1134539) q[1];
sx q[1];
rz(-2.8061562) q[1];
sx q[1];
rz(1.0949275) q[1];
x q[2];
rz(-2.6183073) q[3];
sx q[3];
rz(-1.0417632) q[3];
sx q[3];
rz(2.7494591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8490303) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(3.0714152) q[2];
rz(0.82693806) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(-0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6475911) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(2.8253187) q[0];
rz(-2.0896185) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(2.0064328) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0907744) q[0];
sx q[0];
rz(-2.7047815) q[0];
sx q[0];
rz(0.45849053) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.029149292) q[2];
sx q[2];
rz(-1.8910732) q[2];
sx q[2];
rz(2.0707891) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3865349) q[1];
sx q[1];
rz(-0.48663501) q[1];
sx q[1];
rz(-1.6965894) q[1];
x q[2];
rz(-2.1938571) q[3];
sx q[3];
rz(-1.590168) q[3];
sx q[3];
rz(-1.3458348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6241374) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(-1.5967782) q[2];
rz(-2.4638713) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-0.28963447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90821663) q[0];
sx q[0];
rz(-2.4514618) q[0];
sx q[0];
rz(-0.35183516) q[0];
rz(0.31967638) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(-0.19616729) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0553186) q[0];
sx q[0];
rz(-1.0090758) q[0];
sx q[0];
rz(-2.4713211) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2162131) q[2];
sx q[2];
rz(-2.9973929) q[2];
sx q[2];
rz(-1.2186288) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6876467) q[1];
sx q[1];
rz(-1.5861662) q[1];
sx q[1];
rz(2.3049298) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84368002) q[3];
sx q[3];
rz(-1.0495249) q[3];
sx q[3];
rz(-2.664444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(-0.081136726) q[2];
rz(1.9598512) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(-2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022973013) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(1.8267869) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(0.20028533) q[2];
sx q[2];
rz(-2.9030048) q[2];
sx q[2];
rz(1.5422274) q[2];
rz(1.8004988) q[3];
sx q[3];
rz(-1.5256186) q[3];
sx q[3];
rz(-1.516173) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
