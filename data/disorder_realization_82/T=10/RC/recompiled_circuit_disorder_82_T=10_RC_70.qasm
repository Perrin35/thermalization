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
rz(4.1788221) q[0];
sx q[0];
rz(9.0691789) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(-1.27966) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1188388) q[0];
sx q[0];
rz(-1.7099329) q[0];
sx q[0];
rz(1.820178) q[0];
rz(-pi) q[1];
rz(2.7202397) q[2];
sx q[2];
rz(-1.7755277) q[2];
sx q[2];
rz(0.28819627) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31561786) q[1];
sx q[1];
rz(-2.3992535) q[1];
sx q[1];
rz(2.5151392) q[1];
rz(-pi) q[2];
rz(1.7765331) q[3];
sx q[3];
rz(-1.3418901) q[3];
sx q[3];
rz(-2.8347901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1352284) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(2.485086) q[2];
rz(0.73900977) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083369) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(2.546229) q[0];
rz(3.0796675) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(-0.48746902) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4741164) q[0];
sx q[0];
rz(-1.5063018) q[0];
sx q[0];
rz(1.6518031) q[0];
rz(0.42178085) q[2];
sx q[2];
rz(-1.9739208) q[2];
sx q[2];
rz(-3.1096854) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5408052) q[1];
sx q[1];
rz(-1.4122737) q[1];
sx q[1];
rz(1.3807138) q[1];
rz(-pi) q[2];
rz(1.3120193) q[3];
sx q[3];
rz(-2.8142455) q[3];
sx q[3];
rz(0.60206383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9618824) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(-0.3953735) q[2];
rz(-1.0428492) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(-0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(0.19038598) q[0];
rz(-3.0186675) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(0.22274676) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5771772) q[0];
sx q[0];
rz(-2.5860997) q[0];
sx q[0];
rz(0.76340686) q[0];
x q[1];
rz(2.9163755) q[2];
sx q[2];
rz(-1.9334963) q[2];
sx q[2];
rz(-0.97937102) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0157156) q[1];
sx q[1];
rz(-1.6385957) q[1];
sx q[1];
rz(-2.5281639) q[1];
rz(2.1915216) q[3];
sx q[3];
rz(-0.40682236) q[3];
sx q[3];
rz(-1.0328968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8824076) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(1.9474691) q[2];
rz(1.0549226) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(-0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7524183) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(1.4032723) q[0];
rz(1.6482884) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(-2.2213675) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46566761) q[0];
sx q[0];
rz(-1.9786069) q[0];
sx q[0];
rz(-2.2949335) q[0];
rz(2.9158981) q[2];
sx q[2];
rz(-2.7719471) q[2];
sx q[2];
rz(-0.29633488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0238266) q[1];
sx q[1];
rz(-2.1823366) q[1];
sx q[1];
rz(-2.141342) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3595777) q[3];
sx q[3];
rz(-1.7633071) q[3];
sx q[3];
rz(-1.0958835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.197864) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(-1.8614004) q[2];
rz(-2.3222893) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6827877) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(-1.7368332) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(2.6884902) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3081007) q[0];
sx q[0];
rz(-1.6080603) q[0];
sx q[0];
rz(-1.4971855) q[0];
rz(-0.034899072) q[2];
sx q[2];
rz(-1.7136095) q[2];
sx q[2];
rz(-2.3226483) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52757712) q[1];
sx q[1];
rz(-1.7909389) q[1];
sx q[1];
rz(0.084053587) q[1];
rz(-pi) q[2];
rz(-0.94282486) q[3];
sx q[3];
rz(-2.4403265) q[3];
sx q[3];
rz(0.065746106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0129464) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(1.6298693) q[2];
rz(2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033427514) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(-2.9034555) q[0];
rz(0.37995964) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(1.5135117) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0742764) q[0];
sx q[0];
rz(-2.5551676) q[0];
sx q[0];
rz(-0.61074722) q[0];
x q[1];
rz(-1.1210576) q[2];
sx q[2];
rz(-1.5151086) q[2];
sx q[2];
rz(-2.4186717) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43135333) q[1];
sx q[1];
rz(-0.87190404) q[1];
sx q[1];
rz(-0.61183521) q[1];
rz(-pi) q[2];
rz(0.38576491) q[3];
sx q[3];
rz(-1.1499975) q[3];
sx q[3];
rz(-2.5667131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0825519) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(-1.9980105) q[2];
rz(3.011076) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(0.17091621) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1259574) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(2.7440199) q[0];
rz(-1.3145087) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(-1.1869173) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2448954) q[0];
sx q[0];
rz(-2.5573686) q[0];
sx q[0];
rz(-1.0446192) q[0];
x q[1];
rz(2.6830964) q[2];
sx q[2];
rz(-0.065932238) q[2];
sx q[2];
rz(2.3022848) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.051441593) q[1];
sx q[1];
rz(-1.4370343) q[1];
sx q[1];
rz(2.306288) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39077057) q[3];
sx q[3];
rz(-2.4739389) q[3];
sx q[3];
rz(2.8323176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.879803) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(-1.3558033) q[2];
rz(-1.5073744) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(-0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223406) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(0.85574714) q[0];
rz(3.1198655) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(2.1112679) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4521342) q[0];
sx q[0];
rz(-1.9118714) q[0];
sx q[0];
rz(-0.80490168) q[0];
rz(-1.0439992) q[2];
sx q[2];
rz(-1.6128522) q[2];
sx q[2];
rz(-0.85925697) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9955666) q[1];
sx q[1];
rz(-1.4194173) q[1];
sx q[1];
rz(-1.8712908) q[1];
rz(-pi) q[2];
rz(2.1645032) q[3];
sx q[3];
rz(-2.0168552) q[3];
sx q[3];
rz(0.89531985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.29256233) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(0.070177468) q[2];
rz(0.82693806) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4940015) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(0.31627396) q[0];
rz(-2.0896185) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(-1.1351599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59199698) q[0];
sx q[0];
rz(-1.1816918) q[0];
sx q[0];
rz(-1.7745716) q[0];
rz(3.1124434) q[2];
sx q[2];
rz(-1.8910732) q[2];
sx q[2];
rz(2.0707891) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.70442048) q[1];
sx q[1];
rz(-1.51209) q[1];
sx q[1];
rz(2.0541595) q[1];
rz(0.9477356) q[3];
sx q[3];
rz(-1.5514246) q[3];
sx q[3];
rz(-1.7957578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51745522) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(-1.5967782) q[2];
rz(-0.67772135) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.233376) q[0];
sx q[0];
rz(-2.4514618) q[0];
sx q[0];
rz(-2.7897575) q[0];
rz(-2.8219163) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(0.19616729) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10712121) q[0];
sx q[0];
rz(-0.84566085) q[0];
sx q[0];
rz(0.79191533) q[0];
rz(-pi) q[1];
rz(1.9253795) q[2];
sx q[2];
rz(-2.9973929) q[2];
sx q[2];
rz(-1.2186288) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0386104) q[1];
sx q[1];
rz(-0.83676941) q[1];
sx q[1];
rz(-3.1208913) q[1];
x q[2];
rz(0.85828652) q[3];
sx q[3];
rz(-2.2755816) q[3];
sx q[3];
rz(1.603905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3907884) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(-0.081136726) q[2];
rz(1.1817415) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(2.0666163) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.5224456) q[2];
sx q[2];
rz(-1.3370677) q[2];
sx q[2];
rz(1.7481902) q[2];
rz(1.374791) q[3];
sx q[3];
rz(-2.9075665) q[3];
sx q[3];
rz(-2.8961765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
