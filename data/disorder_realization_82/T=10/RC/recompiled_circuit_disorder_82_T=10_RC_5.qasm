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
rz(1.8619327) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049392603) q[0];
sx q[0];
rz(-0.28486262) q[0];
sx q[0];
rz(-2.0869135) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7202397) q[2];
sx q[2];
rz(-1.7755277) q[2];
sx q[2];
rz(2.8533964) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0919839) q[1];
sx q[1];
rz(-0.99124747) q[1];
sx q[1];
rz(1.0773354) q[1];
rz(-pi) q[2];
rz(2.9079307) q[3];
sx q[3];
rz(-1.7710925) q[3];
sx q[3];
rz(1.3113126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1352284) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(-0.65650666) q[2];
rz(-0.73900977) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(-2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1332557) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(0.59536368) q[0];
rz(-0.061925109) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(2.6541236) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4741164) q[0];
sx q[0];
rz(-1.5063018) q[0];
sx q[0];
rz(1.4897896) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0080645) q[2];
sx q[2];
rz(-1.184706) q[2];
sx q[2];
rz(-1.3646477) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6007874) q[1];
sx q[1];
rz(-1.729319) q[1];
sx q[1];
rz(-1.7608789) q[1];
rz(1.2536212) q[3];
sx q[3];
rz(-1.4884236) q[3];
sx q[3];
rz(-1.9272643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9618824) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(-2.7462192) q[2];
rz(1.0428492) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9449126) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(-2.9512067) q[0];
rz(-3.0186675) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(0.22274676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73244625) q[0];
sx q[0];
rz(-1.1799066) q[0];
sx q[0];
rz(1.976165) q[0];
rz(-pi) q[1];
rz(1.1995302) q[2];
sx q[2];
rz(-1.3604593) q[2];
sx q[2];
rz(-2.631275) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0157156) q[1];
sx q[1];
rz(-1.6385957) q[1];
sx q[1];
rz(2.5281639) q[1];
rz(-pi) q[2];
rz(1.2336897) q[3];
sx q[3];
rz(-1.8030231) q[3];
sx q[3];
rz(3.0984578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8824076) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(-1.9474691) q[2];
rz(-2.0866701) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(-2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38917437) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(1.7383204) q[0];
rz(-1.4933043) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(-0.9202252) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.675925) q[0];
sx q[0];
rz(-1.9786069) q[0];
sx q[0];
rz(-2.2949335) q[0];
x q[1];
rz(0.22569457) q[2];
sx q[2];
rz(-0.36964551) q[2];
sx q[2];
rz(-0.29633488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3355545) q[1];
sx q[1];
rz(-2.028855) q[1];
sx q[1];
rz(-2.4469417) q[1];
rz(-0.50587378) q[3];
sx q[3];
rz(-0.40588356) q[3];
sx q[3];
rz(-0.94569262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.197864) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(-1.8614004) q[2];
rz(2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(-1.357648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.458805) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(2.3732896) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(-0.4531025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74005175) q[0];
sx q[0];
rz(-1.4972367) q[0];
sx q[0];
rz(-0.037365035) q[0];
rz(3.1066936) q[2];
sx q[2];
rz(-1.7136095) q[2];
sx q[2];
rz(0.81894433) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15935005) q[1];
sx q[1];
rz(-2.9061926) q[1];
sx q[1];
rz(1.2118641) q[1];
rz(-pi) q[2];
rz(2.1702607) q[3];
sx q[3];
rz(-1.9595651) q[3];
sx q[3];
rz(2.142981) q[3];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1081651) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(-2.9034555) q[0];
rz(2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.5135117) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0742764) q[0];
sx q[0];
rz(-2.5551676) q[0];
sx q[0];
rz(-2.5308454) q[0];
x q[1];
rz(-2.020535) q[2];
sx q[2];
rz(-1.5151086) q[2];
sx q[2];
rz(2.4186717) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7102393) q[1];
sx q[1];
rz(-2.2696886) q[1];
sx q[1];
rz(-2.5297574) q[1];
rz(-0.87168872) q[3];
sx q[3];
rz(-2.5786434) q[3];
sx q[3];
rz(2.9339919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(1.1435821) q[2];
rz(-0.13051662) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(-0.17091621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-0.54388261) q[1];
sx q[1];
rz(-1.9546753) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884739) q[0];
sx q[0];
rz(-2.0679727) q[0];
sx q[0];
rz(-2.8209646) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.600012) q[2];
sx q[2];
rz(-1.6299106) q[2];
sx q[2];
rz(2.7616449) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4759051) q[1];
sx q[1];
rz(-0.74531065) q[1];
sx q[1];
rz(-1.7687294) q[1];
x q[2];
rz(1.8625453) q[3];
sx q[3];
rz(-0.96127931) q[3];
sx q[3];
rz(0.17385829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.26178965) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(-1.3558033) q[2];
rz(-1.5073744) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223406) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(0.85574714) q[0];
rz(-0.021727173) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(-2.1112679) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5949769) q[0];
sx q[0];
rz(-2.3175276) q[0];
sx q[0];
rz(-2.0440408) q[0];
rz(-0.048642283) q[2];
sx q[2];
rz(-2.0970793) q[2];
sx q[2];
rz(0.7359879) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1134539) q[1];
sx q[1];
rz(-0.33543643) q[1];
sx q[1];
rz(2.0466652) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2780667) q[3];
sx q[3];
rz(-0.72609767) q[3];
sx q[3];
rz(1.2442148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8490303) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(3.0714152) q[2];
rz(2.3146546) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(-0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6475911) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(-0.31627396) q[0];
rz(-2.0896185) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(-1.1351599) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0508182) q[0];
sx q[0];
rz(-2.7047815) q[0];
sx q[0];
rz(-2.6831021) q[0];
x q[1];
rz(-1.6584381) q[2];
sx q[2];
rz(-2.8200375) q[2];
sx q[2];
rz(-1.9784387) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4371722) q[1];
sx q[1];
rz(-1.6295027) q[1];
sx q[1];
rz(-2.0541595) q[1];
x q[2];
rz(1.5376066) q[3];
sx q[3];
rz(-2.5182708) q[3];
sx q[3];
rz(0.25191307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(-1.5967782) q[2];
rz(0.67772135) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-0.28963447) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.233376) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(-2.7897575) q[0];
rz(0.31967638) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(0.19616729) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0553186) q[0];
sx q[0];
rz(-2.1325169) q[0];
sx q[0];
rz(-0.67027153) q[0];
x q[1];
rz(-3.0912193) q[2];
sx q[2];
rz(-1.435624) q[2];
sx q[2];
rz(-0.86063517) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0077121) q[1];
sx q[1];
rz(-0.73426437) q[1];
sx q[1];
rz(1.5937362) q[1];
rz(2.2833061) q[3];
sx q[3];
rz(-2.2755816) q[3];
sx q[3];
rz(1.5376877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(0.081136726) q[2];
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
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.022973013) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(-1.3148057) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(1.5224456) q[2];
sx q[2];
rz(-1.804525) q[2];
sx q[2];
rz(-1.3934025) q[2];
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