OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7689826) q[0];
sx q[0];
rz(-3.09642) q[0];
sx q[0];
rz(0.67155182) q[0];
rz(2.1454732) q[1];
sx q[1];
rz(5.4944333) q[1];
sx q[1];
rz(6.5715437) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.45123) q[0];
sx q[0];
rz(-1.3730064) q[0];
sx q[0];
rz(-3.1044699) q[0];
rz(-2.4774083) q[2];
sx q[2];
rz(-1.3553047) q[2];
sx q[2];
rz(1.8813949) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3934196) q[1];
sx q[1];
rz(-0.44965023) q[1];
sx q[1];
rz(3.054674) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45436556) q[3];
sx q[3];
rz(-1.7011257) q[3];
sx q[3];
rz(2.5443175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4396189) q[2];
sx q[2];
rz(-2.7957323) q[2];
sx q[2];
rz(0.71887476) q[2];
rz(1.6950722) q[3];
sx q[3];
rz(-1.6547763) q[3];
sx q[3];
rz(2.1046624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.059939) q[0];
sx q[0];
rz(-0.43123284) q[0];
sx q[0];
rz(0.28847873) q[0];
rz(-0.55229315) q[1];
sx q[1];
rz(-2.0497132) q[1];
sx q[1];
rz(-1.1757895) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4058454) q[0];
sx q[0];
rz(-1.7205392) q[0];
sx q[0];
rz(1.2558054) q[0];
rz(-pi) q[1];
rz(1.3973049) q[2];
sx q[2];
rz(-1.6855006) q[2];
sx q[2];
rz(2.3115273) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.77681357) q[1];
sx q[1];
rz(-0.59092593) q[1];
sx q[1];
rz(-1.7950115) q[1];
rz(-pi) q[2];
rz(1.5351686) q[3];
sx q[3];
rz(-3.1094915) q[3];
sx q[3];
rz(-0.42313448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6931849) q[2];
sx q[2];
rz(-1.2625445) q[2];
sx q[2];
rz(0.022424879) q[2];
rz(-1.7153995) q[3];
sx q[3];
rz(-1.8636999) q[3];
sx q[3];
rz(2.2414331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4876323) q[0];
sx q[0];
rz(-1.910169) q[0];
sx q[0];
rz(0.76474977) q[0];
rz(-1.359831) q[1];
sx q[1];
rz(-1.9197074) q[1];
sx q[1];
rz(1.2545895) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6343289) q[0];
sx q[0];
rz(-1.6741236) q[0];
sx q[0];
rz(1.4402585) q[0];
x q[1];
rz(2.7255035) q[2];
sx q[2];
rz(-0.35230428) q[2];
sx q[2];
rz(1.2809629) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4855328) q[1];
sx q[1];
rz(-1.1120218) q[1];
sx q[1];
rz(1.8937673) q[1];
x q[2];
rz(-1.6996918) q[3];
sx q[3];
rz(-1.9106955) q[3];
sx q[3];
rz(-1.5101907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4059056) q[2];
sx q[2];
rz(-2.3730706) q[2];
sx q[2];
rz(-3.1206257) q[2];
rz(-1.5812801) q[3];
sx q[3];
rz(-1.0617278) q[3];
sx q[3];
rz(0.90644065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86879325) q[0];
sx q[0];
rz(-0.11196207) q[0];
sx q[0];
rz(-1.772076) q[0];
rz(0.70961332) q[1];
sx q[1];
rz(-1.2779002) q[1];
sx q[1];
rz(-0.69724625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5013789) q[0];
sx q[0];
rz(-1.1637299) q[0];
sx q[0];
rz(-1.908692) q[0];
rz(2.286381) q[2];
sx q[2];
rz(-0.89323275) q[2];
sx q[2];
rz(-1.8032359) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.62909758) q[1];
sx q[1];
rz(-1.8726908) q[1];
sx q[1];
rz(-0.90934335) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3342821) q[3];
sx q[3];
rz(-1.7380889) q[3];
sx q[3];
rz(3.0574556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9559481) q[2];
sx q[2];
rz(-2.1944025) q[2];
sx q[2];
rz(-2.4626125) q[2];
rz(1.125157) q[3];
sx q[3];
rz(-1.9966639) q[3];
sx q[3];
rz(0.2230491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0515902) q[0];
sx q[0];
rz(-1.2035878) q[0];
sx q[0];
rz(-3.0644655) q[0];
rz(-1.1513101) q[1];
sx q[1];
rz(-1.5682033) q[1];
sx q[1];
rz(-0.77879771) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0317694) q[0];
sx q[0];
rz(-1.6187877) q[0];
sx q[0];
rz(-0.035295156) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9735317) q[2];
sx q[2];
rz(-1.2894783) q[2];
sx q[2];
rz(0.14708731) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.59757876) q[1];
sx q[1];
rz(-0.42002267) q[1];
sx q[1];
rz(0.93552621) q[1];
x q[2];
rz(0.16186951) q[3];
sx q[3];
rz(-1.5613982) q[3];
sx q[3];
rz(0.19111253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.516958) q[2];
sx q[2];
rz(-1.4411074) q[2];
sx q[2];
rz(-1.9514294) q[2];
rz(-0.55365753) q[3];
sx q[3];
rz(-0.62479574) q[3];
sx q[3];
rz(0.73738086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5364285) q[0];
sx q[0];
rz(-1.9909415) q[0];
sx q[0];
rz(0.27106699) q[0];
rz(-1.0003264) q[1];
sx q[1];
rz(-1.9258291) q[1];
sx q[1];
rz(-2.2727374) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0603179) q[0];
sx q[0];
rz(-0.79303128) q[0];
sx q[0];
rz(-0.6834553) q[0];
rz(3.0539114) q[2];
sx q[2];
rz(-2.4041135) q[2];
sx q[2];
rz(-0.35754851) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.083243283) q[1];
sx q[1];
rz(-1.740137) q[1];
sx q[1];
rz(1.9237299) q[1];
rz(-pi) q[2];
rz(-0.47021659) q[3];
sx q[3];
rz(-2.3002671) q[3];
sx q[3];
rz(0.39455345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.038593682) q[2];
sx q[2];
rz(-2.8391892) q[2];
sx q[2];
rz(1.137286) q[2];
rz(0.31050995) q[3];
sx q[3];
rz(-0.91028428) q[3];
sx q[3];
rz(0.92946068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.0248658) q[0];
sx q[0];
rz(-1.7127345) q[0];
sx q[0];
rz(-1.0799991) q[0];
rz(0.96177167) q[1];
sx q[1];
rz(-2.5764143) q[1];
sx q[1];
rz(-0.22294179) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59793845) q[0];
sx q[0];
rz(-0.035654457) q[0];
sx q[0];
rz(1.7810506) q[0];
rz(-2.8330363) q[2];
sx q[2];
rz(-0.94538222) q[2];
sx q[2];
rz(-0.61185123) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0698439) q[1];
sx q[1];
rz(-0.78673601) q[1];
sx q[1];
rz(1.9658425) q[1];
rz(-pi) q[2];
rz(1.5490398) q[3];
sx q[3];
rz(-0.59594107) q[3];
sx q[3];
rz(-0.21778743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96233931) q[2];
sx q[2];
rz(-2.8620359) q[2];
sx q[2];
rz(-1.1780098) q[2];
rz(-0.32633215) q[3];
sx q[3];
rz(-0.81063619) q[3];
sx q[3];
rz(-2.0012205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2278263) q[0];
sx q[0];
rz(-2.1837809) q[0];
sx q[0];
rz(0.80818278) q[0];
rz(-2.783964) q[1];
sx q[1];
rz(-1.459815) q[1];
sx q[1];
rz(-1.0825895) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4935166) q[0];
sx q[0];
rz(-0.89663038) q[0];
sx q[0];
rz(-2.9229259) q[0];
rz(-0.2530667) q[2];
sx q[2];
rz(-2.1089206) q[2];
sx q[2];
rz(-2.5122364) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8053003) q[1];
sx q[1];
rz(-2.3545579) q[1];
sx q[1];
rz(-1.8563104) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6851366) q[3];
sx q[3];
rz(-1.136406) q[3];
sx q[3];
rz(1.020177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.064934405) q[2];
sx q[2];
rz(-2.6532463) q[2];
sx q[2];
rz(1.2154382) q[2];
rz(2.2293034) q[3];
sx q[3];
rz(-2.2513697) q[3];
sx q[3];
rz(2.598855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2927581) q[0];
sx q[0];
rz(-0.64105761) q[0];
sx q[0];
rz(-2.7178398) q[0];
rz(-1.1774225) q[1];
sx q[1];
rz(-2.2260428) q[1];
sx q[1];
rz(-0.37568572) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7862512) q[0];
sx q[0];
rz(-1.4999985) q[0];
sx q[0];
rz(-1.7485445) q[0];
x q[1];
rz(-1.5440953) q[2];
sx q[2];
rz(-1.809568) q[2];
sx q[2];
rz(-1.2899866) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60236193) q[1];
sx q[1];
rz(-1.5054387) q[1];
sx q[1];
rz(-1.0687625) q[1];
rz(-pi) q[2];
rz(0.71640941) q[3];
sx q[3];
rz(-1.3103974) q[3];
sx q[3];
rz(1.474787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0287013) q[2];
sx q[2];
rz(-1.4914923) q[2];
sx q[2];
rz(-2.4493307) q[2];
rz(-0.68474692) q[3];
sx q[3];
rz(-2.1919577) q[3];
sx q[3];
rz(-0.14028604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.5570062) q[0];
sx q[0];
rz(-0.96730119) q[0];
sx q[0];
rz(-1.517357) q[0];
rz(1.5380305) q[1];
sx q[1];
rz(-1.3471194) q[1];
sx q[1];
rz(-1.921152) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9346817) q[0];
sx q[0];
rz(-2.5162272) q[0];
sx q[0];
rz(1.6726794) q[0];
rz(-2.6399355) q[2];
sx q[2];
rz(-1.3634472) q[2];
sx q[2];
rz(0.32688552) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8216202) q[1];
sx q[1];
rz(-1.8516685) q[1];
sx q[1];
rz(-0.73001659) q[1];
rz(1.8655928) q[3];
sx q[3];
rz(-2.285897) q[3];
sx q[3];
rz(-3.0218647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5753691) q[2];
sx q[2];
rz(-1.0074002) q[2];
sx q[2];
rz(-2.2929906) q[2];
rz(0.97950116) q[3];
sx q[3];
rz(-1.5666015) q[3];
sx q[3];
rz(0.037467329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.033087) q[0];
sx q[0];
rz(-1.6642234) q[0];
sx q[0];
rz(-2.4432175) q[0];
rz(-2.5820844) q[1];
sx q[1];
rz(-2.5026176) q[1];
sx q[1];
rz(2.7284596) q[1];
rz(-0.32158659) q[2];
sx q[2];
rz(-1.814331) q[2];
sx q[2];
rz(-0.41650256) q[2];
rz(-0.091690334) q[3];
sx q[3];
rz(-2.5602362) q[3];
sx q[3];
rz(0.14701281) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
