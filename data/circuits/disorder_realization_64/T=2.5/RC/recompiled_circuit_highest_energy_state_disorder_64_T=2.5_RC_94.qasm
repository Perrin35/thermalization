OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.1177899) q[0];
sx q[0];
rz(-1.1060214) q[0];
sx q[0];
rz(2.3646781) q[0];
rz(1.39224) q[1];
sx q[1];
rz(-1.3148146) q[1];
sx q[1];
rz(-0.97631747) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7556954) q[0];
sx q[0];
rz(-1.9807666) q[0];
sx q[0];
rz(-2.9958821) q[0];
x q[1];
rz(2.5778779) q[2];
sx q[2];
rz(-1.5140805) q[2];
sx q[2];
rz(-2.4617755) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3562505) q[1];
sx q[1];
rz(-1.1467198) q[1];
sx q[1];
rz(2.8431358) q[1];
rz(-pi) q[2];
rz(0.032194897) q[3];
sx q[3];
rz(-2.1436084) q[3];
sx q[3];
rz(-2.9762852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.50600791) q[2];
sx q[2];
rz(-0.053746544) q[2];
sx q[2];
rz(-2.7347943) q[2];
rz(-0.16945101) q[3];
sx q[3];
rz(-0.5295161) q[3];
sx q[3];
rz(2.0690401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88103831) q[0];
sx q[0];
rz(-2.9091703) q[0];
sx q[0];
rz(3.1378003) q[0];
rz(0.077839851) q[1];
sx q[1];
rz(-0.65996116) q[1];
sx q[1];
rz(-0.30581623) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7821976) q[0];
sx q[0];
rz(-2.5919624) q[0];
sx q[0];
rz(-1.2122985) q[0];
x q[1];
rz(1.7051093) q[2];
sx q[2];
rz(-0.82342734) q[2];
sx q[2];
rz(0.01625492) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7485436) q[1];
sx q[1];
rz(-2.3094588) q[1];
sx q[1];
rz(-2.119675) q[1];
rz(3.0806957) q[3];
sx q[3];
rz(-1.506926) q[3];
sx q[3];
rz(2.8012432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1514312) q[2];
sx q[2];
rz(-1.7502681) q[2];
sx q[2];
rz(2.4988417) q[2];
rz(0.07240545) q[3];
sx q[3];
rz(-2.0824771) q[3];
sx q[3];
rz(-1.36093) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088242315) q[0];
sx q[0];
rz(-0.14277661) q[0];
sx q[0];
rz(2.9852168) q[0];
rz(0.0414255) q[1];
sx q[1];
rz(-0.62774575) q[1];
sx q[1];
rz(1.5511537) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0254733) q[0];
sx q[0];
rz(-1.0295233) q[0];
sx q[0];
rz(1.812029) q[0];
rz(-0.088038283) q[2];
sx q[2];
rz(-2.4836342) q[2];
sx q[2];
rz(-0.41706271) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2385905) q[1];
sx q[1];
rz(-2.6005473) q[1];
sx q[1];
rz(-0.21958406) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82156397) q[3];
sx q[3];
rz(-2.2872777) q[3];
sx q[3];
rz(2.2301073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40858832) q[2];
sx q[2];
rz(-0.53885794) q[2];
sx q[2];
rz(3.1206701) q[2];
rz(-0.18713348) q[3];
sx q[3];
rz(-2.9360866) q[3];
sx q[3];
rz(-0.11370295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1784096) q[0];
sx q[0];
rz(-0.33247501) q[0];
sx q[0];
rz(-2.7030429) q[0];
rz(-1.5248388) q[1];
sx q[1];
rz(-0.33477819) q[1];
sx q[1];
rz(-2.8964892) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49798508) q[0];
sx q[0];
rz(-1.6542572) q[0];
sx q[0];
rz(0.85026922) q[0];
rz(-2.7510277) q[2];
sx q[2];
rz(-2.7465944) q[2];
sx q[2];
rz(-0.67827889) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7286543) q[1];
sx q[1];
rz(-1.0333583) q[1];
sx q[1];
rz(-2.9793903) q[1];
rz(-2.80527) q[3];
sx q[3];
rz(-0.78236474) q[3];
sx q[3];
rz(2.4049644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54030067) q[2];
sx q[2];
rz(-2.7056594) q[2];
sx q[2];
rz(0.33176804) q[2];
rz(2.6541384) q[3];
sx q[3];
rz(-2.09477) q[3];
sx q[3];
rz(0.99307466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29193923) q[0];
sx q[0];
rz(-1.4537469) q[0];
sx q[0];
rz(2.3680903) q[0];
rz(1.1812814) q[1];
sx q[1];
rz(-0.14207323) q[1];
sx q[1];
rz(1.7519417) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8256102) q[0];
sx q[0];
rz(-1.9709236) q[0];
sx q[0];
rz(-1.9670301) q[0];
x q[1];
rz(-2.2393786) q[2];
sx q[2];
rz(-0.83219516) q[2];
sx q[2];
rz(2.8813643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4378499) q[1];
sx q[1];
rz(-1.4885167) q[1];
sx q[1];
rz(0.47474307) q[1];
x q[2];
rz(-1.3895274) q[3];
sx q[3];
rz(-1.3366404) q[3];
sx q[3];
rz(-0.93269809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9224077) q[2];
sx q[2];
rz(-1.2097404) q[2];
sx q[2];
rz(0.47214559) q[2];
rz(-1.8426497) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(-0.81056547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2209114) q[0];
sx q[0];
rz(-0.26316106) q[0];
sx q[0];
rz(0.26350185) q[0];
rz(-2.0384516) q[1];
sx q[1];
rz(-1.3145072) q[1];
sx q[1];
rz(-2.7679494) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.061122594) q[0];
sx q[0];
rz(-1.6310638) q[0];
sx q[0];
rz(-1.4979657) q[0];
rz(-0.84381004) q[2];
sx q[2];
rz(-1.4792974) q[2];
sx q[2];
rz(-1.3841656) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9571324) q[1];
sx q[1];
rz(-1.2533979) q[1];
sx q[1];
rz(0.98021345) q[1];
rz(-pi) q[2];
rz(0.75833894) q[3];
sx q[3];
rz(-2.3985574) q[3];
sx q[3];
rz(-0.61279994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0282447) q[2];
sx q[2];
rz(-2.9714606) q[2];
sx q[2];
rz(-2.587758) q[2];
rz(1.7438186) q[3];
sx q[3];
rz(-0.60269409) q[3];
sx q[3];
rz(0.3127313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58436191) q[0];
sx q[0];
rz(-1.0065684) q[0];
sx q[0];
rz(1.9118017) q[0];
rz(2.9025485) q[1];
sx q[1];
rz(-1.5115279) q[1];
sx q[1];
rz(2.8412433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8431138) q[0];
sx q[0];
rz(-1.4550147) q[0];
sx q[0];
rz(0.16364574) q[0];
rz(2.5744252) q[2];
sx q[2];
rz(-2.3544899) q[2];
sx q[2];
rz(1.2445104) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3546037) q[1];
sx q[1];
rz(-1.5711849) q[1];
sx q[1];
rz(1.555519) q[1];
rz(-3.1170397) q[3];
sx q[3];
rz(-0.85806393) q[3];
sx q[3];
rz(2.7803382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0099237) q[2];
sx q[2];
rz(-1.4574304) q[2];
sx q[2];
rz(2.8966676) q[2];
rz(2.6217672) q[3];
sx q[3];
rz(-2.2782785) q[3];
sx q[3];
rz(0.68827099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7898665) q[0];
sx q[0];
rz(-2.7930197) q[0];
sx q[0];
rz(-1.1630195) q[0];
rz(3.0746958) q[1];
sx q[1];
rz(-1.6480548) q[1];
sx q[1];
rz(1.012872) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5236008) q[0];
sx q[0];
rz(-2.3396684) q[0];
sx q[0];
rz(-2.2973934) q[0];
x q[1];
rz(-1.0032907) q[2];
sx q[2];
rz(-1.7219647) q[2];
sx q[2];
rz(-1.0790107) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7808395) q[1];
sx q[1];
rz(-1.2560802) q[1];
sx q[1];
rz(-0.34784045) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49533923) q[3];
sx q[3];
rz(-2.3916349) q[3];
sx q[3];
rz(1.0487674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11671994) q[2];
sx q[2];
rz(-2.1634384) q[2];
sx q[2];
rz(-0.32279521) q[2];
rz(-2.5358477) q[3];
sx q[3];
rz(-0.79234684) q[3];
sx q[3];
rz(-2.7927223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054319687) q[0];
sx q[0];
rz(-2.9948586) q[0];
sx q[0];
rz(-0.012454575) q[0];
rz(-0.74673486) q[1];
sx q[1];
rz(-2.2181999) q[1];
sx q[1];
rz(-2.8616203) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6387647) q[0];
sx q[0];
rz(-1.4584714) q[0];
sx q[0];
rz(3.0905484) q[0];
rz(-1.1367646) q[2];
sx q[2];
rz(-1.2723337) q[2];
sx q[2];
rz(3.0500183) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.758773) q[1];
sx q[1];
rz(-2.3291322) q[1];
sx q[1];
rz(0.42940213) q[1];
rz(1.9334698) q[3];
sx q[3];
rz(-1.8020523) q[3];
sx q[3];
rz(1.711292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.81194699) q[2];
sx q[2];
rz(-2.3880366) q[2];
sx q[2];
rz(-2.9296056) q[2];
rz(2.3181465) q[3];
sx q[3];
rz(-1.6971089) q[3];
sx q[3];
rz(-0.25920355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9987746) q[0];
sx q[0];
rz(-3.0840254) q[0];
sx q[0];
rz(2.4488191) q[0];
rz(2.5686) q[1];
sx q[1];
rz(-1.7968105) q[1];
sx q[1];
rz(2.7105892) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8739024) q[0];
sx q[0];
rz(-2.5866716) q[0];
sx q[0];
rz(3.0303427) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0153322) q[2];
sx q[2];
rz(-2.1241786) q[2];
sx q[2];
rz(-1.044342) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6083518) q[1];
sx q[1];
rz(-2.2679288) q[1];
sx q[1];
rz(1.9014386) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1206843) q[3];
sx q[3];
rz(-1.8326933) q[3];
sx q[3];
rz(0.49347116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4219605) q[2];
sx q[2];
rz(-2.8675291) q[2];
sx q[2];
rz(2.5813622) q[2];
rz(-0.46960056) q[3];
sx q[3];
rz(-2.738939) q[3];
sx q[3];
rz(-0.71389055) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2847168) q[0];
sx q[0];
rz(-1.7354043) q[0];
sx q[0];
rz(2.0866557) q[0];
rz(-2.3003385) q[1];
sx q[1];
rz(-1.1019191) q[1];
sx q[1];
rz(-0.046774653) q[1];
rz(-1.9404491) q[2];
sx q[2];
rz(-1.2751725) q[2];
sx q[2];
rz(-1.0775492) q[2];
rz(1.5359405) q[3];
sx q[3];
rz(-2.4823454) q[3];
sx q[3];
rz(-2.1128826) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
