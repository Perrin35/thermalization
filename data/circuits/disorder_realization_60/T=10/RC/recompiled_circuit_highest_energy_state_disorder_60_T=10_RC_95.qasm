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
rz(0.76437104) q[0];
sx q[0];
rz(1.8005014) q[0];
sx q[0];
rz(10.330248) q[0];
rz(1.5068997) q[1];
sx q[1];
rz(-2.1807179) q[1];
sx q[1];
rz(-1.6406055) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1849821) q[0];
sx q[0];
rz(-2.2652103) q[0];
sx q[0];
rz(-2.882769) q[0];
rz(1.2086966) q[2];
sx q[2];
rz(-2.675867) q[2];
sx q[2];
rz(0.036368792) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6480744) q[1];
sx q[1];
rz(-1.5474802) q[1];
sx q[1];
rz(3.0991395) q[1];
rz(-2.7952173) q[3];
sx q[3];
rz(-1.0358255) q[3];
sx q[3];
rz(2.3807824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2291439) q[2];
sx q[2];
rz(-1.7488166) q[2];
sx q[2];
rz(2.8196715) q[2];
rz(0.95494444) q[3];
sx q[3];
rz(-2.3117282) q[3];
sx q[3];
rz(-1.9979075) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8973273) q[0];
sx q[0];
rz(-2.0877512) q[0];
sx q[0];
rz(0.51189297) q[0];
rz(1.5533252) q[1];
sx q[1];
rz(-0.48738185) q[1];
sx q[1];
rz(2.0194676) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13817912) q[0];
sx q[0];
rz(-1.2597359) q[0];
sx q[0];
rz(-1.6775234) q[0];
x q[1];
rz(-2.000359) q[2];
sx q[2];
rz(-1.6074174) q[2];
sx q[2];
rz(-3.1319194) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0341238) q[1];
sx q[1];
rz(-1.9975047) q[1];
sx q[1];
rz(-0.36860768) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53776017) q[3];
sx q[3];
rz(-0.7825635) q[3];
sx q[3];
rz(-1.3937221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.7684218) q[2];
sx q[2];
rz(-1.3139407) q[2];
sx q[2];
rz(-2.6785417) q[2];
rz(-0.57524663) q[3];
sx q[3];
rz(-1.3653711) q[3];
sx q[3];
rz(0.099253207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7333882) q[0];
sx q[0];
rz(-0.84771228) q[0];
sx q[0];
rz(0.43011618) q[0];
rz(-2.6759713) q[1];
sx q[1];
rz(-0.72048134) q[1];
sx q[1];
rz(0.63708416) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4642536) q[0];
sx q[0];
rz(-1.5298944) q[0];
sx q[0];
rz(3.1275355) q[0];
x q[1];
rz(-1.5445721) q[2];
sx q[2];
rz(-0.78557149) q[2];
sx q[2];
rz(2.3674813) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8608537) q[1];
sx q[1];
rz(-2.2855084) q[1];
sx q[1];
rz(1.2398941) q[1];
rz(-pi) q[2];
rz(-2.5175321) q[3];
sx q[3];
rz(-2.5284323) q[3];
sx q[3];
rz(-1.1414736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35885262) q[2];
sx q[2];
rz(-2.095486) q[2];
sx q[2];
rz(2.4165238) q[2];
rz(1.8105761) q[3];
sx q[3];
rz(-2.1264117) q[3];
sx q[3];
rz(2.5031808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8722039) q[0];
sx q[0];
rz(-1.6870455) q[0];
sx q[0];
rz(1.697502) q[0];
rz(-3.0929502) q[1];
sx q[1];
rz(-1.3261869) q[1];
sx q[1];
rz(-2.8526502) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3061482) q[0];
sx q[0];
rz(-1.939538) q[0];
sx q[0];
rz(-2.8822417) q[0];
rz(-1.4202576) q[2];
sx q[2];
rz(-1.5988013) q[2];
sx q[2];
rz(1.3893407) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1370586) q[1];
sx q[1];
rz(-1.6797795) q[1];
sx q[1];
rz(-1.7458899) q[1];
rz(-pi) q[2];
rz(3.0722202) q[3];
sx q[3];
rz(-2.4751304) q[3];
sx q[3];
rz(1.3795167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.284953) q[2];
sx q[2];
rz(-2.608947) q[2];
sx q[2];
rz(2.7613769) q[2];
rz(2.5034261) q[3];
sx q[3];
rz(-1.2297945) q[3];
sx q[3];
rz(-2.0130472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83288348) q[0];
sx q[0];
rz(-2.5888011) q[0];
sx q[0];
rz(2.047245) q[0];
rz(-2.265918) q[1];
sx q[1];
rz(-0.62962571) q[1];
sx q[1];
rz(-1.099115) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5364781) q[0];
sx q[0];
rz(-2.7301359) q[0];
sx q[0];
rz(2.7921929) q[0];
rz(-pi) q[1];
rz(1.9149295) q[2];
sx q[2];
rz(-1.2467071) q[2];
sx q[2];
rz(-1.7050755) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.588394) q[1];
sx q[1];
rz(-2.1171682) q[1];
sx q[1];
rz(-2.9562034) q[1];
x q[2];
rz(-1.7775675) q[3];
sx q[3];
rz(-2.4768157) q[3];
sx q[3];
rz(2.1867276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2431474) q[2];
sx q[2];
rz(-1.7869608) q[2];
sx q[2];
rz(0.97664991) q[2];
rz(0.53168932) q[3];
sx q[3];
rz(-0.44638005) q[3];
sx q[3];
rz(-0.71438742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0584745) q[0];
sx q[0];
rz(-1.1019022) q[0];
sx q[0];
rz(1.2716768) q[0];
rz(-0.42287982) q[1];
sx q[1];
rz(-1.5300749) q[1];
sx q[1];
rz(-1.5333102) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6278224) q[0];
sx q[0];
rz(-1.5830399) q[0];
sx q[0];
rz(-1.5041385) q[0];
rz(-1.4752059) q[2];
sx q[2];
rz(-3.0431192) q[2];
sx q[2];
rz(-2.5661039) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9257838) q[1];
sx q[1];
rz(-1.5152351) q[1];
sx q[1];
rz(0.30345281) q[1];
rz(-pi) q[2];
rz(1.815237) q[3];
sx q[3];
rz(-2.3688753) q[3];
sx q[3];
rz(-2.6561007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5292458) q[2];
sx q[2];
rz(-1.8491448) q[2];
sx q[2];
rz(0.38522729) q[2];
rz(-3.0569844) q[3];
sx q[3];
rz(-0.43360964) q[3];
sx q[3];
rz(0.87578526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92152921) q[0];
sx q[0];
rz(-1.055287) q[0];
sx q[0];
rz(0.39749843) q[0];
rz(2.6037604) q[1];
sx q[1];
rz(-2.718524) q[1];
sx q[1];
rz(0.0040815512) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794681) q[0];
sx q[0];
rz(-1.4433658) q[0];
sx q[0];
rz(0.80029933) q[0];
rz(0.40553826) q[2];
sx q[2];
rz(-2.2161762) q[2];
sx q[2];
rz(-2.1359512) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.34985182) q[1];
sx q[1];
rz(-2.428973) q[1];
sx q[1];
rz(-2.5373759) q[1];
rz(-pi) q[2];
rz(-2.2479731) q[3];
sx q[3];
rz(-2.1359518) q[3];
sx q[3];
rz(2.2806185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5092545) q[2];
sx q[2];
rz(-1.1804487) q[2];
sx q[2];
rz(-2.9285367) q[2];
rz(-0.80900711) q[3];
sx q[3];
rz(-3.1000948) q[3];
sx q[3];
rz(1.2808778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1236561) q[0];
sx q[0];
rz(-1.9391215) q[0];
sx q[0];
rz(-1.1630455) q[0];
rz(0.2991547) q[1];
sx q[1];
rz(-1.9367633) q[1];
sx q[1];
rz(0.97602731) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3523063) q[0];
sx q[0];
rz(-1.8160161) q[0];
sx q[0];
rz(-1.4568491) q[0];
rz(0.15248044) q[2];
sx q[2];
rz(-1.7465542) q[2];
sx q[2];
rz(0.7712785) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27777003) q[1];
sx q[1];
rz(-0.43102095) q[1];
sx q[1];
rz(-0.96993877) q[1];
rz(2.7175432) q[3];
sx q[3];
rz(-1.3433787) q[3];
sx q[3];
rz(0.018674803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8396847) q[2];
sx q[2];
rz(-2.8161616) q[2];
sx q[2];
rz(-0.1304661) q[2];
rz(1.8179551) q[3];
sx q[3];
rz(-1.2490844) q[3];
sx q[3];
rz(-0.29449335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.2015304) q[0];
sx q[0];
rz(-0.37684965) q[0];
sx q[0];
rz(-1.6424204) q[0];
rz(-2.6014853) q[1];
sx q[1];
rz(-2.3441548) q[1];
sx q[1];
rz(1.3444208) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8136052) q[0];
sx q[0];
rz(-1.9623956) q[0];
sx q[0];
rz(2.8285976) q[0];
rz(-0.62144582) q[2];
sx q[2];
rz(-2.0315731) q[2];
sx q[2];
rz(-0.64475343) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0186179) q[1];
sx q[1];
rz(-0.6280762) q[1];
sx q[1];
rz(-2.0043892) q[1];
rz(2.2910883) q[3];
sx q[3];
rz(-2.0652186) q[3];
sx q[3];
rz(1.1098343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53238791) q[2];
sx q[2];
rz(-2.8049073) q[2];
sx q[2];
rz(1.6877635) q[2];
rz(0.42803556) q[3];
sx q[3];
rz(-1.5902218) q[3];
sx q[3];
rz(-1.9868896) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5880599) q[0];
sx q[0];
rz(-2.689671) q[0];
sx q[0];
rz(-1.6499299) q[0];
rz(2.5904169) q[1];
sx q[1];
rz(-1.0124413) q[1];
sx q[1];
rz(0.59250441) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4649749) q[0];
sx q[0];
rz(-1.3335557) q[0];
sx q[0];
rz(2.854191) q[0];
rz(-pi) q[1];
rz(0.5430605) q[2];
sx q[2];
rz(-2.0351699) q[2];
sx q[2];
rz(-2.4891702) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7043276) q[1];
sx q[1];
rz(-2.6710837) q[1];
sx q[1];
rz(-0.60948845) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6423934) q[3];
sx q[3];
rz(-1.6147524) q[3];
sx q[3];
rz(1.6658962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7117915) q[2];
sx q[2];
rz(-2.1000803) q[2];
sx q[2];
rz(-0.05833021) q[2];
rz(-2.77099) q[3];
sx q[3];
rz(-0.27675089) q[3];
sx q[3];
rz(-3.1378194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8534828) q[0];
sx q[0];
rz(-1.3567038) q[0];
sx q[0];
rz(0.10467341) q[0];
rz(1.7755605) q[1];
sx q[1];
rz(-0.80148253) q[1];
sx q[1];
rz(-2.186224) q[1];
rz(-2.7544207) q[2];
sx q[2];
rz(-2.5324814) q[2];
sx q[2];
rz(-1.0003288) q[2];
rz(2.3351135) q[3];
sx q[3];
rz(-2.4462593) q[3];
sx q[3];
rz(-1.367955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
