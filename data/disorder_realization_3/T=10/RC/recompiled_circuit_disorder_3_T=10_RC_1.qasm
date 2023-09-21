OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(-0.48506919) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(5.8689868) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038923351) q[0];
sx q[0];
rz(-1.9063213) q[0];
sx q[0];
rz(-1.9468007) q[0];
rz(1.0010927) q[2];
sx q[2];
rz(-1.6506519) q[2];
sx q[2];
rz(2.9542838) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1297076) q[1];
sx q[1];
rz(-1.1457448) q[1];
sx q[1];
rz(3.0711864) q[1];
x q[2];
rz(0.36177735) q[3];
sx q[3];
rz(-2.8594115) q[3];
sx q[3];
rz(-2.2892945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2177314) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(3.1100173) q[2];
rz(1.2565553) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(-2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.3027705) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(2.9717428) q[0];
rz(0.70392144) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(-0.53952113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2182506) q[0];
sx q[0];
rz(-1.6174416) q[0];
sx q[0];
rz(1.1211066) q[0];
x q[1];
rz(1.9794481) q[2];
sx q[2];
rz(-2.2505629) q[2];
sx q[2];
rz(-2.0603927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62944618) q[1];
sx q[1];
rz(-1.0735895) q[1];
sx q[1];
rz(0.21558233) q[1];
x q[2];
rz(-2.2803454) q[3];
sx q[3];
rz(-1.1303139) q[3];
sx q[3];
rz(-1.4351821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4743621) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(-2.898522) q[2];
rz(0.66611755) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(1.8977785) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77984017) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(0.4483805) q[0];
rz(-1.386863) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(-0.2562491) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51107823) q[0];
sx q[0];
rz(-2.179638) q[0];
sx q[0];
rz(1.8947253) q[0];
x q[1];
rz(-0.82661144) q[2];
sx q[2];
rz(-2.5275143) q[2];
sx q[2];
rz(-0.6073063) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5970522) q[1];
sx q[1];
rz(-1.2433194) q[1];
sx q[1];
rz(-0.87045963) q[1];
rz(-0.17685299) q[3];
sx q[3];
rz(-1.6079418) q[3];
sx q[3];
rz(2.0166486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3391352) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(-2.5668872) q[2];
rz(-1.7859219) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(-1.2333966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.213585) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(-0.25948778) q[0];
rz(-1.9909987) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(2.4096699) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5314732) q[0];
sx q[0];
rz(-2.0391132) q[0];
sx q[0];
rz(-2.5584695) q[0];
x q[1];
rz(-2.4243083) q[2];
sx q[2];
rz(-1.6589763) q[2];
sx q[2];
rz(2.7862273) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.96958292) q[1];
sx q[1];
rz(-1.828555) q[1];
sx q[1];
rz(1.7715363) q[1];
rz(-pi) q[2];
rz(0.51283522) q[3];
sx q[3];
rz(-0.25179112) q[3];
sx q[3];
rz(0.73392111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(0.15110061) q[2];
rz(2.5949196) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(-0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916409) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(-1.8792101) q[0];
rz(-1.6732015) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(0.79777065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7724458) q[0];
sx q[0];
rz(-1.6909084) q[0];
sx q[0];
rz(-3.1390879) q[0];
x q[1];
rz(-1.7902137) q[2];
sx q[2];
rz(-1.8674388) q[2];
sx q[2];
rz(-2.9079633) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0610173) q[1];
sx q[1];
rz(-1.4772381) q[1];
sx q[1];
rz(0.94228014) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9005152) q[3];
sx q[3];
rz(-1.9922678) q[3];
sx q[3];
rz(-0.7082522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.795934) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(-0.30203715) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(-0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7323332) q[0];
sx q[0];
rz(-0.091826037) q[0];
sx q[0];
rz(-1.9858032) q[0];
rz(2.0571158) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(-0.070080431) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1228186) q[0];
sx q[0];
rz(-0.2304603) q[0];
sx q[0];
rz(-0.83968681) q[0];
rz(-2.5885133) q[2];
sx q[2];
rz(-2.2773909) q[2];
sx q[2];
rz(-1.4097708) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.97092123) q[1];
sx q[1];
rz(-2.1340003) q[1];
sx q[1];
rz(1.7621653) q[1];
x q[2];
rz(-3.1061884) q[3];
sx q[3];
rz(-1.4962713) q[3];
sx q[3];
rz(-2.7217334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8391116) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-2.5202259) q[2];
rz(1.7012043) q[3];
sx q[3];
rz(-0.50783235) q[3];
sx q[3];
rz(0.5733718) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0550585) q[0];
sx q[0];
rz(-1.7250412) q[0];
sx q[0];
rz(2.281718) q[0];
rz(-1.2043918) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(3.133657) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7746349) q[0];
sx q[0];
rz(-1.3152221) q[0];
sx q[0];
rz(-0.5704244) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2631049) q[2];
sx q[2];
rz(-1.9108859) q[2];
sx q[2];
rz(-2.0932587) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8994645) q[1];
sx q[1];
rz(-0.90485307) q[1];
sx q[1];
rz(-0.66023402) q[1];
rz(-pi) q[2];
rz(-2.2357335) q[3];
sx q[3];
rz(-2.6568036) q[3];
sx q[3];
rz(1.0672027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4456711) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(2.725214) q[2];
rz(-1.3683866) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798582) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(-2.7767048) q[0];
rz(-2.2015613) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(-1.6392802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68740326) q[0];
sx q[0];
rz(-2.8912376) q[0];
sx q[0];
rz(-0.040391163) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6697568) q[2];
sx q[2];
rz(-0.89315692) q[2];
sx q[2];
rz(3.133528) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2316206) q[1];
sx q[1];
rz(-0.93712229) q[1];
sx q[1];
rz(2.2909067) q[1];
rz(-1.7765462) q[3];
sx q[3];
rz(-1.4806517) q[3];
sx q[3];
rz(2.3004325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49729785) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(-1.0650744) q[2];
rz(-0.30125695) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(1.6206954) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97312462) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(-0.43564963) q[0];
rz(-1.7565953) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(-2.7246144) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57396736) q[0];
sx q[0];
rz(-1.6084533) q[0];
sx q[0];
rz(-2.2249939) q[0];
rz(2.142799) q[2];
sx q[2];
rz(-0.98888328) q[2];
sx q[2];
rz(1.7172161) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3560564) q[1];
sx q[1];
rz(-2.0691263) q[1];
sx q[1];
rz(-2.9836125) q[1];
rz(1.7008408) q[3];
sx q[3];
rz(-1.7107043) q[3];
sx q[3];
rz(3.1062982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0372662) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(-0.22582516) q[2];
rz(-0.2078235) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.726783) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(-1.5243994) q[0];
rz(-2.1879451) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(-1.3226002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029862558) q[0];
sx q[0];
rz(-0.39806453) q[0];
sx q[0];
rz(1.8882621) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8906891) q[2];
sx q[2];
rz(-2.0076027) q[2];
sx q[2];
rz(-1.3383588) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1982556) q[1];
sx q[1];
rz(-0.88611929) q[1];
sx q[1];
rz(-0.68590045) q[1];
x q[2];
rz(2.2262276) q[3];
sx q[3];
rz(-1.4487106) q[3];
sx q[3];
rz(-1.9742427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3502729) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(2.1255169) q[2];
rz(-1.919205) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(2.5861752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9983457) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(-1.3636419) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(-1.580711) q[2];
sx q[2];
rz(-2.0799939) q[2];
sx q[2];
rz(-1.6891198) q[2];
rz(0.090311269) q[3];
sx q[3];
rz(-1.3406546) q[3];
sx q[3];
rz(-1.8484074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];