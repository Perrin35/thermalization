OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47867632) q[0];
sx q[0];
rz(-1.588568) q[0];
sx q[0];
rz(-1.2583865) q[0];
rz(-1.8796743) q[1];
sx q[1];
rz(-0.51845297) q[1];
sx q[1];
rz(0.060103091) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10190554) q[0];
sx q[0];
rz(-1.5400258) q[0];
sx q[0];
rz(0.090307856) q[0];
rz(-pi) q[1];
rz(0.2807887) q[2];
sx q[2];
rz(-2.6323689) q[2];
sx q[2];
rz(1.6564684) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2138252) q[1];
sx q[1];
rz(-2.7672184) q[1];
sx q[1];
rz(-1.8273749) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46334907) q[3];
sx q[3];
rz(-1.0245205) q[3];
sx q[3];
rz(2.0093371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.95187676) q[2];
sx q[2];
rz(-2.3592301) q[2];
sx q[2];
rz(-2.961109) q[2];
rz(-0.30432025) q[3];
sx q[3];
rz(-0.92420998) q[3];
sx q[3];
rz(2.9144104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75905281) q[0];
sx q[0];
rz(-2.5681684) q[0];
sx q[0];
rz(3.0376814) q[0];
rz(0.78481627) q[1];
sx q[1];
rz(-1.3094614) q[1];
sx q[1];
rz(0.58194247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4020445) q[0];
sx q[0];
rz(-0.53240896) q[0];
sx q[0];
rz(-3.0976712) q[0];
x q[1];
rz(1.9230827) q[2];
sx q[2];
rz(-1.0880044) q[2];
sx q[2];
rz(2.5329075) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84060639) q[1];
sx q[1];
rz(-2.0475629) q[1];
sx q[1];
rz(0.029164) q[1];
x q[2];
rz(-1.1725015) q[3];
sx q[3];
rz(-1.1886667) q[3];
sx q[3];
rz(2.7492461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4265784) q[2];
sx q[2];
rz(-0.38807401) q[2];
sx q[2];
rz(-2.1132054) q[2];
rz(-2.5905124) q[3];
sx q[3];
rz(-1.9469399) q[3];
sx q[3];
rz(2.6330131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55117637) q[0];
sx q[0];
rz(-1.6352147) q[0];
sx q[0];
rz(1.2580385) q[0];
rz(2.538077) q[1];
sx q[1];
rz(-1.3110833) q[1];
sx q[1];
rz(-0.18361941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1574137) q[0];
sx q[0];
rz(-1.8350198) q[0];
sx q[0];
rz(1.9963032) q[0];
x q[1];
rz(0.48086353) q[2];
sx q[2];
rz(-1.6893975) q[2];
sx q[2];
rz(3.1080217) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2752876) q[1];
sx q[1];
rz(-0.86967378) q[1];
sx q[1];
rz(-0.57544586) q[1];
rz(1.5294936) q[3];
sx q[3];
rz(-1.777711) q[3];
sx q[3];
rz(-2.4504091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52643481) q[2];
sx q[2];
rz(-1.7650471) q[2];
sx q[2];
rz(-0.60205013) q[2];
rz(2.708882) q[3];
sx q[3];
rz(-1.5940462) q[3];
sx q[3];
rz(1.6656779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.3848307) q[0];
sx q[0];
rz(-2.6669406) q[0];
sx q[0];
rz(-0.62622825) q[0];
rz(-1.8734044) q[1];
sx q[1];
rz(-1.405193) q[1];
sx q[1];
rz(-2.7558806) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69067467) q[0];
sx q[0];
rz(-1.1407778) q[0];
sx q[0];
rz(-2.3818092) q[0];
x q[1];
rz(1.1168295) q[2];
sx q[2];
rz(-1.0359284) q[2];
sx q[2];
rz(-1.7810389) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.205769) q[1];
sx q[1];
rz(-0.21142928) q[1];
sx q[1];
rz(0.32306674) q[1];
rz(-1.3703501) q[3];
sx q[3];
rz(-0.55209898) q[3];
sx q[3];
rz(-0.31262661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4422153) q[2];
sx q[2];
rz(-1.135004) q[2];
sx q[2];
rz(2.8507612) q[2];
rz(1.95131) q[3];
sx q[3];
rz(-0.61605993) q[3];
sx q[3];
rz(-1.0440089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.46665835) q[0];
sx q[0];
rz(-0.062677296) q[0];
sx q[0];
rz(-1.4516996) q[0];
rz(0.85707227) q[1];
sx q[1];
rz(-1.8430201) q[1];
sx q[1];
rz(-2.7136386) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4820517) q[0];
sx q[0];
rz(-0.45705802) q[0];
sx q[0];
rz(-1.8607451) q[0];
rz(-2.0796135) q[2];
sx q[2];
rz(-2.2428838) q[2];
sx q[2];
rz(0.88269688) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9792069) q[1];
sx q[1];
rz(-2.2507994) q[1];
sx q[1];
rz(1.1080145) q[1];
rz(-pi) q[2];
rz(-2.4492743) q[3];
sx q[3];
rz(-2.025617) q[3];
sx q[3];
rz(-2.2934283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9383303) q[2];
sx q[2];
rz(-2.0945956) q[2];
sx q[2];
rz(2.9948575) q[2];
rz(-2.9444368) q[3];
sx q[3];
rz(-1.4898172) q[3];
sx q[3];
rz(-2.3728235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7305304) q[0];
sx q[0];
rz(-2.125232) q[0];
sx q[0];
rz(-2.9732669) q[0];
rz(-2.7492145) q[1];
sx q[1];
rz(-0.43586755) q[1];
sx q[1];
rz(-1.4515152) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1010676) q[0];
sx q[0];
rz(-2.7042964) q[0];
sx q[0];
rz(1.5202724) q[0];
x q[1];
rz(-2.8807202) q[2];
sx q[2];
rz(-2.2745273) q[2];
sx q[2];
rz(2.1056625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5785006) q[1];
sx q[1];
rz(-2.3804286) q[1];
sx q[1];
rz(-2.9369686) q[1];
x q[2];
rz(-2.3528162) q[3];
sx q[3];
rz(-1.8982072) q[3];
sx q[3];
rz(-2.56853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.77439857) q[2];
sx q[2];
rz(-2.6937679) q[2];
sx q[2];
rz(1.1773342) q[2];
rz(-0.92159739) q[3];
sx q[3];
rz(-2.2955743) q[3];
sx q[3];
rz(0.53068501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0717936) q[0];
sx q[0];
rz(-1.8853747) q[0];
sx q[0];
rz(-1.0569093) q[0];
rz(-2.9131556) q[1];
sx q[1];
rz(-0.7799131) q[1];
sx q[1];
rz(-0.25100073) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0989444) q[0];
sx q[0];
rz(-1.8269208) q[0];
sx q[0];
rz(0.29447181) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86871882) q[2];
sx q[2];
rz(-2.0755322) q[2];
sx q[2];
rz(-2.2622893) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7556762) q[1];
sx q[1];
rz(-1.7337013) q[1];
sx q[1];
rz(-1.7079123) q[1];
rz(-2.6282309) q[3];
sx q[3];
rz(-1.2436449) q[3];
sx q[3];
rz(2.1029496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.935427) q[2];
sx q[2];
rz(-2.3212104) q[2];
sx q[2];
rz(2.7247735) q[2];
rz(-2.5721512) q[3];
sx q[3];
rz(-1.1007525) q[3];
sx q[3];
rz(-2.4429564) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57132974) q[0];
sx q[0];
rz(-1.2705734) q[0];
sx q[0];
rz(2.6369693) q[0];
rz(-2.8847671) q[1];
sx q[1];
rz(-0.80373126) q[1];
sx q[1];
rz(1.1508734) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52792943) q[0];
sx q[0];
rz(-1.3656034) q[0];
sx q[0];
rz(2.0733207) q[0];
rz(-pi) q[1];
rz(2.0382284) q[2];
sx q[2];
rz(-0.69004493) q[2];
sx q[2];
rz(-2.0375117) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.098365) q[1];
sx q[1];
rz(-2.4154764) q[1];
sx q[1];
rz(-2.4642706) q[1];
x q[2];
rz(-0.34567771) q[3];
sx q[3];
rz(-2.7015903) q[3];
sx q[3];
rz(2.1328762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5095832) q[2];
sx q[2];
rz(-2.7183967) q[2];
sx q[2];
rz(2.7022341) q[2];
rz(-1.1257233) q[3];
sx q[3];
rz(-1.537354) q[3];
sx q[3];
rz(-1.2996947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87089649) q[0];
sx q[0];
rz(-0.82042158) q[0];
sx q[0];
rz(0.46723715) q[0];
rz(2.1029419) q[1];
sx q[1];
rz(-2.6324582) q[1];
sx q[1];
rz(1.4899303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91083807) q[0];
sx q[0];
rz(-1.0808792) q[0];
sx q[0];
rz(1.7696437) q[0];
rz(-1.4760255) q[2];
sx q[2];
rz(-1.4247923) q[2];
sx q[2];
rz(2.3340747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0584152) q[1];
sx q[1];
rz(-1.5442368) q[1];
sx q[1];
rz(0.81574622) q[1];
rz(2.5760004) q[3];
sx q[3];
rz(-1.9533691) q[3];
sx q[3];
rz(2.3747834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44522875) q[2];
sx q[2];
rz(-0.35564056) q[2];
sx q[2];
rz(1.5544372) q[2];
rz(-2.1571531) q[3];
sx q[3];
rz(-2.2319904) q[3];
sx q[3];
rz(1.4136774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5928891) q[0];
sx q[0];
rz(-2.2991572) q[0];
sx q[0];
rz(2.5133251) q[0];
rz(2.9382622) q[1];
sx q[1];
rz(-1.1140099) q[1];
sx q[1];
rz(-2.5471953) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5448611) q[0];
sx q[0];
rz(-0.98698915) q[0];
sx q[0];
rz(-0.24703276) q[0];
x q[1];
rz(-2.6028676) q[2];
sx q[2];
rz(-1.7598514) q[2];
sx q[2];
rz(-0.75720398) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0600784) q[1];
sx q[1];
rz(-1.8412672) q[1];
sx q[1];
rz(-2.8287776) q[1];
rz(-pi) q[2];
rz(0.61305586) q[3];
sx q[3];
rz(-1.5125015) q[3];
sx q[3];
rz(-1.6976274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2079042) q[2];
sx q[2];
rz(-2.4173357) q[2];
sx q[2];
rz(0.17624632) q[2];
rz(-1.3147973) q[3];
sx q[3];
rz(-0.81614554) q[3];
sx q[3];
rz(-0.76752457) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15277302) q[0];
sx q[0];
rz(-1.7012699) q[0];
sx q[0];
rz(-1.758601) q[0];
rz(3.026961) q[1];
sx q[1];
rz(-0.75420598) q[1];
sx q[1];
rz(2.4892714) q[1];
rz(0.23723142) q[2];
sx q[2];
rz(-1.5330412) q[2];
sx q[2];
rz(1.588149) q[2];
rz(0.38269855) q[3];
sx q[3];
rz(-1.1094339) q[3];
sx q[3];
rz(-1.4239428) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
