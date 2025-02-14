OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.86201) q[0];
sx q[0];
rz(5.7481264) q[0];
sx q[0];
rz(9.6587107) q[0];
rz(0.84199953) q[1];
sx q[1];
rz(-1.7276126) q[1];
sx q[1];
rz(1.6230621) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6739963) q[0];
sx q[0];
rz(-2.5270871) q[0];
sx q[0];
rz(2.5835681) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9286381) q[2];
sx q[2];
rz(-1.1446272) q[2];
sx q[2];
rz(-1.0077241) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78511158) q[1];
sx q[1];
rz(-0.82725305) q[1];
sx q[1];
rz(-0.77036304) q[1];
x q[2];
rz(2.150334) q[3];
sx q[3];
rz(-2.391444) q[3];
sx q[3];
rz(2.0969491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6003549) q[2];
sx q[2];
rz(-1.7011832) q[2];
sx q[2];
rz(2.3360628) q[2];
rz(-0.39189288) q[3];
sx q[3];
rz(-1.7854179) q[3];
sx q[3];
rz(2.384757) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5517752) q[0];
sx q[0];
rz(-2.4601695) q[0];
sx q[0];
rz(-2.7031194) q[0];
rz(1.27502) q[1];
sx q[1];
rz(-1.4507111) q[1];
sx q[1];
rz(3.0335887) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94919449) q[0];
sx q[0];
rz(-0.11365232) q[0];
sx q[0];
rz(-1.3495007) q[0];
x q[1];
rz(-0.89620583) q[2];
sx q[2];
rz(-0.87088481) q[2];
sx q[2];
rz(-0.43255478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9671772) q[1];
sx q[1];
rz(-1.4870054) q[1];
sx q[1];
rz(-1.3166974) q[1];
x q[2];
rz(-2.9505902) q[3];
sx q[3];
rz(-2.0536978) q[3];
sx q[3];
rz(-2.1306899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7763623) q[2];
sx q[2];
rz(-0.47351101) q[2];
sx q[2];
rz(-0.11621172) q[2];
rz(-2.727437) q[3];
sx q[3];
rz(-2.0678346) q[3];
sx q[3];
rz(2.3381086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6666343) q[0];
sx q[0];
rz(-1.3131498) q[0];
sx q[0];
rz(0.83589244) q[0];
rz(-1.6235141) q[1];
sx q[1];
rz(-1.514879) q[1];
sx q[1];
rz(-1.2380884) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94592124) q[0];
sx q[0];
rz(-0.64650457) q[0];
sx q[0];
rz(2.2158428) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5290501) q[2];
sx q[2];
rz(-2.7296737) q[2];
sx q[2];
rz(-1.6430133) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.064621335) q[1];
sx q[1];
rz(-0.13084743) q[1];
sx q[1];
rz(-2.0276643) q[1];
rz(0.79583056) q[3];
sx q[3];
rz(-2.748877) q[3];
sx q[3];
rz(2.8098729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.28430024) q[2];
sx q[2];
rz(-2.9967873) q[2];
sx q[2];
rz(0.71883744) q[2];
rz(2.5607064) q[3];
sx q[3];
rz(-1.7234756) q[3];
sx q[3];
rz(0.7287997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.523664) q[0];
sx q[0];
rz(-1.5939465) q[0];
sx q[0];
rz(-1.1267598) q[0];
rz(-1.843938) q[1];
sx q[1];
rz(-1.490386) q[1];
sx q[1];
rz(1.0430956) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1820871) q[0];
sx q[0];
rz(-0.27217492) q[0];
sx q[0];
rz(-1.683781) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0342233) q[2];
sx q[2];
rz(-1.787428) q[2];
sx q[2];
rz(-0.41484268) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9981737) q[1];
sx q[1];
rz(-2.7980248) q[1];
sx q[1];
rz(-2.0433389) q[1];
x q[2];
rz(0.65535069) q[3];
sx q[3];
rz(-1.6387625) q[3];
sx q[3];
rz(0.36253906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.8635233) q[2];
sx q[2];
rz(-1.7570644) q[2];
sx q[2];
rz(-1.0370022) q[2];
rz(0.95748025) q[3];
sx q[3];
rz(-1.0629531) q[3];
sx q[3];
rz(2.3069042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1400414) q[0];
sx q[0];
rz(-0.20714864) q[0];
sx q[0];
rz(-3.0539883) q[0];
rz(0.43830782) q[1];
sx q[1];
rz(-2.0594845) q[1];
sx q[1];
rz(-1.8278488) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7676108) q[0];
sx q[0];
rz(-1.0896519) q[0];
sx q[0];
rz(0.26480459) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1042287) q[2];
sx q[2];
rz(-1.4722344) q[2];
sx q[2];
rz(-2.7600206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.61699762) q[1];
sx q[1];
rz(-1.4819078) q[1];
sx q[1];
rz(-1.7434381) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.143712) q[3];
sx q[3];
rz(-2.8393203) q[3];
sx q[3];
rz(-0.59461601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6614512) q[2];
sx q[2];
rz(-1.9037312) q[2];
sx q[2];
rz(-0.56337774) q[2];
rz(1.9208469) q[3];
sx q[3];
rz(-1.7505587) q[3];
sx q[3];
rz(-2.5105072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(3.0193943) q[0];
sx q[0];
rz(-0.65827426) q[0];
sx q[0];
rz(-0.011938183) q[0];
rz(2.7519233) q[1];
sx q[1];
rz(-0.7020815) q[1];
sx q[1];
rz(2.8964002) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6547873) q[0];
sx q[0];
rz(-1.1410895) q[0];
sx q[0];
rz(-0.91581099) q[0];
x q[1];
rz(1.7952314) q[2];
sx q[2];
rz(-1.8230453) q[2];
sx q[2];
rz(2.4750965) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4575544) q[1];
sx q[1];
rz(-2.0054617) q[1];
sx q[1];
rz(-2.5358389) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61556365) q[3];
sx q[3];
rz(-2.6608753) q[3];
sx q[3];
rz(1.2445039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.37990722) q[2];
sx q[2];
rz(-1.8319538) q[2];
sx q[2];
rz(-1.3690108) q[2];
rz(0.53064972) q[3];
sx q[3];
rz(-0.68415087) q[3];
sx q[3];
rz(2.0556889) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.948792) q[0];
sx q[0];
rz(-1.3230319) q[0];
sx q[0];
rz(-2.8444994) q[0];
rz(1.5104712) q[1];
sx q[1];
rz(-2.511697) q[1];
sx q[1];
rz(-2.7105601) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6490899) q[0];
sx q[0];
rz(-1.8592872) q[0];
sx q[0];
rz(-1.7899075) q[0];
x q[1];
rz(-2.9210429) q[2];
sx q[2];
rz(-0.23633453) q[2];
sx q[2];
rz(0.55921171) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7546852) q[1];
sx q[1];
rz(-1.0656989) q[1];
sx q[1];
rz(2.7688857) q[1];
x q[2];
rz(-2.0991574) q[3];
sx q[3];
rz(-1.8260806) q[3];
sx q[3];
rz(1.135889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.83583528) q[2];
sx q[2];
rz(-2.367986) q[2];
sx q[2];
rz(-2.8744899) q[2];
rz(3.109572) q[3];
sx q[3];
rz(-1.9710385) q[3];
sx q[3];
rz(-0.10572461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778462) q[0];
sx q[0];
rz(-1.2910605) q[0];
sx q[0];
rz(0.63999501) q[0];
rz(-0.85482875) q[1];
sx q[1];
rz(-1.4850441) q[1];
sx q[1];
rz(-1.7291501) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4025637) q[0];
sx q[0];
rz(-1.1629346) q[0];
sx q[0];
rz(-0.5824851) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.090254) q[2];
sx q[2];
rz(-2.0624537) q[2];
sx q[2];
rz(-1.7887539) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1302274) q[1];
sx q[1];
rz(-1.9107358) q[1];
sx q[1];
rz(1.4814427) q[1];
rz(2.2858587) q[3];
sx q[3];
rz(-0.54045709) q[3];
sx q[3];
rz(0.021322535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60014805) q[2];
sx q[2];
rz(-1.4358127) q[2];
sx q[2];
rz(-2.7195462) q[2];
rz(-0.47354928) q[3];
sx q[3];
rz(-2.519042) q[3];
sx q[3];
rz(0.1951018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37050978) q[0];
sx q[0];
rz(-2.2042553) q[0];
sx q[0];
rz(-1.7899845) q[0];
rz(-0.31556684) q[1];
sx q[1];
rz(-1.4548929) q[1];
sx q[1];
rz(1.1531166) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4723268) q[0];
sx q[0];
rz(-1.6231704) q[0];
sx q[0];
rz(-3.0906424) q[0];
x q[1];
rz(3.0598203) q[2];
sx q[2];
rz(-2.6047629) q[2];
sx q[2];
rz(2.5935136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2994581) q[1];
sx q[1];
rz(-0.99255764) q[1];
sx q[1];
rz(-1.7229573) q[1];
x q[2];
rz(1.159129) q[3];
sx q[3];
rz(-1.921706) q[3];
sx q[3];
rz(-2.0716435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6844668) q[2];
sx q[2];
rz(-1.212333) q[2];
sx q[2];
rz(-2.7678164) q[2];
rz(-1.2260381) q[3];
sx q[3];
rz(-2.2235179) q[3];
sx q[3];
rz(1.8019684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95947295) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(1.8208338) q[0];
rz(0.93718115) q[1];
sx q[1];
rz(-1.8056185) q[1];
sx q[1];
rz(-2.6722867) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.824911) q[0];
sx q[0];
rz(-1.4535722) q[0];
sx q[0];
rz(-3.0646695) q[0];
rz(1.2861757) q[2];
sx q[2];
rz(-0.60065833) q[2];
sx q[2];
rz(-0.17135581) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0476956) q[1];
sx q[1];
rz(-0.17054955) q[1];
sx q[1];
rz(-2.6490414) q[1];
rz(-0.68113459) q[3];
sx q[3];
rz(-0.38060846) q[3];
sx q[3];
rz(-1.3342302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4361973) q[2];
sx q[2];
rz(-0.86279482) q[2];
sx q[2];
rz(-1.2104642) q[2];
rz(-2.5734899) q[3];
sx q[3];
rz(-1.3848687) q[3];
sx q[3];
rz(-0.97584045) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2198915) q[0];
sx q[0];
rz(-0.85048631) q[0];
sx q[0];
rz(1.462107) q[0];
rz(-0.89314356) q[1];
sx q[1];
rz(-1.3124663) q[1];
sx q[1];
rz(1.5193473) q[1];
rz(2.2371815) q[2];
sx q[2];
rz(-0.36199649) q[2];
sx q[2];
rz(0.71142759) q[2];
rz(-2.1410971) q[3];
sx q[3];
rz(-1.9420997) q[3];
sx q[3];
rz(1.6122769) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
