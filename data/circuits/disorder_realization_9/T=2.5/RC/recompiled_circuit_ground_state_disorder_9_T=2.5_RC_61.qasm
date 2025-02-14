OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0661434) q[0];
sx q[0];
rz(-2.0976522) q[0];
sx q[0];
rz(-0.010308417) q[0];
rz(-0.97996867) q[1];
sx q[1];
rz(-1.6966532) q[1];
sx q[1];
rz(0.57656062) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.071237) q[0];
sx q[0];
rz(-2.8094387) q[0];
sx q[0];
rz(-2.7624056) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61277436) q[2];
sx q[2];
rz(-0.98721993) q[2];
sx q[2];
rz(0.60223168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3763872) q[1];
sx q[1];
rz(-2.58316) q[1];
sx q[1];
rz(0.33513481) q[1];
x q[2];
rz(1.3862085) q[3];
sx q[3];
rz(-2.3999676) q[3];
sx q[3];
rz(-3.1390934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6171241) q[2];
sx q[2];
rz(-1.4602129) q[2];
sx q[2];
rz(-1.8560393) q[2];
rz(-1.589795) q[3];
sx q[3];
rz(-2.0701305) q[3];
sx q[3];
rz(-1.0514528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0153506) q[0];
sx q[0];
rz(-1.9168357) q[0];
sx q[0];
rz(2.8958564) q[0];
rz(-2.083678) q[1];
sx q[1];
rz(-1.3163047) q[1];
sx q[1];
rz(-0.33448514) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6545821) q[0];
sx q[0];
rz(-1.0791057) q[0];
sx q[0];
rz(-2.0915178) q[0];
rz(-pi) q[1];
rz(2.7203015) q[2];
sx q[2];
rz(-2.3347046) q[2];
sx q[2];
rz(0.56299201) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3547438) q[1];
sx q[1];
rz(-1.9450099) q[1];
sx q[1];
rz(2.0207062) q[1];
rz(0.40615079) q[3];
sx q[3];
rz(-2.5037919) q[3];
sx q[3];
rz(2.8075308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1841396) q[2];
sx q[2];
rz(-0.89749557) q[2];
sx q[2];
rz(1.3857566) q[2];
rz(2.1615248) q[3];
sx q[3];
rz(-0.46531427) q[3];
sx q[3];
rz(3.0906299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053213) q[0];
sx q[0];
rz(-0.956981) q[0];
sx q[0];
rz(-1.7514239) q[0];
rz(-2.7888489) q[1];
sx q[1];
rz(-0.91068641) q[1];
sx q[1];
rz(2.5211451) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.016997) q[0];
sx q[0];
rz(-3.0244138) q[0];
sx q[0];
rz(-2.4524053) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0376126) q[2];
sx q[2];
rz(-0.98538387) q[2];
sx q[2];
rz(2.3217161) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31252334) q[1];
sx q[1];
rz(-1.2770953) q[1];
sx q[1];
rz(1.9109265) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51692112) q[3];
sx q[3];
rz(-2.1955588) q[3];
sx q[3];
rz(-2.3464835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1227526) q[2];
sx q[2];
rz(-0.41308013) q[2];
sx q[2];
rz(1.2472461) q[2];
rz(-0.16417575) q[3];
sx q[3];
rz(-1.7069867) q[3];
sx q[3];
rz(-2.5119787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4267047) q[0];
sx q[0];
rz(-2.1737104) q[0];
sx q[0];
rz(-1.8545275) q[0];
rz(0.60028589) q[1];
sx q[1];
rz(-1.7900107) q[1];
sx q[1];
rz(-2.2379025) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716887) q[0];
sx q[0];
rz(-1.6780417) q[0];
sx q[0];
rz(1.6374169) q[0];
x q[1];
rz(1.5026592) q[2];
sx q[2];
rz(-1.9499511) q[2];
sx q[2];
rz(0.013163002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0541573) q[1];
sx q[1];
rz(-1.8241) q[1];
sx q[1];
rz(-1.3852055) q[1];
rz(-pi) q[2];
rz(2.4900808) q[3];
sx q[3];
rz(-0.50262302) q[3];
sx q[3];
rz(0.16829695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51957447) q[2];
sx q[2];
rz(-2.3081686) q[2];
sx q[2];
rz(0.77345094) q[2];
rz(-1.2889688) q[3];
sx q[3];
rz(-0.52923146) q[3];
sx q[3];
rz(2.4066063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2499823) q[0];
sx q[0];
rz(-2.1365428) q[0];
sx q[0];
rz(-2.6494359) q[0];
rz(-0.53897578) q[1];
sx q[1];
rz(-1.69918) q[1];
sx q[1];
rz(-1.0999365) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66456074) q[0];
sx q[0];
rz(-2.7072621) q[0];
sx q[0];
rz(0.6566027) q[0];
rz(-3.1308168) q[2];
sx q[2];
rz(-0.95913619) q[2];
sx q[2];
rz(-1.7623368) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.392994) q[1];
sx q[1];
rz(-1.6637319) q[1];
sx q[1];
rz(-2.3307021) q[1];
x q[2];
rz(2.5650129) q[3];
sx q[3];
rz(-0.55906536) q[3];
sx q[3];
rz(-0.41043974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34053549) q[2];
sx q[2];
rz(-0.33581442) q[2];
sx q[2];
rz(2.578793) q[2];
rz(-2.4417012) q[3];
sx q[3];
rz(-1.8507345) q[3];
sx q[3];
rz(2.8904397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3444779) q[0];
sx q[0];
rz(-2.4176702) q[0];
sx q[0];
rz(2.7864454) q[0];
rz(1.2190602) q[1];
sx q[1];
rz(-1.9025758) q[1];
sx q[1];
rz(2.6893137) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1817976) q[0];
sx q[0];
rz(-0.99645146) q[0];
sx q[0];
rz(2.6622165) q[0];
rz(-pi) q[1];
rz(-2.7554871) q[2];
sx q[2];
rz(-1.9506644) q[2];
sx q[2];
rz(0.42797336) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.86231316) q[1];
sx q[1];
rz(-1.24839) q[1];
sx q[1];
rz(-3.1161948) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5127453) q[3];
sx q[3];
rz(-1.4553242) q[3];
sx q[3];
rz(-2.2001298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63048116) q[2];
sx q[2];
rz(-0.52383542) q[2];
sx q[2];
rz(-0.13709489) q[2];
rz(-1.5824205) q[3];
sx q[3];
rz(-1.987792) q[3];
sx q[3];
rz(2.3649575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1906076) q[0];
sx q[0];
rz(-0.76222104) q[0];
sx q[0];
rz(1.5964339) q[0];
rz(0.51721382) q[1];
sx q[1];
rz(-1.1511753) q[1];
sx q[1];
rz(-2.1545765) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68970976) q[0];
sx q[0];
rz(-1.5919884) q[0];
sx q[0];
rz(3.0107493) q[0];
x q[1];
rz(-0.9632684) q[2];
sx q[2];
rz(-1.844256) q[2];
sx q[2];
rz(2.4475522) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4728022) q[1];
sx q[1];
rz(-1.3061532) q[1];
sx q[1];
rz(-1.3331148) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8214885) q[3];
sx q[3];
rz(-2.8134228) q[3];
sx q[3];
rz(1.8712957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78642693) q[2];
sx q[2];
rz(-1.0309018) q[2];
sx q[2];
rz(-0.008358566) q[2];
rz(-2.9110294) q[3];
sx q[3];
rz(-1.2059261) q[3];
sx q[3];
rz(1.4768538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01345988) q[0];
sx q[0];
rz(-0.020543329) q[0];
sx q[0];
rz(1.6101366) q[0];
rz(1.6186591) q[1];
sx q[1];
rz(-1.5459272) q[1];
sx q[1];
rz(-1.5628372) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3561479) q[0];
sx q[0];
rz(-1.4250989) q[0];
sx q[0];
rz(2.5375318) q[0];
rz(-pi) q[1];
rz(2.7885116) q[2];
sx q[2];
rz(-0.39829474) q[2];
sx q[2];
rz(0.19206599) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2819689) q[1];
sx q[1];
rz(-1.109237) q[1];
sx q[1];
rz(-2.0412316) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39852972) q[3];
sx q[3];
rz(-0.70793286) q[3];
sx q[3];
rz(1.7205451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3030887) q[2];
sx q[2];
rz(-1.9946626) q[2];
sx q[2];
rz(-3.1040891) q[2];
rz(2.904902) q[3];
sx q[3];
rz(-0.37875566) q[3];
sx q[3];
rz(1.8481351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26838747) q[0];
sx q[0];
rz(-0.78998843) q[0];
sx q[0];
rz(1.0733806) q[0];
rz(2.7997596) q[1];
sx q[1];
rz(-0.57917246) q[1];
sx q[1];
rz(-2.5198708) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2916002) q[0];
sx q[0];
rz(-1.1860285) q[0];
sx q[0];
rz(1.2906533) q[0];
x q[1];
rz(-2.9831977) q[2];
sx q[2];
rz(-0.6662874) q[2];
sx q[2];
rz(1.0384384) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3860064) q[1];
sx q[1];
rz(-1.2666128) q[1];
sx q[1];
rz(-2.6261397) q[1];
rz(2.5863566) q[3];
sx q[3];
rz(-0.16517565) q[3];
sx q[3];
rz(3.0495838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6217893) q[2];
sx q[2];
rz(-2.9534464) q[2];
sx q[2];
rz(0.56358799) q[2];
rz(1.311519) q[3];
sx q[3];
rz(-1.2389641) q[3];
sx q[3];
rz(-0.81609503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1821063) q[0];
sx q[0];
rz(-0.86903787) q[0];
sx q[0];
rz(-0.8748138) q[0];
rz(-1.733571) q[1];
sx q[1];
rz(-0.66649109) q[1];
sx q[1];
rz(2.5061238) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1615026) q[0];
sx q[0];
rz(-2.3005565) q[0];
sx q[0];
rz(-1.1671216) q[0];
x q[1];
rz(-2.8671625) q[2];
sx q[2];
rz(-1.1688902) q[2];
sx q[2];
rz(-2.2257471) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8567896) q[1];
sx q[1];
rz(-1.9724047) q[1];
sx q[1];
rz(1.0076341) q[1];
x q[2];
rz(2.1346774) q[3];
sx q[3];
rz(-2.0479408) q[3];
sx q[3];
rz(1.3729707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94329876) q[2];
sx q[2];
rz(-1.0903) q[2];
sx q[2];
rz(-0.80149209) q[2];
rz(1.9258026) q[3];
sx q[3];
rz(-2.2299168) q[3];
sx q[3];
rz(0.041778684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6112678) q[0];
sx q[0];
rz(-1.4401191) q[0];
sx q[0];
rz(-2.8339207) q[0];
rz(1.2904185) q[1];
sx q[1];
rz(-0.94480521) q[1];
sx q[1];
rz(0.31029846) q[1];
rz(1.7062832) q[2];
sx q[2];
rz(-1.3904962) q[2];
sx q[2];
rz(2.4615859) q[2];
rz(1.4361868) q[3];
sx q[3];
rz(-0.91840875) q[3];
sx q[3];
rz(2.8847532) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
