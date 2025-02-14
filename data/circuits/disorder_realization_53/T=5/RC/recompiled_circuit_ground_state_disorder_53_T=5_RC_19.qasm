OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5767515) q[0];
sx q[0];
rz(-0.54902005) q[0];
sx q[0];
rz(-2.9615871) q[0];
rz(4.0408673) q[1];
sx q[1];
rz(4.0087357) q[1];
sx q[1];
rz(8.0198159) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7265018) q[0];
sx q[0];
rz(-0.15299061) q[0];
sx q[0];
rz(0.57463519) q[0];
rz(-pi) q[1];
rz(1.7424217) q[2];
sx q[2];
rz(-2.1213946) q[2];
sx q[2];
rz(-1.0792749) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3480543) q[1];
sx q[1];
rz(-2.5691433) q[1];
sx q[1];
rz(-0.94041108) q[1];
rz(2.5710315) q[3];
sx q[3];
rz(-0.56260144) q[3];
sx q[3];
rz(-2.1747052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11294242) q[2];
sx q[2];
rz(-0.74256998) q[2];
sx q[2];
rz(-2.5600625) q[2];
rz(0.90315008) q[3];
sx q[3];
rz(-1.6553469) q[3];
sx q[3];
rz(0.36710292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1309995) q[0];
sx q[0];
rz(-1.7342664) q[0];
sx q[0];
rz(2.0043688) q[0];
rz(1.8164903) q[1];
sx q[1];
rz(-2.4749327) q[1];
sx q[1];
rz(2.4773662) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8839845) q[0];
sx q[0];
rz(-0.31766787) q[0];
sx q[0];
rz(-0.57967107) q[0];
rz(0.84404805) q[2];
sx q[2];
rz(-1.6327452) q[2];
sx q[2];
rz(0.88050851) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6146734) q[1];
sx q[1];
rz(-2.0179948) q[1];
sx q[1];
rz(1.1695444) q[1];
rz(-pi) q[2];
rz(-3.0647633) q[3];
sx q[3];
rz(-1.3424918) q[3];
sx q[3];
rz(1.5692776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5372411) q[2];
sx q[2];
rz(-1.9862572) q[2];
sx q[2];
rz(-1.7604766) q[2];
rz(0.85754496) q[3];
sx q[3];
rz(-2.2159135) q[3];
sx q[3];
rz(-0.05923567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9581167) q[0];
sx q[0];
rz(-2.1689132) q[0];
sx q[0];
rz(-2.3532975) q[0];
rz(-2.6745785) q[1];
sx q[1];
rz(-0.6643509) q[1];
sx q[1];
rz(-2.3036387) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77694511) q[0];
sx q[0];
rz(-2.6968156) q[0];
sx q[0];
rz(0.096707896) q[0];
x q[1];
rz(-1.6347209) q[2];
sx q[2];
rz(-2.7294559) q[2];
sx q[2];
rz(-1.9024379) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3242692) q[1];
sx q[1];
rz(-1.2887635) q[1];
sx q[1];
rz(-3.0685496) q[1];
rz(2.8280452) q[3];
sx q[3];
rz(-1.6308074) q[3];
sx q[3];
rz(1.0911694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9517842) q[2];
sx q[2];
rz(-2.8010247) q[2];
sx q[2];
rz(1.5365441) q[2];
rz(2.8168081) q[3];
sx q[3];
rz(-1.4835417) q[3];
sx q[3];
rz(-1.8713846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.58295163) q[0];
sx q[0];
rz(-2.1138209) q[0];
sx q[0];
rz(-0.047274832) q[0];
rz(1.524823) q[1];
sx q[1];
rz(-0.44175092) q[1];
sx q[1];
rz(1.5025274) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3887848) q[0];
sx q[0];
rz(-0.54600096) q[0];
sx q[0];
rz(-2.5939119) q[0];
x q[1];
rz(1.7743913) q[2];
sx q[2];
rz(-2.0171391) q[2];
sx q[2];
rz(2.8466895) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33519519) q[1];
sx q[1];
rz(-1.8180221) q[1];
sx q[1];
rz(2.8901495) q[1];
x q[2];
rz(2.8144275) q[3];
sx q[3];
rz(-1.7594271) q[3];
sx q[3];
rz(-1.908055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23102418) q[2];
sx q[2];
rz(-0.39414057) q[2];
sx q[2];
rz(-0.91050371) q[2];
rz(-0.85601151) q[3];
sx q[3];
rz(-1.7249707) q[3];
sx q[3];
rz(-0.12106171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1055792) q[0];
sx q[0];
rz(-1.3904904) q[0];
sx q[0];
rz(3.0730096) q[0];
rz(2.4225281) q[1];
sx q[1];
rz(-1.1545352) q[1];
sx q[1];
rz(2.7209435) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42297573) q[0];
sx q[0];
rz(-2.291922) q[0];
sx q[0];
rz(2.3628985) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44488944) q[2];
sx q[2];
rz(-0.25212461) q[2];
sx q[2];
rz(-2.0201403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8068841) q[1];
sx q[1];
rz(-1.9568048) q[1];
sx q[1];
rz(-1.5287148) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4937872) q[3];
sx q[3];
rz(-1.8595549) q[3];
sx q[3];
rz(-1.6132068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5258096) q[2];
sx q[2];
rz(-2.1685648) q[2];
sx q[2];
rz(2.3619377) q[2];
rz(-0.097213216) q[3];
sx q[3];
rz(-1.8153056) q[3];
sx q[3];
rz(-1.9450845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6355316) q[0];
sx q[0];
rz(-2.4748635) q[0];
sx q[0];
rz(-2.7147103) q[0];
rz(-1.4510669) q[1];
sx q[1];
rz(-1.1208813) q[1];
sx q[1];
rz(2.5371187) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35854912) q[0];
sx q[0];
rz(-0.32354) q[0];
sx q[0];
rz(-1.4586252) q[0];
rz(-pi) q[1];
rz(-1.2859743) q[2];
sx q[2];
rz(-1.3249791) q[2];
sx q[2];
rz(2.071683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8661202) q[1];
sx q[1];
rz(-0.19323142) q[1];
sx q[1];
rz(1.92255) q[1];
x q[2];
rz(2.5786798) q[3];
sx q[3];
rz(-1.7640224) q[3];
sx q[3];
rz(2.4667127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.84588593) q[2];
sx q[2];
rz(-2.2658331) q[2];
sx q[2];
rz(0.24629822) q[2];
rz(1.0036428) q[3];
sx q[3];
rz(-1.6946038) q[3];
sx q[3];
rz(0.083855696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.87171444) q[0];
sx q[0];
rz(-3.0666252) q[0];
sx q[0];
rz(0.20498928) q[0];
rz(-2.9221453) q[1];
sx q[1];
rz(-1.7013902) q[1];
sx q[1];
rz(-0.27935371) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.646307) q[0];
sx q[0];
rz(-1.6834462) q[0];
sx q[0];
rz(-1.1239284) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.749239) q[2];
sx q[2];
rz(-2.8861336) q[2];
sx q[2];
rz(-0.58641543) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6705888) q[1];
sx q[1];
rz(-1.6274954) q[1];
sx q[1];
rz(-0.55135552) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0950915) q[3];
sx q[3];
rz(-1.8439652) q[3];
sx q[3];
rz(-2.7577443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4651935) q[2];
sx q[2];
rz(-2.3602844) q[2];
sx q[2];
rz(-2.264273) q[2];
rz(-1.3389795) q[3];
sx q[3];
rz(-0.78290144) q[3];
sx q[3];
rz(-1.7507929) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1124697) q[0];
sx q[0];
rz(-0.59337297) q[0];
sx q[0];
rz(2.8068338) q[0];
rz(-2.8107457) q[1];
sx q[1];
rz(-2.0956764) q[1];
sx q[1];
rz(2.5312993) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2318678) q[0];
sx q[0];
rz(-1.2522765) q[0];
sx q[0];
rz(-0.18778778) q[0];
x q[1];
rz(1.1534821) q[2];
sx q[2];
rz(-1.252982) q[2];
sx q[2];
rz(-1.2714047) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0706705) q[1];
sx q[1];
rz(-1.0713995) q[1];
sx q[1];
rz(-1.3424177) q[1];
x q[2];
rz(-0.957358) q[3];
sx q[3];
rz(-1.4612373) q[3];
sx q[3];
rz(0.027519634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.200013) q[2];
sx q[2];
rz(-1.0369077) q[2];
sx q[2];
rz(1.3882136) q[2];
rz(0.16418223) q[3];
sx q[3];
rz(-1.0378342) q[3];
sx q[3];
rz(-2.8431456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2309017) q[0];
sx q[0];
rz(-1.1158442) q[0];
sx q[0];
rz(-0.13634613) q[0];
rz(3.0935822) q[1];
sx q[1];
rz(-1.0142356) q[1];
sx q[1];
rz(-1.2887352) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25094095) q[0];
sx q[0];
rz(-1.3539292) q[0];
sx q[0];
rz(-1.8801539) q[0];
x q[1];
rz(1.1103915) q[2];
sx q[2];
rz(-2.6458323) q[2];
sx q[2];
rz(2.9637703) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37367108) q[1];
sx q[1];
rz(-1.8264378) q[1];
sx q[1];
rz(1.7152183) q[1];
rz(0.8701635) q[3];
sx q[3];
rz(-1.2618801) q[3];
sx q[3];
rz(-0.66701057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62782225) q[2];
sx q[2];
rz(-0.0063535293) q[2];
sx q[2];
rz(-2.9610236) q[2];
rz(-2.0020961) q[3];
sx q[3];
rz(-1.5151016) q[3];
sx q[3];
rz(0.12038055) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1861495) q[0];
sx q[0];
rz(-2.1685765) q[0];
sx q[0];
rz(-1.0261616) q[0];
rz(-1.8852662) q[1];
sx q[1];
rz(-1.8812814) q[1];
sx q[1];
rz(0.06591448) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72425264) q[0];
sx q[0];
rz(-2.1542962) q[0];
sx q[0];
rz(2.3704607) q[0];
rz(-pi) q[1];
rz(0.49123571) q[2];
sx q[2];
rz(-0.35235534) q[2];
sx q[2];
rz(1.6879326) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0225214) q[1];
sx q[1];
rz(-2.1069659) q[1];
sx q[1];
rz(-1.8462935) q[1];
x q[2];
rz(-2.7700636) q[3];
sx q[3];
rz(-2.5209628) q[3];
sx q[3];
rz(0.21658235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5790448) q[2];
sx q[2];
rz(-1.0479835) q[2];
sx q[2];
rz(-0.5272131) q[2];
rz(2.9910917) q[3];
sx q[3];
rz(-2.7121057) q[3];
sx q[3];
rz(2.785717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81446205) q[0];
sx q[0];
rz(-2.3476063) q[0];
sx q[0];
rz(0.052477947) q[0];
rz(1.3606701) q[1];
sx q[1];
rz(-2.4130029) q[1];
sx q[1];
rz(-1.3389814) q[1];
rz(2.0122779) q[2];
sx q[2];
rz(-1.2766576) q[2];
sx q[2];
rz(1.192391) q[2];
rz(-0.99450022) q[3];
sx q[3];
rz(-2.4296843) q[3];
sx q[3];
rz(-1.4083023) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
