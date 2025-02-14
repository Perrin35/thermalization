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
rz(0.89927468) q[1];
sx q[1];
rz(-0.86714309) q[1];
sx q[1];
rz(1.4049621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1464935) q[0];
sx q[0];
rz(-1.4425264) q[0];
sx q[0];
rz(1.6544106) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5843732) q[2];
sx q[2];
rz(-1.4247334) q[2];
sx q[2];
rz(-2.5596325) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6331455) q[1];
sx q[1];
rz(-2.0237013) q[1];
sx q[1];
rz(-0.36302723) q[1];
x q[2];
rz(-1.8990535) q[3];
sx q[3];
rz(-2.0363288) q[3];
sx q[3];
rz(1.6159323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0286502) q[2];
sx q[2];
rz(-2.3990227) q[2];
sx q[2];
rz(0.58153018) q[2];
rz(0.90315008) q[3];
sx q[3];
rz(-1.6553469) q[3];
sx q[3];
rz(-2.7744897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1309995) q[0];
sx q[0];
rz(-1.4073263) q[0];
sx q[0];
rz(1.1372239) q[0];
rz(-1.8164903) q[1];
sx q[1];
rz(-2.4749327) q[1];
sx q[1];
rz(-2.4773662) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8695852) q[0];
sx q[0];
rz(-1.39886) q[0];
sx q[0];
rz(-0.26845064) q[0];
rz(-3.0587836) q[2];
sx q[2];
rz(-0.84575221) q[2];
sx q[2];
rz(0.74527455) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8926516) q[1];
sx q[1];
rz(-0.59160691) q[1];
sx q[1];
rz(-0.68343648) q[1];
rz(1.799753) q[3];
sx q[3];
rz(-1.6456283) q[3];
sx q[3];
rz(-0.018939806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5372411) q[2];
sx q[2];
rz(-1.9862572) q[2];
sx q[2];
rz(1.7604766) q[2];
rz(0.85754496) q[3];
sx q[3];
rz(-2.2159135) q[3];
sx q[3];
rz(-0.05923567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.183476) q[0];
sx q[0];
rz(-0.97267946) q[0];
sx q[0];
rz(2.3532975) q[0];
rz(2.6745785) q[1];
sx q[1];
rz(-0.6643509) q[1];
sx q[1];
rz(-0.83795396) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70650202) q[0];
sx q[0];
rz(-1.52924) q[0];
sx q[0];
rz(2.6986319) q[0];
x q[1];
rz(1.6347209) q[2];
sx q[2];
rz(-2.7294559) q[2];
sx q[2];
rz(-1.2391547) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0744255) q[1];
sx q[1];
rz(-0.29109368) q[1];
sx q[1];
rz(-1.8175248) q[1];
rz(-2.9492018) q[3];
sx q[3];
rz(-2.8225401) q[3];
sx q[3];
rz(-2.8448679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.18980846) q[2];
sx q[2];
rz(-0.34056792) q[2];
sx q[2];
rz(1.6050485) q[2];
rz(0.32478452) q[3];
sx q[3];
rz(-1.6580509) q[3];
sx q[3];
rz(1.270208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
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
rz(-1.6390653) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3887848) q[0];
sx q[0];
rz(-0.54600096) q[0];
sx q[0];
rz(-2.5939119) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45452228) q[2];
sx q[2];
rz(-1.3873867) q[2];
sx q[2];
rz(1.3647788) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6671424) q[1];
sx q[1];
rz(-2.7908142) q[1];
sx q[1];
rz(-0.79255484) q[1];
rz(0.32716515) q[3];
sx q[3];
rz(-1.7594271) q[3];
sx q[3];
rz(1.908055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23102418) q[2];
sx q[2];
rz(-0.39414057) q[2];
sx q[2];
rz(2.2310889) q[2];
rz(0.85601151) q[3];
sx q[3];
rz(-1.7249707) q[3];
sx q[3];
rz(-3.0205309) q[3];
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
rz(pi/2) q[3];
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
rz(-1.1055792) q[0];
sx q[0];
rz(-1.3904904) q[0];
sx q[0];
rz(-0.068583071) q[0];
rz(-0.71906459) q[1];
sx q[1];
rz(-1.9870575) q[1];
sx q[1];
rz(-2.7209435) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7375355) q[0];
sx q[0];
rz(-2.1348663) q[0];
sx q[0];
rz(0.89667908) q[0];
rz(-pi) q[1];
rz(0.44488944) q[2];
sx q[2];
rz(-0.25212461) q[2];
sx q[2];
rz(2.0201403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2202363) q[1];
sx q[1];
rz(-1.6097798) q[1];
sx q[1];
rz(-2.7552752) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4937872) q[3];
sx q[3];
rz(-1.2820377) q[3];
sx q[3];
rz(-1.6132068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61578304) q[2];
sx q[2];
rz(-2.1685648) q[2];
sx q[2];
rz(0.77965492) q[2];
rz(3.0443794) q[3];
sx q[3];
rz(-1.8153056) q[3];
sx q[3];
rz(1.1965082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.6355316) q[0];
sx q[0];
rz(-0.66672915) q[0];
sx q[0];
rz(2.7147103) q[0];
rz(-1.4510669) q[1];
sx q[1];
rz(-2.0207113) q[1];
sx q[1];
rz(-2.5371187) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0357416) q[0];
sx q[0];
rz(-1.5352016) q[0];
sx q[0];
rz(-1.2491519) q[0];
rz(2.8858917) q[2];
sx q[2];
rz(-1.2947645) q[2];
sx q[2];
rz(-2.7118341) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6410445) q[1];
sx q[1];
rz(-1.637008) q[1];
sx q[1];
rz(-1.3891298) q[1];
rz(-pi) q[2];
rz(-1.3434308) q[3];
sx q[3];
rz(-2.1220088) q[3];
sx q[3];
rz(0.77533305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.84588593) q[2];
sx q[2];
rz(-2.2658331) q[2];
sx q[2];
rz(2.8952944) q[2];
rz(-2.1379499) q[3];
sx q[3];
rz(-1.6946038) q[3];
sx q[3];
rz(-3.057737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87171444) q[0];
sx q[0];
rz(-0.074967472) q[0];
sx q[0];
rz(-0.20498928) q[0];
rz(0.21944731) q[1];
sx q[1];
rz(-1.4402025) q[1];
sx q[1];
rz(0.27935371) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.646307) q[0];
sx q[0];
rz(-1.6834462) q[0];
sx q[0];
rz(-1.1239284) q[0];
x q[1];
rz(-3.0952697) q[2];
sx q[2];
rz(-1.319482) q[2];
sx q[2];
rz(0.77071079) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0079035) q[1];
sx q[1];
rz(-0.55396307) q[1];
sx q[1];
rz(0.10793198) q[1];
x q[2];
rz(2.1198307) q[3];
sx q[3];
rz(-0.543299) q[3];
sx q[3];
rz(-0.70453139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67639914) q[2];
sx q[2];
rz(-0.78130829) q[2];
sx q[2];
rz(0.87731963) q[2];
rz(1.8026132) q[3];
sx q[3];
rz(-2.3586912) q[3];
sx q[3];
rz(1.7507929) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1124697) q[0];
sx q[0];
rz(-0.59337297) q[0];
sx q[0];
rz(0.33475885) q[0];
rz(-0.33084694) q[1];
sx q[1];
rz(-2.0956764) q[1];
sx q[1];
rz(0.61029339) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7772737) q[0];
sx q[0];
rz(-2.7734682) q[0];
sx q[0];
rz(-2.0859531) q[0];
rz(-1.1534821) q[2];
sx q[2];
rz(-1.8886107) q[2];
sx q[2];
rz(1.870188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.6107294) q[1];
sx q[1];
rz(-1.7708774) q[1];
sx q[1];
rz(-2.6310573) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0078681) q[3];
sx q[3];
rz(-0.96157295) q[3];
sx q[3];
rz(-1.5214868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.200013) q[2];
sx q[2];
rz(-1.0369077) q[2];
sx q[2];
rz(-1.3882136) q[2];
rz(2.9774104) q[3];
sx q[3];
rz(-2.1037585) q[3];
sx q[3];
rz(-2.8431456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91069094) q[0];
sx q[0];
rz(-2.0257484) q[0];
sx q[0];
rz(-0.13634613) q[0];
rz(0.048010437) q[1];
sx q[1];
rz(-1.0142356) q[1];
sx q[1];
rz(-1.8528574) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.890402) q[0];
sx q[0];
rz(-1.2689166) q[0];
sx q[0];
rz(2.9142778) q[0];
rz(0.23581712) q[2];
sx q[2];
rz(-1.1305446) q[2];
sx q[2];
rz(-2.4503675) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9077085) q[1];
sx q[1];
rz(-1.4310992) q[1];
sx q[1];
rz(-0.25821347) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7461122) q[3];
sx q[3];
rz(-2.2321475) q[3];
sx q[3];
rz(-1.1547853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5137704) q[2];
sx q[2];
rz(-3.1352391) q[2];
sx q[2];
rz(-2.9610236) q[2];
rz(2.0020961) q[3];
sx q[3];
rz(-1.5151016) q[3];
sx q[3];
rz(-0.12038055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1861495) q[0];
sx q[0];
rz(-2.1685765) q[0];
sx q[0];
rz(-1.0261616) q[0];
rz(1.8852662) q[1];
sx q[1];
rz(-1.8812814) q[1];
sx q[1];
rz(-0.06591448) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33082066) q[0];
sx q[0];
rz(-2.2123824) q[0];
sx q[0];
rz(0.75831428) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49123571) q[2];
sx q[2];
rz(-0.35235534) q[2];
sx q[2];
rz(-1.6879326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5951426) q[1];
sx q[1];
rz(-1.3347581) q[1];
sx q[1];
rz(-2.5883893) q[1];
x q[2];
rz(-1.3168748) q[3];
sx q[3];
rz(-0.99813509) q[3];
sx q[3];
rz(0.66324011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5790448) q[2];
sx q[2];
rz(-2.0936091) q[2];
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
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3271306) q[0];
sx q[0];
rz(-0.7939864) q[0];
sx q[0];
rz(-3.0891147) q[0];
rz(1.3606701) q[1];
sx q[1];
rz(-2.4130029) q[1];
sx q[1];
rz(-1.3389814) q[1];
rz(-0.95407964) q[2];
sx q[2];
rz(-0.5250684) q[2];
sx q[2];
rz(-2.9697408) q[2];
rz(0.94447847) q[3];
sx q[3];
rz(-1.2068268) q[3];
sx q[3];
rz(-2.5218388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
