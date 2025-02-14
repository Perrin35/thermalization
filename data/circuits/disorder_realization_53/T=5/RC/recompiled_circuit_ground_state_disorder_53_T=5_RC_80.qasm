OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.56484115) q[0];
sx q[0];
rz(-2.5925726) q[0];
sx q[0];
rz(-0.18000552) q[0];
rz(0.89927468) q[1];
sx q[1];
rz(-0.86714309) q[1];
sx q[1];
rz(1.4049621) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9950992) q[0];
sx q[0];
rz(-1.6990663) q[0];
sx q[0];
rz(1.6544106) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27132224) q[2];
sx q[2];
rz(-2.567511) q[2];
sx q[2];
rz(0.75936356) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3480543) q[1];
sx q[1];
rz(-0.57244937) q[1];
sx q[1];
rz(2.2011816) q[1];
rz(0.48790055) q[3];
sx q[3];
rz(-1.8630233) q[3];
sx q[3];
rz(3.0350181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0286502) q[2];
sx q[2];
rz(-0.74256998) q[2];
sx q[2];
rz(-0.58153018) q[2];
rz(0.90315008) q[3];
sx q[3];
rz(-1.6553469) q[3];
sx q[3];
rz(-2.7744897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0105932) q[0];
sx q[0];
rz(-1.4073263) q[0];
sx q[0];
rz(-2.0043688) q[0];
rz(1.8164903) q[1];
sx q[1];
rz(-0.66665998) q[1];
sx q[1];
rz(0.66422647) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8839845) q[0];
sx q[0];
rz(-0.31766787) q[0];
sx q[0];
rz(-0.57967107) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2975446) q[2];
sx q[2];
rz(-1.6327452) q[2];
sx q[2];
rz(-2.2610841) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8926516) q[1];
sx q[1];
rz(-0.59160691) q[1];
sx q[1];
rz(2.4581562) q[1];
rz(-pi) q[2];
x q[2];
rz(1.799753) q[3];
sx q[3];
rz(-1.6456283) q[3];
sx q[3];
rz(3.1226528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.60435158) q[2];
sx q[2];
rz(-1.1553355) q[2];
sx q[2];
rz(1.381116) q[2];
rz(0.85754496) q[3];
sx q[3];
rz(-0.92567912) q[3];
sx q[3];
rz(-3.082357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9581167) q[0];
sx q[0];
rz(-2.1689132) q[0];
sx q[0];
rz(-0.78829515) q[0];
rz(0.46701416) q[1];
sx q[1];
rz(-2.4772418) q[1];
sx q[1];
rz(-0.83795396) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4350906) q[0];
sx q[0];
rz(-1.6123527) q[0];
sx q[0];
rz(-2.6986319) q[0];
rz(0.027919876) q[2];
sx q[2];
rz(-1.1595524) q[2];
sx q[2];
rz(1.9721861) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8173234) q[1];
sx q[1];
rz(-1.8528291) q[1];
sx q[1];
rz(-3.0685496) q[1];
rz(-pi) q[2];
rz(0.31354745) q[3];
sx q[3];
rz(-1.6308074) q[3];
sx q[3];
rz(-1.0911694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9517842) q[2];
sx q[2];
rz(-0.34056792) q[2];
sx q[2];
rz(-1.5365441) q[2];
rz(0.32478452) q[3];
sx q[3];
rz(-1.4835417) q[3];
sx q[3];
rz(-1.270208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58295163) q[0];
sx q[0];
rz(-2.1138209) q[0];
sx q[0];
rz(0.047274832) q[0];
rz(-1.524823) q[1];
sx q[1];
rz(-2.6998417) q[1];
sx q[1];
rz(1.5025274) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3887848) q[0];
sx q[0];
rz(-0.54600096) q[0];
sx q[0];
rz(2.5939119) q[0];
x q[1];
rz(-1.3672013) q[2];
sx q[2];
rz(-2.0171391) q[2];
sx q[2];
rz(2.8466895) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8063975) q[1];
sx q[1];
rz(-1.3235705) q[1];
sx q[1];
rz(0.25144318) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8144275) q[3];
sx q[3];
rz(-1.7594271) q[3];
sx q[3];
rz(1.908055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9105685) q[2];
sx q[2];
rz(-2.7474521) q[2];
sx q[2];
rz(2.2310889) q[2];
rz(-2.2855811) q[3];
sx q[3];
rz(-1.7249707) q[3];
sx q[3];
rz(-3.0205309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.0360134) q[0];
sx q[0];
rz(-1.7511022) q[0];
sx q[0];
rz(-3.0730096) q[0];
rz(-2.4225281) q[1];
sx q[1];
rz(-1.1545352) q[1];
sx q[1];
rz(-2.7209435) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7186169) q[0];
sx q[0];
rz(-0.84967063) q[0];
sx q[0];
rz(0.77869418) q[0];
x q[1];
rz(2.6967032) q[2];
sx q[2];
rz(-2.889468) q[2];
sx q[2];
rz(2.0201403) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2202363) q[1];
sx q[1];
rz(-1.6097798) q[1];
sx q[1];
rz(2.7552752) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64780541) q[3];
sx q[3];
rz(-1.2820377) q[3];
sx q[3];
rz(-1.5283858) q[3];
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
rz(0.097213216) q[3];
sx q[3];
rz(-1.8153056) q[3];
sx q[3];
rz(-1.1965082) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6355316) q[0];
sx q[0];
rz(-0.66672915) q[0];
sx q[0];
rz(0.42688236) q[0];
rz(1.6905258) q[1];
sx q[1];
rz(-2.0207113) q[1];
sx q[1];
rz(-2.5371187) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7830435) q[0];
sx q[0];
rz(-0.32354) q[0];
sx q[0];
rz(-1.6829674) q[0];
rz(-0.25570095) q[2];
sx q[2];
rz(-1.2947645) q[2];
sx q[2];
rz(-2.7118341) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6410445) q[1];
sx q[1];
rz(-1.637008) q[1];
sx q[1];
rz(-1.3891298) q[1];
rz(-pi) q[2];
rz(-0.56291286) q[3];
sx q[3];
rz(-1.3775702) q[3];
sx q[3];
rz(-2.4667127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2957067) q[2];
sx q[2];
rz(-2.2658331) q[2];
sx q[2];
rz(-2.8952944) q[2];
rz(2.1379499) q[3];
sx q[3];
rz(-1.4469888) q[3];
sx q[3];
rz(0.083855696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2698782) q[0];
sx q[0];
rz(-0.074967472) q[0];
sx q[0];
rz(2.9366034) q[0];
rz(0.21944731) q[1];
sx q[1];
rz(-1.4402025) q[1];
sx q[1];
rz(0.27935371) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0216952) q[0];
sx q[0];
rz(-1.1269635) q[0];
sx q[0];
rz(0.12479513) q[0];
x q[1];
rz(1.8223693) q[2];
sx q[2];
rz(-1.5259296) q[2];
sx q[2];
rz(0.81161273) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6705888) q[1];
sx q[1];
rz(-1.5140972) q[1];
sx q[1];
rz(2.5902371) q[1];
rz(-pi) q[2];
x q[2];
rz(1.021762) q[3];
sx q[3];
rz(-2.5982937) q[3];
sx q[3];
rz(-0.70453139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67639914) q[2];
sx q[2];
rz(-0.78130829) q[2];
sx q[2];
rz(-0.87731963) q[2];
rz(1.8026132) q[3];
sx q[3];
rz(-2.3586912) q[3];
sx q[3];
rz(-1.3907998) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.029123) q[0];
sx q[0];
rz(-2.5482197) q[0];
sx q[0];
rz(2.8068338) q[0];
rz(2.8107457) q[1];
sx q[1];
rz(-2.0956764) q[1];
sx q[1];
rz(-2.5312993) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7205104) q[0];
sx q[0];
rz(-1.7490343) q[0];
sx q[0];
rz(-1.2469638) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3454204) q[2];
sx q[2];
rz(-1.1755786) q[2];
sx q[2];
rz(0.43705979) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0706705) q[1];
sx q[1];
rz(-1.0713995) q[1];
sx q[1];
rz(-1.7991749) q[1];
x q[2];
rz(-0.13372453) q[3];
sx q[3];
rz(-2.1800197) q[3];
sx q[3];
rz(-1.5214868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.94157964) q[2];
sx q[2];
rz(-2.104685) q[2];
sx q[2];
rz(1.3882136) q[2];
rz(2.9774104) q[3];
sx q[3];
rz(-1.0378342) q[3];
sx q[3];
rz(2.8431456) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2309017) q[0];
sx q[0];
rz(-1.1158442) q[0];
sx q[0];
rz(-3.0052465) q[0];
rz(-3.0935822) q[1];
sx q[1];
rz(-1.0142356) q[1];
sx q[1];
rz(-1.8528574) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2511906) q[0];
sx q[0];
rz(-1.8726761) q[0];
sx q[0];
rz(-0.22731486) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23581712) q[2];
sx q[2];
rz(-2.011048) q[2];
sx q[2];
rz(0.69122512) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37367108) q[1];
sx q[1];
rz(-1.3151549) q[1];
sx q[1];
rz(-1.7152183) q[1];
rz(-0.39548042) q[3];
sx q[3];
rz(-0.90944511) q[3];
sx q[3];
rz(-1.1547853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62782225) q[2];
sx q[2];
rz(-0.0063535293) q[2];
sx q[2];
rz(-2.9610236) q[2];
rz(-1.1394966) q[3];
sx q[3];
rz(-1.6264911) q[3];
sx q[3];
rz(0.12038055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9554431) q[0];
sx q[0];
rz(-0.97301617) q[0];
sx q[0];
rz(2.1154311) q[0];
rz(1.8852662) q[1];
sx q[1];
rz(-1.8812814) q[1];
sx q[1];
rz(-0.06591448) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.41734) q[0];
sx q[0];
rz(-2.1542962) q[0];
sx q[0];
rz(-2.3704607) q[0];
rz(-1.3990551) q[2];
sx q[2];
rz(-1.8800003) q[2];
sx q[2];
rz(-1.9717252) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7552398) q[1];
sx q[1];
rz(-2.5449979) q[1];
sx q[1];
rz(-0.42930023) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58760175) q[3];
sx q[3];
rz(-1.7835225) q[3];
sx q[3];
rz(-1.0472681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5790448) q[2];
sx q[2];
rz(-1.0479835) q[2];
sx q[2];
rz(2.6143796) q[2];
rz(-0.15050091) q[3];
sx q[3];
rz(-2.7121057) q[3];
sx q[3];
rz(-0.35587564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-2.3271306) q[0];
sx q[0];
rz(-0.7939864) q[0];
sx q[0];
rz(-3.0891147) q[0];
rz(1.7809226) q[1];
sx q[1];
rz(-0.7285898) q[1];
sx q[1];
rz(1.8026112) q[1];
rz(-2.187513) q[2];
sx q[2];
rz(-2.6165243) q[2];
sx q[2];
rz(0.17185186) q[2];
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
