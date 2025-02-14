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
rz(3.6906127) q[0];
sx q[0];
rz(9.2447724) q[0];
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
rz(1.9950992) q[0];
sx q[0];
rz(-1.6990663) q[0];
sx q[0];
rz(-1.6544106) q[0];
x q[1];
rz(1.7424217) q[2];
sx q[2];
rz(-2.1213946) q[2];
sx q[2];
rz(2.0623178) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22706341) q[1];
sx q[1];
rz(-1.8957912) q[1];
sx q[1];
rz(-2.0507647) q[1];
rz(2.6536921) q[3];
sx q[3];
rz(-1.8630233) q[3];
sx q[3];
rz(-3.0350181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.11294242) q[2];
sx q[2];
rz(-2.3990227) q[2];
sx q[2];
rz(-2.5600625) q[2];
rz(2.2384426) q[3];
sx q[3];
rz(-1.4862458) q[3];
sx q[3];
rz(-2.7744897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0105932) q[0];
sx q[0];
rz(-1.7342664) q[0];
sx q[0];
rz(-1.1372239) q[0];
rz(-1.3251023) q[1];
sx q[1];
rz(-2.4749327) q[1];
sx q[1];
rz(2.4773662) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8695852) q[0];
sx q[0];
rz(-1.39886) q[0];
sx q[0];
rz(-0.26845064) q[0];
rz(2.2975446) q[2];
sx q[2];
rz(-1.6327452) q[2];
sx q[2];
rz(2.2610841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8926516) q[1];
sx q[1];
rz(-0.59160691) q[1];
sx q[1];
rz(0.68343648) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.076829314) q[3];
sx q[3];
rz(-1.3424918) q[3];
sx q[3];
rz(1.5723151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5372411) q[2];
sx q[2];
rz(-1.9862572) q[2];
sx q[2];
rz(-1.381116) q[2];
rz(2.2840477) q[3];
sx q[3];
rz(-0.92567912) q[3];
sx q[3];
rz(-0.05923567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9581167) q[0];
sx q[0];
rz(-2.1689132) q[0];
sx q[0];
rz(-0.78829515) q[0];
rz(2.6745785) q[1];
sx q[1];
rz(-2.4772418) q[1];
sx q[1];
rz(-2.3036387) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2575921) q[0];
sx q[0];
rz(-1.1282451) q[0];
sx q[0];
rz(1.5248067) q[0];
rz(1.6347209) q[2];
sx q[2];
rz(-2.7294559) q[2];
sx q[2];
rz(1.9024379) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0671671) q[1];
sx q[1];
rz(-0.29109368) q[1];
sx q[1];
rz(1.3240678) q[1];
x q[2];
rz(0.19239088) q[3];
sx q[3];
rz(-0.31905252) q[3];
sx q[3];
rz(-0.2967248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.18980846) q[2];
sx q[2];
rz(-0.34056792) q[2];
sx q[2];
rz(1.6050485) q[2];
rz(-0.32478452) q[3];
sx q[3];
rz(-1.6580509) q[3];
sx q[3];
rz(-1.270208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58295163) q[0];
sx q[0];
rz(-2.1138209) q[0];
sx q[0];
rz(0.047274832) q[0];
rz(-1.6167697) q[1];
sx q[1];
rz(-0.44175092) q[1];
sx q[1];
rz(-1.6390653) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0086347) q[0];
sx q[0];
rz(-2.0300976) q[0];
sx q[0];
rz(-1.8772222) q[0];
rz(-pi) q[1];
rz(-0.39975651) q[2];
sx q[2];
rz(-0.48772487) q[2];
sx q[2];
rz(2.9903811) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8432118) q[1];
sx q[1];
rz(-1.8144368) q[1];
sx q[1];
rz(1.8257105) q[1];
rz(-1.7697236) q[3];
sx q[3];
rz(-1.2496447) q[3];
sx q[3];
rz(-0.40080788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.23102418) q[2];
sx q[2];
rz(-0.39414057) q[2];
sx q[2];
rz(2.2310889) q[2];
rz(2.2855811) q[3];
sx q[3];
rz(-1.4166219) q[3];
sx q[3];
rz(0.12106171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-0.4206492) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7375355) q[0];
sx q[0];
rz(-2.1348663) q[0];
sx q[0];
rz(-0.89667908) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6967032) q[2];
sx q[2];
rz(-2.889468) q[2];
sx q[2];
rz(-2.0201403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9213564) q[1];
sx q[1];
rz(-1.6097798) q[1];
sx q[1];
rz(2.7552752) q[1];
x q[2];
rz(0.4574538) q[3];
sx q[3];
rz(-2.4408999) q[3];
sx q[3];
rz(2.7392741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5258096) q[2];
sx q[2];
rz(-0.97302786) q[2];
sx q[2];
rz(-2.3619377) q[2];
rz(-0.097213216) q[3];
sx q[3];
rz(-1.3262871) q[3];
sx q[3];
rz(1.9450845) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6355316) q[0];
sx q[0];
rz(-0.66672915) q[0];
sx q[0];
rz(0.42688236) q[0];
rz(1.4510669) q[1];
sx q[1];
rz(-1.1208813) q[1];
sx q[1];
rz(0.60447398) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6647897) q[0];
sx q[0];
rz(-1.8922297) q[0];
sx q[0];
rz(3.1040756) q[0];
rz(-0.25570095) q[2];
sx q[2];
rz(-1.8468282) q[2];
sx q[2];
rz(2.7118341) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.082401245) q[1];
sx q[1];
rz(-1.7520604) q[1];
sx q[1];
rz(-3.0742765) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7901561) q[3];
sx q[3];
rz(-2.5498423) q[3];
sx q[3];
rz(-1.9502673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84588593) q[2];
sx q[2];
rz(-2.2658331) q[2];
sx q[2];
rz(-0.24629822) q[2];
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
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87171444) q[0];
sx q[0];
rz(-3.0666252) q[0];
sx q[0];
rz(-0.20498928) q[0];
rz(-0.21944731) q[1];
sx q[1];
rz(-1.7013902) q[1];
sx q[1];
rz(0.27935371) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.646307) q[0];
sx q[0];
rz(-1.4581465) q[0];
sx q[0];
rz(-2.0176642) q[0];
rz(-pi) q[1];
rz(1.749239) q[2];
sx q[2];
rz(-0.25545909) q[2];
sx q[2];
rz(-0.58641543) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.134628) q[1];
sx q[1];
rz(-1.0204287) q[1];
sx q[1];
rz(-1.637332) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0465011) q[3];
sx q[3];
rz(-1.8439652) q[3];
sx q[3];
rz(-0.38384837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.67639914) q[2];
sx q[2];
rz(-0.78130829) q[2];
sx q[2];
rz(0.87731963) q[2];
rz(1.3389795) q[3];
sx q[3];
rz(-0.78290144) q[3];
sx q[3];
rz(-1.3907998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1124697) q[0];
sx q[0];
rz(-2.5482197) q[0];
sx q[0];
rz(0.33475885) q[0];
rz(-0.33084694) q[1];
sx q[1];
rz(-1.0459162) q[1];
sx q[1];
rz(-0.61029339) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36431896) q[0];
sx q[0];
rz(-0.36812447) q[0];
sx q[0];
rz(1.0556396) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3454204) q[2];
sx q[2];
rz(-1.9660141) q[2];
sx q[2];
rz(0.43705979) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5308633) q[1];
sx q[1];
rz(-1.3707152) q[1];
sx q[1];
rz(-0.51053534) q[1];
rz(-pi) q[2];
x q[2];
rz(0.957358) q[3];
sx q[3];
rz(-1.6803553) q[3];
sx q[3];
rz(-3.114073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.94157964) q[2];
sx q[2];
rz(-2.104685) q[2];
sx q[2];
rz(-1.3882136) q[2];
rz(0.16418223) q[3];
sx q[3];
rz(-1.0378342) q[3];
sx q[3];
rz(-2.8431456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
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
rz(-2.127357) q[1];
sx q[1];
rz(1.2887352) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8906517) q[0];
sx q[0];
rz(-1.3539292) q[0];
sx q[0];
rz(1.8801539) q[0];
rz(1.1103915) q[2];
sx q[2];
rz(-2.6458323) q[2];
sx q[2];
rz(-0.17782234) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7679216) q[1];
sx q[1];
rz(-1.3151549) q[1];
sx q[1];
rz(1.4263743) q[1];
rz(-0.39548042) q[3];
sx q[3];
rz(-2.2321475) q[3];
sx q[3];
rz(-1.9868074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62782225) q[2];
sx q[2];
rz(-3.1352391) q[2];
sx q[2];
rz(-0.18056907) q[2];
rz(-2.0020961) q[3];
sx q[3];
rz(-1.6264911) q[3];
sx q[3];
rz(3.0212121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9554431) q[0];
sx q[0];
rz(-0.97301617) q[0];
sx q[0];
rz(-2.1154311) q[0];
rz(-1.2563264) q[1];
sx q[1];
rz(-1.2603113) q[1];
sx q[1];
rz(-3.0756782) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.810772) q[0];
sx q[0];
rz(-2.2123824) q[0];
sx q[0];
rz(-2.3832784) q[0];
x q[1];
rz(2.6503569) q[2];
sx q[2];
rz(-0.35235534) q[2];
sx q[2];
rz(1.45366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5951426) q[1];
sx q[1];
rz(-1.3347581) q[1];
sx q[1];
rz(2.5883893) q[1];
rz(-pi) q[2];
rz(-2.7700636) q[3];
sx q[3];
rz(-0.62062988) q[3];
sx q[3];
rz(2.9250103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5790448) q[2];
sx q[2];
rz(-1.0479835) q[2];
sx q[2];
rz(0.5272131) q[2];
rz(0.15050091) q[3];
sx q[3];
rz(-2.7121057) q[3];
sx q[3];
rz(-2.785717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81446205) q[0];
sx q[0];
rz(-0.7939864) q[0];
sx q[0];
rz(-3.0891147) q[0];
rz(1.3606701) q[1];
sx q[1];
rz(-2.4130029) q[1];
sx q[1];
rz(-1.3389814) q[1];
rz(-2.8182975) q[2];
sx q[2];
rz(-1.9920762) q[2];
sx q[2];
rz(2.6270234) q[2];
rz(-0.94447847) q[3];
sx q[3];
rz(-1.9347659) q[3];
sx q[3];
rz(0.61975382) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
