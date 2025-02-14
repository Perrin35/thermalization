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
rz(-1.7366306) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1464935) q[0];
sx q[0];
rz(-1.4425264) q[0];
sx q[0];
rz(1.487182) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7424217) q[2];
sx q[2];
rz(-2.1213946) q[2];
sx q[2];
rz(2.0623178) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9145292) q[1];
sx q[1];
rz(-1.8957912) q[1];
sx q[1];
rz(-1.090828) q[1];
rz(-pi) q[2];
rz(-1.8990535) q[3];
sx q[3];
rz(-1.1052638) q[3];
sx q[3];
rz(1.5256603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11294242) q[2];
sx q[2];
rz(-0.74256998) q[2];
sx q[2];
rz(2.5600625) q[2];
rz(-0.90315008) q[3];
sx q[3];
rz(-1.6553469) q[3];
sx q[3];
rz(-0.36710292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0105932) q[0];
sx q[0];
rz(-1.7342664) q[0];
sx q[0];
rz(-1.1372239) q[0];
rz(1.3251023) q[1];
sx q[1];
rz(-2.4749327) q[1];
sx q[1];
rz(0.66422647) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25760812) q[0];
sx q[0];
rz(-2.8239248) q[0];
sx q[0];
rz(-0.57967107) q[0];
x q[1];
rz(-0.082809049) q[2];
sx q[2];
rz(-2.2958404) q[2];
sx q[2];
rz(-2.3963181) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8926516) q[1];
sx q[1];
rz(-0.59160691) q[1];
sx q[1];
rz(-2.4581562) q[1];
rz(-pi) q[2];
rz(1.799753) q[3];
sx q[3];
rz(-1.6456283) q[3];
sx q[3];
rz(-0.018939806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5372411) q[2];
sx q[2];
rz(-1.9862572) q[2];
sx q[2];
rz(-1.381116) q[2];
rz(0.85754496) q[3];
sx q[3];
rz(-2.2159135) q[3];
sx q[3];
rz(3.082357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9581167) q[0];
sx q[0];
rz(-0.97267946) q[0];
sx q[0];
rz(-0.78829515) q[0];
rz(-0.46701416) q[1];
sx q[1];
rz(-2.4772418) q[1];
sx q[1];
rz(-2.3036387) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70650202) q[0];
sx q[0];
rz(-1.52924) q[0];
sx q[0];
rz(-2.6986319) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5068717) q[2];
sx q[2];
rz(-0.41213671) q[2];
sx q[2];
rz(1.2391547) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.26688901) q[1];
sx q[1];
rz(-1.6409487) q[1];
sx q[1];
rz(-1.288049) q[1];
rz(-pi) q[2];
rz(-0.31354745) q[3];
sx q[3];
rz(-1.6308074) q[3];
sx q[3];
rz(1.0911694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.18980846) q[2];
sx q[2];
rz(-0.34056792) q[2];
sx q[2];
rz(-1.5365441) q[2];
rz(-0.32478452) q[3];
sx q[3];
rz(-1.6580509) q[3];
sx q[3];
rz(1.8713846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.558641) q[0];
sx q[0];
rz(-1.0277717) q[0];
sx q[0];
rz(3.0943178) q[0];
rz(1.524823) q[1];
sx q[1];
rz(-2.6998417) q[1];
sx q[1];
rz(-1.5025274) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29848443) q[0];
sx q[0];
rz(-1.2969979) q[0];
sx q[0];
rz(-0.47852935) q[0];
x q[1];
rz(2.6870704) q[2];
sx q[2];
rz(-1.3873867) q[2];
sx q[2];
rz(-1.3647788) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8432118) q[1];
sx q[1];
rz(-1.3271558) q[1];
sx q[1];
rz(1.3158821) q[1];
rz(-pi) q[2];
rz(-1.7697236) q[3];
sx q[3];
rz(-1.891948) q[3];
sx q[3];
rz(-2.7407848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9105685) q[2];
sx q[2];
rz(-0.39414057) q[2];
sx q[2];
rz(-2.2310889) q[2];
rz(-2.2855811) q[3];
sx q[3];
rz(-1.4166219) q[3];
sx q[3];
rz(3.0205309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1055792) q[0];
sx q[0];
rz(-1.3904904) q[0];
sx q[0];
rz(-3.0730096) q[0];
rz(-0.71906459) q[1];
sx q[1];
rz(-1.9870575) q[1];
sx q[1];
rz(0.4206492) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7186169) q[0];
sx q[0];
rz(-2.291922) q[0];
sx q[0];
rz(2.3628985) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6967032) q[2];
sx q[2];
rz(-0.25212461) q[2];
sx q[2];
rz(2.0201403) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6955086) q[1];
sx q[1];
rz(-0.38818103) q[1];
sx q[1];
rz(0.10315) q[1];
rz(-1.2141905) q[3];
sx q[3];
rz(-2.187629) q[3];
sx q[3];
rz(2.9717556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5258096) q[2];
sx q[2];
rz(-0.97302786) q[2];
sx q[2];
rz(-0.77965492) q[2];
rz(-3.0443794) q[3];
sx q[3];
rz(-1.8153056) q[3];
sx q[3];
rz(1.9450845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5060611) q[0];
sx q[0];
rz(-0.66672915) q[0];
sx q[0];
rz(2.7147103) q[0];
rz(1.4510669) q[1];
sx q[1];
rz(-2.0207113) q[1];
sx q[1];
rz(2.5371187) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7830435) q[0];
sx q[0];
rz(-0.32354) q[0];
sx q[0];
rz(-1.6829674) q[0];
rz(-1.8556183) q[2];
sx q[2];
rz(-1.3249791) q[2];
sx q[2];
rz(-2.071683) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5005481) q[1];
sx q[1];
rz(-1.637008) q[1];
sx q[1];
rz(1.7524629) q[1];
rz(-pi) q[2];
rz(1.3434308) q[3];
sx q[3];
rz(-2.1220088) q[3];
sx q[3];
rz(-0.77533305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84588593) q[2];
sx q[2];
rz(-0.87575951) q[2];
sx q[2];
rz(0.24629822) q[2];
rz(-1.0036428) q[3];
sx q[3];
rz(-1.6946038) q[3];
sx q[3];
rz(3.057737) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2698782) q[0];
sx q[0];
rz(-0.074967472) q[0];
sx q[0];
rz(-2.9366034) q[0];
rz(2.9221453) q[1];
sx q[1];
rz(-1.7013902) q[1];
sx q[1];
rz(0.27935371) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.646307) q[0];
sx q[0];
rz(-1.4581465) q[0];
sx q[0];
rz(-1.1239284) q[0];
rz(-pi) q[1];
rz(-1.749239) q[2];
sx q[2];
rz(-0.25545909) q[2];
sx q[2];
rz(-2.5551772) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47100386) q[1];
sx q[1];
rz(-1.5140972) q[1];
sx q[1];
rz(-0.55135552) q[1];
rz(-2.1198307) q[3];
sx q[3];
rz(-2.5982937) q[3];
sx q[3];
rz(-0.70453139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4651935) q[2];
sx q[2];
rz(-0.78130829) q[2];
sx q[2];
rz(2.264273) q[2];
rz(-1.3389795) q[3];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1124697) q[0];
sx q[0];
rz(-2.5482197) q[0];
sx q[0];
rz(-2.8068338) q[0];
rz(0.33084694) q[1];
sx q[1];
rz(-2.0956764) q[1];
sx q[1];
rz(-0.61029339) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7205104) q[0];
sx q[0];
rz(-1.3925584) q[0];
sx q[0];
rz(-1.8946289) q[0];
rz(-pi) q[1];
rz(2.2525983) q[2];
sx q[2];
rz(-0.51883139) q[2];
sx q[2];
rz(-2.8270222) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.522532) q[1];
sx q[1];
rz(-2.5964964) q[1];
sx q[1];
rz(-2.7482102) q[1];
x q[2];
rz(0.957358) q[3];
sx q[3];
rz(-1.6803553) q[3];
sx q[3];
rz(0.027519634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.200013) q[2];
sx q[2];
rz(-1.0369077) q[2];
sx q[2];
rz(-1.3882136) q[2];
rz(-2.9774104) q[3];
sx q[3];
rz(-2.1037585) q[3];
sx q[3];
rz(2.8431456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2309017) q[0];
sx q[0];
rz(-1.1158442) q[0];
sx q[0];
rz(-3.0052465) q[0];
rz(-0.048010437) q[1];
sx q[1];
rz(-1.0142356) q[1];
sx q[1];
rz(1.8528574) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2292176) q[0];
sx q[0];
rz(-0.37579094) q[0];
sx q[0];
rz(2.1972607) q[0];
x q[1];
rz(-0.23581712) q[2];
sx q[2];
rz(-1.1305446) q[2];
sx q[2];
rz(2.4503675) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14827189) q[1];
sx q[1];
rz(-2.848756) q[1];
sx q[1];
rz(-0.50334986) q[1];
rz(-1.1111497) q[3];
sx q[3];
rz(-0.75503317) q[3];
sx q[3];
rz(-0.5577969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.62782225) q[2];
sx q[2];
rz(-0.0063535293) q[2];
sx q[2];
rz(2.9610236) q[2];
rz(2.0020961) q[3];
sx q[3];
rz(-1.5151016) q[3];
sx q[3];
rz(-0.12038055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(1.1861495) q[0];
sx q[0];
rz(-0.97301617) q[0];
sx q[0];
rz(1.0261616) q[0];
rz(1.8852662) q[1];
sx q[1];
rz(-1.8812814) q[1];
sx q[1];
rz(3.0756782) q[1];
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
rz(-0.77113192) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8280667) q[2];
sx q[2];
rz(-1.7343177) q[2];
sx q[2];
rz(2.7933957) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5951426) q[1];
sx q[1];
rz(-1.3347581) q[1];
sx q[1];
rz(-0.55320338) q[1];
rz(-pi) q[2];
rz(-2.5539909) q[3];
sx q[3];
rz(-1.3580702) q[3];
sx q[3];
rz(1.0472681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56254783) q[2];
sx q[2];
rz(-1.0479835) q[2];
sx q[2];
rz(2.6143796) q[2];
rz(-2.9910917) q[3];
sx q[3];
rz(-0.42948693) q[3];
sx q[3];
rz(-0.35587564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
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
rz(-2.0122779) q[2];
sx q[2];
rz(-1.8649351) q[2];
sx q[2];
rz(-1.9492016) q[2];
rz(-2.1470924) q[3];
sx q[3];
rz(-0.71190833) q[3];
sx q[3];
rz(1.7332903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
