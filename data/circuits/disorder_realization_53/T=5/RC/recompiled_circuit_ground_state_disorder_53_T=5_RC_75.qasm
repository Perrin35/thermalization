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
rz(0.18000552) q[0];
rz(-2.242318) q[1];
sx q[1];
rz(-2.2744496) q[1];
sx q[1];
rz(-1.4049621) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4150909) q[0];
sx q[0];
rz(-0.15299061) q[0];
sx q[0];
rz(0.57463519) q[0];
x q[1];
rz(-0.55721941) q[2];
sx q[2];
rz(-1.7168593) q[2];
sx q[2];
rz(0.58196011) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22706341) q[1];
sx q[1];
rz(-1.8957912) q[1];
sx q[1];
rz(-2.0507647) q[1];
rz(-pi) q[2];
rz(-0.57056114) q[3];
sx q[3];
rz(-2.5789912) q[3];
sx q[3];
rz(2.1747052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11294242) q[2];
sx q[2];
rz(-2.3990227) q[2];
sx q[2];
rz(0.58153018) q[2];
rz(-0.90315008) q[3];
sx q[3];
rz(-1.4862458) q[3];
sx q[3];
rz(0.36710292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
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
rz(-1.7342664) q[0];
sx q[0];
rz(2.0043688) q[0];
rz(1.3251023) q[1];
sx q[1];
rz(-0.66665998) q[1];
sx q[1];
rz(2.4773662) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8839845) q[0];
sx q[0];
rz(-0.31766787) q[0];
sx q[0];
rz(-0.57967107) q[0];
x q[1];
rz(2.2975446) q[2];
sx q[2];
rz(-1.6327452) q[2];
sx q[2];
rz(-0.88050851) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9162619) q[1];
sx q[1];
rz(-1.2109149) q[1];
sx q[1];
rz(2.6612984) q[1];
rz(-1.8898403) q[3];
sx q[3];
rz(-0.24067146) q[3];
sx q[3];
rz(1.2414207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.183476) q[0];
sx q[0];
rz(-2.1689132) q[0];
sx q[0];
rz(0.78829515) q[0];
rz(2.6745785) q[1];
sx q[1];
rz(-2.4772418) q[1];
sx q[1];
rz(0.83795396) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88400052) q[0];
sx q[0];
rz(-2.0133475) q[0];
sx q[0];
rz(1.5248067) q[0];
x q[1];
rz(1.6347209) q[2];
sx q[2];
rz(-2.7294559) q[2];
sx q[2];
rz(1.9024379) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0671671) q[1];
sx q[1];
rz(-2.850499) q[1];
sx q[1];
rz(1.3240678) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8280452) q[3];
sx q[3];
rz(-1.6308074) q[3];
sx q[3];
rz(-1.0911694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.18980846) q[2];
sx q[2];
rz(-0.34056792) q[2];
sx q[2];
rz(-1.6050485) q[2];
rz(-2.8168081) q[3];
sx q[3];
rz(-1.6580509) q[3];
sx q[3];
rz(1.270208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.558641) q[0];
sx q[0];
rz(-2.1138209) q[0];
sx q[0];
rz(3.0943178) q[0];
rz(-1.6167697) q[1];
sx q[1];
rz(-0.44175092) q[1];
sx q[1];
rz(-1.6390653) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29848443) q[0];
sx q[0];
rz(-1.2969979) q[0];
sx q[0];
rz(0.47852935) q[0];
x q[1];
rz(-1.3672013) q[2];
sx q[2];
rz(-2.0171391) q[2];
sx q[2];
rz(2.8466895) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8432118) q[1];
sx q[1];
rz(-1.3271558) q[1];
sx q[1];
rz(1.3158821) q[1];
x q[2];
rz(1.371869) q[3];
sx q[3];
rz(-1.2496447) q[3];
sx q[3];
rz(2.7407848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9105685) q[2];
sx q[2];
rz(-2.7474521) q[2];
sx q[2];
rz(-2.2310889) q[2];
rz(2.2855811) q[3];
sx q[3];
rz(-1.7249707) q[3];
sx q[3];
rz(-0.12106171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0360134) q[0];
sx q[0];
rz(-1.7511022) q[0];
sx q[0];
rz(-0.068583071) q[0];
rz(0.71906459) q[1];
sx q[1];
rz(-1.9870575) q[1];
sx q[1];
rz(2.7209435) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42297573) q[0];
sx q[0];
rz(-2.291922) q[0];
sx q[0];
rz(-0.77869418) q[0];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.44608403) q[1];
sx q[1];
rz(-0.38818103) q[1];
sx q[1];
rz(-3.0384427) q[1];
x q[2];
rz(-2.6841389) q[3];
sx q[3];
rz(-0.70069271) q[3];
sx q[3];
rz(0.40231857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61578304) q[2];
sx q[2];
rz(-0.97302786) q[2];
sx q[2];
rz(0.77965492) q[2];
rz(-3.0443794) q[3];
sx q[3];
rz(-1.3262871) q[3];
sx q[3];
rz(1.1965082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-2.0207113) q[1];
sx q[1];
rz(-2.5371187) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.105851) q[0];
sx q[0];
rz(-1.5352016) q[0];
sx q[0];
rz(1.2491519) q[0];
rz(-pi) q[1];
rz(2.2996713) q[2];
sx q[2];
rz(-2.7675602) q[2];
sx q[2];
rz(1.9472515) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.27547243) q[1];
sx q[1];
rz(-2.9483612) q[1];
sx q[1];
rz(1.2190427) q[1];
x q[2];
rz(-0.35143654) q[3];
sx q[3];
rz(-2.5498423) q[3];
sx q[3];
rz(-1.9502673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2957067) q[2];
sx q[2];
rz(-2.2658331) q[2];
sx q[2];
rz(-0.24629822) q[2];
rz(1.0036428) q[3];
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
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87171444) q[0];
sx q[0];
rz(-3.0666252) q[0];
sx q[0];
rz(-0.20498928) q[0];
rz(2.9221453) q[1];
sx q[1];
rz(-1.4402025) q[1];
sx q[1];
rz(-0.27935371) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.646307) q[0];
sx q[0];
rz(-1.6834462) q[0];
sx q[0];
rz(1.1239284) q[0];
x q[1];
rz(0.046322919) q[2];
sx q[2];
rz(-1.8221107) q[2];
sx q[2];
rz(-0.77071079) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1336891) q[1];
sx q[1];
rz(-2.5876296) q[1];
sx q[1];
rz(3.0336607) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30531136) q[3];
sx q[3];
rz(-2.0274912) q[3];
sx q[3];
rz(1.8165464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4651935) q[2];
sx q[2];
rz(-2.3602844) q[2];
sx q[2];
rz(-2.264273) q[2];
rz(-1.3389795) q[3];
sx q[3];
rz(-0.78290144) q[3];
sx q[3];
rz(1.3907998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.029123) q[0];
sx q[0];
rz(-2.5482197) q[0];
sx q[0];
rz(0.33475885) q[0];
rz(-0.33084694) q[1];
sx q[1];
rz(-2.0956764) q[1];
sx q[1];
rz(0.61029339) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90972483) q[0];
sx q[0];
rz(-1.2522765) q[0];
sx q[0];
rz(2.9538049) q[0];
x q[1];
rz(2.2525983) q[2];
sx q[2];
rz(-0.51883139) q[2];
sx q[2];
rz(-2.8270222) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.6107294) q[1];
sx q[1];
rz(-1.3707152) q[1];
sx q[1];
rz(2.6310573) q[1];
rz(-0.957358) q[3];
sx q[3];
rz(-1.4612373) q[3];
sx q[3];
rz(-3.114073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94157964) q[2];
sx q[2];
rz(-2.104685) q[2];
sx q[2];
rz(1.753379) q[2];
rz(0.16418223) q[3];
sx q[3];
rz(-2.1037585) q[3];
sx q[3];
rz(2.8431456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91069094) q[0];
sx q[0];
rz(-2.0257484) q[0];
sx q[0];
rz(-3.0052465) q[0];
rz(-3.0935822) q[1];
sx q[1];
rz(-2.127357) q[1];
sx q[1];
rz(1.8528574) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2511906) q[0];
sx q[0];
rz(-1.8726761) q[0];
sx q[0];
rz(0.22731486) q[0];
rz(-pi) q[1];
rz(-2.0219649) q[2];
sx q[2];
rz(-1.7837614) q[2];
sx q[2];
rz(-2.1599744) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9077085) q[1];
sx q[1];
rz(-1.7104935) q[1];
sx q[1];
rz(-2.8833792) q[1];
rz(-pi) q[2];
rz(1.1111497) q[3];
sx q[3];
rz(-2.3865595) q[3];
sx q[3];
rz(-0.5577969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5137704) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9554431) q[0];
sx q[0];
rz(-0.97301617) q[0];
sx q[0];
rz(1.0261616) q[0];
rz(-1.8852662) q[1];
sx q[1];
rz(-1.2603113) q[1];
sx q[1];
rz(3.0756782) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3381461) q[0];
sx q[0];
rz(-2.1915276) q[0];
sx q[0];
rz(2.3148763) q[0];
rz(-2.8280667) q[2];
sx q[2];
rz(-1.7343177) q[2];
sx q[2];
rz(-2.7933957) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5464501) q[1];
sx q[1];
rz(-1.8068346) q[1];
sx q[1];
rz(0.55320338) q[1];
x q[2];
rz(0.37152901) q[3];
sx q[3];
rz(-2.5209628) q[3];
sx q[3];
rz(0.21658235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56254783) q[2];
sx q[2];
rz(-2.0936091) q[2];
sx q[2];
rz(-0.5272131) q[2];
rz(2.9910917) q[3];
sx q[3];
rz(-0.42948693) q[3];
sx q[3];
rz(0.35587564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.95407964) q[2];
sx q[2];
rz(-0.5250684) q[2];
sx q[2];
rz(-2.9697408) q[2];
rz(2.1971142) q[3];
sx q[3];
rz(-1.9347659) q[3];
sx q[3];
rz(0.61975382) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
