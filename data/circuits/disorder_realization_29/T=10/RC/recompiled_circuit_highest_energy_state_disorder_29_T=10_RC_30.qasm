OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2588377) q[0];
sx q[0];
rz(-0.14482276) q[0];
sx q[0];
rz(-0.75135279) q[0];
rz(-1.6169647) q[1];
sx q[1];
rz(5.3137988) q[1];
sx q[1];
rz(9.245524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0517901) q[0];
sx q[0];
rz(-1.2915011) q[0];
sx q[0];
rz(-2.0816878) q[0];
rz(-pi) q[1];
rz(3.1044311) q[2];
sx q[2];
rz(-2.4183309) q[2];
sx q[2];
rz(-0.0497555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7764531) q[1];
sx q[1];
rz(-1.6426597) q[1];
sx q[1];
rz(0.30993575) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.059285941) q[3];
sx q[3];
rz(-2.2036457) q[3];
sx q[3];
rz(2.2135753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.90193191) q[2];
sx q[2];
rz(-1.8472981) q[2];
sx q[2];
rz(-2.530976) q[2];
rz(2.2527952) q[3];
sx q[3];
rz(-2.461268) q[3];
sx q[3];
rz(0.73386598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4438542) q[0];
sx q[0];
rz(-0.30174169) q[0];
sx q[0];
rz(1.3296211) q[0];
rz(1.4854206) q[1];
sx q[1];
rz(-1.679436) q[1];
sx q[1];
rz(2.2775473) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5651503) q[0];
sx q[0];
rz(-2.8684542) q[0];
sx q[0];
rz(1.1821724) q[0];
rz(-pi) q[1];
rz(2.8301622) q[2];
sx q[2];
rz(-2.0144147) q[2];
sx q[2];
rz(1.6964427) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73693528) q[1];
sx q[1];
rz(-2.9184249) q[1];
sx q[1];
rz(-0.88921247) q[1];
rz(-pi) q[2];
rz(2.3584189) q[3];
sx q[3];
rz(-2.0370954) q[3];
sx q[3];
rz(-0.1649905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0835421) q[2];
sx q[2];
rz(-2.4594049) q[2];
sx q[2];
rz(2.9502499) q[2];
rz(-0.30609104) q[3];
sx q[3];
rz(-1.4602665) q[3];
sx q[3];
rz(2.2050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3092344) q[0];
sx q[0];
rz(-1.2795871) q[0];
sx q[0];
rz(0.424463) q[0];
rz(1.3407432) q[1];
sx q[1];
rz(-0.91573358) q[1];
sx q[1];
rz(-0.65139604) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1196647) q[0];
sx q[0];
rz(-2.9777973) q[0];
sx q[0];
rz(-1.5271565) q[0];
rz(-pi) q[1];
rz(1.3878787) q[2];
sx q[2];
rz(-1.9486893) q[2];
sx q[2];
rz(-1.479508) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60247682) q[1];
sx q[1];
rz(-2.0664363) q[1];
sx q[1];
rz(1.7784761) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3440787) q[3];
sx q[3];
rz(-1.7245088) q[3];
sx q[3];
rz(-2.8938229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1059619) q[2];
sx q[2];
rz(-2.0023846) q[2];
sx q[2];
rz(0.6967217) q[2];
rz(2.4524955) q[3];
sx q[3];
rz(-1.1545811) q[3];
sx q[3];
rz(-2.885163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43831393) q[0];
sx q[0];
rz(-2.8499481) q[0];
sx q[0];
rz(-2.1412204) q[0];
rz(-1.6814303) q[1];
sx q[1];
rz(-1.5337475) q[1];
sx q[1];
rz(-1.9283074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54736894) q[0];
sx q[0];
rz(-1.692019) q[0];
sx q[0];
rz(1.8531043) q[0];
rz(-pi) q[1];
rz(2.4126265) q[2];
sx q[2];
rz(-1.2104958) q[2];
sx q[2];
rz(-3.136022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91649969) q[1];
sx q[1];
rz(-1.2429534) q[1];
sx q[1];
rz(1.6986548) q[1];
x q[2];
rz(1.4198076) q[3];
sx q[3];
rz(-0.54364341) q[3];
sx q[3];
rz(1.0796987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0393684) q[2];
sx q[2];
rz(-0.82304707) q[2];
sx q[2];
rz(-0.064112045) q[2];
rz(2.4168849) q[3];
sx q[3];
rz(-1.9331845) q[3];
sx q[3];
rz(-0.39976111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940014) q[0];
sx q[0];
rz(-1.8444909) q[0];
sx q[0];
rz(-0.79175788) q[0];
rz(1.1119615) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(-1.8066822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2091141) q[0];
sx q[0];
rz(-0.9830342) q[0];
sx q[0];
rz(-0.57508075) q[0];
rz(-2.56078) q[2];
sx q[2];
rz(-1.6585095) q[2];
sx q[2];
rz(0.49077362) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68685442) q[1];
sx q[1];
rz(-1.6461186) q[1];
sx q[1];
rz(-1.3999062) q[1];
x q[2];
rz(2.9369257) q[3];
sx q[3];
rz(-2.137326) q[3];
sx q[3];
rz(0.03290225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2935334) q[2];
sx q[2];
rz(-2.2424825) q[2];
sx q[2];
rz(2.000957) q[2];
rz(3.0349777) q[3];
sx q[3];
rz(-1.5299608) q[3];
sx q[3];
rz(-0.30985668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3265729) q[0];
sx q[0];
rz(-2.6928379) q[0];
sx q[0];
rz(0.003224592) q[0];
rz(-3.0942753) q[1];
sx q[1];
rz(-1.686217) q[1];
sx q[1];
rz(1.7914194) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4989717) q[0];
sx q[0];
rz(-1.823171) q[0];
sx q[0];
rz(-3.0809771) q[0];
rz(-pi) q[1];
rz(-2.2590738) q[2];
sx q[2];
rz(-0.79791683) q[2];
sx q[2];
rz(1.5882815) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9474831) q[1];
sx q[1];
rz(-2.1305741) q[1];
sx q[1];
rz(0.67621381) q[1];
rz(-1.7058701) q[3];
sx q[3];
rz(-0.22264847) q[3];
sx q[3];
rz(2.6378353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.99001592) q[2];
sx q[2];
rz(-0.185597) q[2];
sx q[2];
rz(3.0604176) q[2];
rz(-0.89933991) q[3];
sx q[3];
rz(-0.81493655) q[3];
sx q[3];
rz(-1.2871294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4261632) q[0];
sx q[0];
rz(-1.3112712) q[0];
sx q[0];
rz(0.50450605) q[0];
rz(0.99507487) q[1];
sx q[1];
rz(-1.0202531) q[1];
sx q[1];
rz(-0.02034932) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5581455) q[0];
sx q[0];
rz(-0.37675315) q[0];
sx q[0];
rz(1.3135629) q[0];
rz(-pi) q[1];
rz(-1.3624914) q[2];
sx q[2];
rz(-1.4464738) q[2];
sx q[2];
rz(-0.97036874) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6183747) q[1];
sx q[1];
rz(-1.3058387) q[1];
sx q[1];
rz(2.9441903) q[1];
rz(-pi) q[2];
rz(1.6571088) q[3];
sx q[3];
rz(-1.619297) q[3];
sx q[3];
rz(-2.8750471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4703579) q[2];
sx q[2];
rz(-1.8378259) q[2];
sx q[2];
rz(-1.7669558) q[2];
rz(1.5291519) q[3];
sx q[3];
rz(-1.0498472) q[3];
sx q[3];
rz(0.092718743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2328211) q[0];
sx q[0];
rz(-3.0290373) q[0];
sx q[0];
rz(1.0821279) q[0];
rz(2.1807561) q[1];
sx q[1];
rz(-1.6583534) q[1];
sx q[1];
rz(-0.94295162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0884224) q[0];
sx q[0];
rz(-0.57500792) q[0];
sx q[0];
rz(2.1413598) q[0];
rz(-pi) q[1];
rz(-0.93859886) q[2];
sx q[2];
rz(-1.3812806) q[2];
sx q[2];
rz(2.0507484) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.3245736) q[1];
sx q[1];
rz(-2.1789411) q[1];
sx q[1];
rz(-2.2961246) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2632964) q[3];
sx q[3];
rz(-1.4715183) q[3];
sx q[3];
rz(0.57534522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1711787) q[2];
sx q[2];
rz(-1.3033988) q[2];
sx q[2];
rz(1.1478434) q[2];
rz(-2.7796699) q[3];
sx q[3];
rz(-1.7636834) q[3];
sx q[3];
rz(-1.3981147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73087937) q[0];
sx q[0];
rz(-2.78237) q[0];
sx q[0];
rz(0.10243375) q[0];
rz(2.5275285) q[1];
sx q[1];
rz(-1.026261) q[1];
sx q[1];
rz(0.817743) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2376643) q[0];
sx q[0];
rz(-0.92887628) q[0];
sx q[0];
rz(-0.68591811) q[0];
rz(-pi) q[1];
rz(3.0258437) q[2];
sx q[2];
rz(-2.2384746) q[2];
sx q[2];
rz(0.69224651) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9407258) q[1];
sx q[1];
rz(-1.5518909) q[1];
sx q[1];
rz(-2.3000642) q[1];
x q[2];
rz(2.6800167) q[3];
sx q[3];
rz(-0.99050039) q[3];
sx q[3];
rz(3.1212774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7927336) q[2];
sx q[2];
rz(-1.0772971) q[2];
sx q[2];
rz(0.31361541) q[2];
rz(-2.6775728) q[3];
sx q[3];
rz(-0.84657621) q[3];
sx q[3];
rz(2.2217506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2713276) q[0];
sx q[0];
rz(-2.3509404) q[0];
sx q[0];
rz(0.23649293) q[0];
rz(2.927921) q[1];
sx q[1];
rz(-2.2387319) q[1];
sx q[1];
rz(2.9170091) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29447039) q[0];
sx q[0];
rz(-1.0884388) q[0];
sx q[0];
rz(0.76089528) q[0];
rz(-pi) q[1];
rz(1.2055254) q[2];
sx q[2];
rz(-2.1238616) q[2];
sx q[2];
rz(1.6286563) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27222363) q[1];
sx q[1];
rz(-2.096513) q[1];
sx q[1];
rz(-1.438526) q[1];
rz(-pi) q[2];
rz(-2.2891232) q[3];
sx q[3];
rz(-1.895088) q[3];
sx q[3];
rz(-1.7773903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1647722) q[2];
sx q[2];
rz(-0.87146622) q[2];
sx q[2];
rz(0.15178794) q[2];
rz(3.1101036) q[3];
sx q[3];
rz(-1.7446691) q[3];
sx q[3];
rz(-0.12846863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.0781773) q[0];
sx q[0];
rz(-1.6178394) q[0];
sx q[0];
rz(2.7375426) q[0];
rz(1.0833441) q[1];
sx q[1];
rz(-1.5768408) q[1];
sx q[1];
rz(1.5595938) q[1];
rz(0.2486817) q[2];
sx q[2];
rz(-2.9338825) q[2];
sx q[2];
rz(2.6354229) q[2];
rz(2.8121171) q[3];
sx q[3];
rz(-1.4737045) q[3];
sx q[3];
rz(1.2274689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
