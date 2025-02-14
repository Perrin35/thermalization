OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0786809) q[0];
sx q[0];
rz(-0.39830387) q[0];
sx q[0];
rz(2.8226573) q[0];
rz(-0.4660663) q[1];
sx q[1];
rz(4.5748413) q[1];
sx q[1];
rz(9.6984697) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84209937) q[0];
sx q[0];
rz(-0.84171879) q[0];
sx q[0];
rz(-1.0925348) q[0];
rz(0.75184043) q[2];
sx q[2];
rz(-0.32773563) q[2];
sx q[2];
rz(0.40204266) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4727386) q[1];
sx q[1];
rz(-1.7484807) q[1];
sx q[1];
rz(-0.11858003) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1496278) q[3];
sx q[3];
rz(-0.85794373) q[3];
sx q[3];
rz(1.9496535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8493001) q[2];
sx q[2];
rz(-2.2165074) q[2];
sx q[2];
rz(-1.383847) q[2];
rz(-2.2960831) q[3];
sx q[3];
rz(-3.129831) q[3];
sx q[3];
rz(0.99172926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.200835) q[0];
sx q[0];
rz(-2.2469914) q[0];
sx q[0];
rz(1.0351329) q[0];
rz(1.9738522) q[1];
sx q[1];
rz(-0.20226856) q[1];
sx q[1];
rz(0.77441961) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7709633) q[0];
sx q[0];
rz(-1.9815195) q[0];
sx q[0];
rz(-2.109101) q[0];
rz(0.90986386) q[2];
sx q[2];
rz(-1.8468886) q[2];
sx q[2];
rz(-0.7731437) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.32774156) q[1];
sx q[1];
rz(-0.5264414) q[1];
sx q[1];
rz(1.3415643) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.033942698) q[3];
sx q[3];
rz(-1.5205624) q[3];
sx q[3];
rz(-1.4242581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3356129) q[2];
sx q[2];
rz(-0.077294417) q[2];
sx q[2];
rz(-2.198931) q[2];
rz(-3.0521657) q[3];
sx q[3];
rz(-0.663203) q[3];
sx q[3];
rz(-2.8360143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5494004) q[0];
sx q[0];
rz(-2.5553199) q[0];
sx q[0];
rz(2.9495682) q[0];
rz(-0.61344719) q[1];
sx q[1];
rz(-2.6934721) q[1];
sx q[1];
rz(0.014785756) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3131994) q[0];
sx q[0];
rz(-1.9678741) q[0];
sx q[0];
rz(-0.81531378) q[0];
rz(-pi) q[1];
rz(2.4528265) q[2];
sx q[2];
rz(-1.3564132) q[2];
sx q[2];
rz(-0.87300473) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.046349) q[1];
sx q[1];
rz(-0.63094891) q[1];
sx q[1];
rz(1.9722523) q[1];
rz(0.1802894) q[3];
sx q[3];
rz(-0.064924463) q[3];
sx q[3];
rz(-2.8466259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6141367) q[2];
sx q[2];
rz(-1.5469896) q[2];
sx q[2];
rz(-0.0097489348) q[2];
rz(2.5629432) q[3];
sx q[3];
rz(-2.0895683) q[3];
sx q[3];
rz(-0.78154045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0233199) q[0];
sx q[0];
rz(-2.3319722) q[0];
sx q[0];
rz(0.6672346) q[0];
rz(2.1369333) q[1];
sx q[1];
rz(-0.15548448) q[1];
sx q[1];
rz(1.7612693) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3074493) q[0];
sx q[0];
rz(-0.47795313) q[0];
sx q[0];
rz(-1.5909999) q[0];
rz(-pi) q[1];
rz(1.097408) q[2];
sx q[2];
rz(-1.2815099) q[2];
sx q[2];
rz(-2.8662967) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8246824) q[1];
sx q[1];
rz(-1.0774634) q[1];
sx q[1];
rz(1.6073336) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9033086) q[3];
sx q[3];
rz(-1.0589979) q[3];
sx q[3];
rz(2.3121367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.37463793) q[2];
sx q[2];
rz(-1.1344323) q[2];
sx q[2];
rz(-0.053442027) q[2];
rz(2.5080569) q[3];
sx q[3];
rz(-0.41826785) q[3];
sx q[3];
rz(-0.0088648349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2706914) q[0];
sx q[0];
rz(-2.5636261) q[0];
sx q[0];
rz(2.0722678) q[0];
rz(1.2879734) q[1];
sx q[1];
rz(-0.063871495) q[1];
sx q[1];
rz(0.091648253) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1250138) q[0];
sx q[0];
rz(-1.582394) q[0];
sx q[0];
rz(1.5931104) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32817082) q[2];
sx q[2];
rz(-2.0053021) q[2];
sx q[2];
rz(-0.19319867) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.55404592) q[1];
sx q[1];
rz(-2.5623706) q[1];
sx q[1];
rz(2.278797) q[1];
rz(-pi) q[2];
rz(1.5318457) q[3];
sx q[3];
rz(-0.80128968) q[3];
sx q[3];
rz(1.0391446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2962467) q[2];
sx q[2];
rz(-0.77719921) q[2];
sx q[2];
rz(-3.0261611) q[2];
rz(0.7655862) q[3];
sx q[3];
rz(-1.4976394) q[3];
sx q[3];
rz(-2.7287741) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63824832) q[0];
sx q[0];
rz(-0.64953506) q[0];
sx q[0];
rz(-2.5576538) q[0];
rz(3.0931547) q[1];
sx q[1];
rz(-2.9190639) q[1];
sx q[1];
rz(-2.4979874) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38356919) q[0];
sx q[0];
rz(-1.8789046) q[0];
sx q[0];
rz(-0.049806101) q[0];
x q[1];
rz(-2.6765049) q[2];
sx q[2];
rz(-2.1784867) q[2];
sx q[2];
rz(-1.6179784) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4847624) q[1];
sx q[1];
rz(-1.4387913) q[1];
sx q[1];
rz(-0.85797347) q[1];
rz(-2.3619283) q[3];
sx q[3];
rz(-1.7105967) q[3];
sx q[3];
rz(-1.9817022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.863997) q[2];
sx q[2];
rz(-2.3433351) q[2];
sx q[2];
rz(-2.4994728) q[2];
rz(-3.0025205) q[3];
sx q[3];
rz(-0.13834794) q[3];
sx q[3];
rz(-1.4596435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.36215255) q[0];
sx q[0];
rz(-2.7208949) q[0];
sx q[0];
rz(-2.5162589) q[0];
rz(-2.9026237) q[1];
sx q[1];
rz(-2.9049554) q[1];
sx q[1];
rz(-1.7854569) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5820393) q[0];
sx q[0];
rz(-1.6252717) q[0];
sx q[0];
rz(0.031886851) q[0];
rz(-1.2797592) q[2];
sx q[2];
rz(-1.2087529) q[2];
sx q[2];
rz(-2.9933528) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6125638) q[1];
sx q[1];
rz(-2.7269438) q[1];
sx q[1];
rz(2.0842444) q[1];
x q[2];
rz(0.21276833) q[3];
sx q[3];
rz(-1.4179061) q[3];
sx q[3];
rz(-2.2365856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8465392) q[2];
sx q[2];
rz(-1.876538) q[2];
sx q[2];
rz(-2.1758101) q[2];
rz(-2.8948808) q[3];
sx q[3];
rz(-1.2652946) q[3];
sx q[3];
rz(-1.770796) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5024289) q[0];
sx q[0];
rz(-2.736709) q[0];
sx q[0];
rz(-0.13614458) q[0];
rz(0.73774058) q[1];
sx q[1];
rz(-2.9164011) q[1];
sx q[1];
rz(-1.9053316) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8566187) q[0];
sx q[0];
rz(-0.98561937) q[0];
sx q[0];
rz(-2.7057458) q[0];
rz(-pi) q[1];
rz(-1.3278373) q[2];
sx q[2];
rz(-2.4907673) q[2];
sx q[2];
rz(-1.5200159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.511179) q[1];
sx q[1];
rz(-1.8466338) q[1];
sx q[1];
rz(-0.44075573) q[1];
rz(-pi) q[2];
rz(1.1499722) q[3];
sx q[3];
rz(-2.4462786) q[3];
sx q[3];
rz(-3.1041404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4514734) q[2];
sx q[2];
rz(-1.5855007) q[2];
sx q[2];
rz(2.1915009) q[2];
rz(-2.7443366) q[3];
sx q[3];
rz(-0.49644956) q[3];
sx q[3];
rz(0.44601405) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5600679) q[0];
sx q[0];
rz(-2.4604605) q[0];
sx q[0];
rz(-0.11823046) q[0];
rz(1.0826348) q[1];
sx q[1];
rz(-1.2018459) q[1];
sx q[1];
rz(-1.3523678) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3782458) q[0];
sx q[0];
rz(-1.1781426) q[0];
sx q[0];
rz(1.4598085) q[0];
rz(-pi) q[1];
rz(-0.544205) q[2];
sx q[2];
rz(-1.4883071) q[2];
sx q[2];
rz(-1.5832886) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0491421) q[1];
sx q[1];
rz(-1.0901206) q[1];
sx q[1];
rz(0.1436101) q[1];
rz(2.1264137) q[3];
sx q[3];
rz(-2.5930282) q[3];
sx q[3];
rz(0.11840222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.81789404) q[2];
sx q[2];
rz(-0.69043058) q[2];
sx q[2];
rz(-2.3696118) q[2];
rz(-3.0585994) q[3];
sx q[3];
rz(-1.8738184) q[3];
sx q[3];
rz(-2.6193589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25987396) q[0];
sx q[0];
rz(-0.75749713) q[0];
sx q[0];
rz(-2.3990193) q[0];
rz(-1.8450129) q[1];
sx q[1];
rz(-0.80755889) q[1];
sx q[1];
rz(-0.47320941) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6542235) q[0];
sx q[0];
rz(-2.8608436) q[0];
sx q[0];
rz(1.8868179) q[0];
x q[1];
rz(2.8933018) q[2];
sx q[2];
rz(-1.3220698) q[2];
sx q[2];
rz(-1.3922214) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.73903436) q[1];
sx q[1];
rz(-1.1087045) q[1];
sx q[1];
rz(-0.38886221) q[1];
x q[2];
rz(-1.3172174) q[3];
sx q[3];
rz(-0.84159708) q[3];
sx q[3];
rz(-2.1104476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5903198) q[2];
sx q[2];
rz(-1.4164305) q[2];
sx q[2];
rz(1.2335221) q[2];
rz(2.7963855) q[3];
sx q[3];
rz(-0.39185169) q[3];
sx q[3];
rz(-2.3046618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.681916) q[0];
sx q[0];
rz(-1.7727333) q[0];
sx q[0];
rz(1.75417) q[0];
rz(0.058902901) q[1];
sx q[1];
rz(-1.7210996) q[1];
sx q[1];
rz(0.65239418) q[1];
rz(0.38517135) q[2];
sx q[2];
rz(-1.9829911) q[2];
sx q[2];
rz(0.34172716) q[2];
rz(2.6860438) q[3];
sx q[3];
rz(-0.98916097) q[3];
sx q[3];
rz(-2.6670585) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
