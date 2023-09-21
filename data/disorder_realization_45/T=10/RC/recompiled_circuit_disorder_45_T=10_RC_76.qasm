OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.835445) q[0];
sx q[0];
rz(2.4581576) q[0];
sx q[0];
rz(12.0876) q[0];
rz(0.03102826) q[1];
sx q[1];
rz(5.1217084) q[1];
sx q[1];
rz(6.9245467) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35533479) q[0];
sx q[0];
rz(-1.9541652) q[0];
sx q[0];
rz(-1.3134365) q[0];
rz(-pi) q[1];
rz(1.1429206) q[2];
sx q[2];
rz(-1.9170205) q[2];
sx q[2];
rz(1.8510173) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76170834) q[1];
sx q[1];
rz(-2.0239081) q[1];
sx q[1];
rz(-2.3945827) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6061329) q[3];
sx q[3];
rz(-0.70264953) q[3];
sx q[3];
rz(-0.040576064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.550094) q[2];
sx q[2];
rz(-1.9248362) q[2];
sx q[2];
rz(-2.6634898) q[2];
rz(-1.6889307) q[3];
sx q[3];
rz(-2.0958459) q[3];
sx q[3];
rz(-12/(7*pi)) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1801382) q[0];
sx q[0];
rz(-0.77471662) q[0];
sx q[0];
rz(-2.1226728) q[0];
rz(1.6628751) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(0.63308024) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2791726) q[0];
sx q[0];
rz(-2.3775568) q[0];
sx q[0];
rz(-3.0252302) q[0];
rz(2.1585629) q[2];
sx q[2];
rz(-2.2289742) q[2];
sx q[2];
rz(-0.36511974) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9004904) q[1];
sx q[1];
rz(-1.9111269) q[1];
sx q[1];
rz(1.9963005) q[1];
x q[2];
rz(-2.8301864) q[3];
sx q[3];
rz(-2.2283471) q[3];
sx q[3];
rz(2.1962375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7818266) q[2];
sx q[2];
rz(-2.7414913) q[2];
sx q[2];
rz(3.0241372) q[2];
rz(-2.84058) q[3];
sx q[3];
rz(-1.8071226) q[3];
sx q[3];
rz(1.3628179) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85787073) q[0];
sx q[0];
rz(-2.4283333) q[0];
sx q[0];
rz(-0.088407956) q[0];
rz(-1.1075426) q[1];
sx q[1];
rz(-0.91845599) q[1];
sx q[1];
rz(-0.69082469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3120964) q[0];
sx q[0];
rz(-1.7236241) q[0];
sx q[0];
rz(-1.6468847) q[0];
rz(-pi) q[1];
rz(-0.76491852) q[2];
sx q[2];
rz(-1.7419445) q[2];
sx q[2];
rz(0.82676065) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55331992) q[1];
sx q[1];
rz(-1.5968482) q[1];
sx q[1];
rz(-1.7528898) q[1];
rz(-2.7615943) q[3];
sx q[3];
rz(-1.2424386) q[3];
sx q[3];
rz(0.85193714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2174125) q[2];
sx q[2];
rz(-1.8841382) q[2];
sx q[2];
rz(0.10647354) q[2];
rz(-1.7051833) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(-1.4829372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0728264) q[0];
sx q[0];
rz(-2.9077353) q[0];
sx q[0];
rz(-2.6191214) q[0];
rz(-0.32896313) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(-0.55363384) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75008167) q[0];
sx q[0];
rz(-2.2405853) q[0];
sx q[0];
rz(-2.1704587) q[0];
x q[1];
rz(0.76866863) q[2];
sx q[2];
rz(-0.93066356) q[2];
sx q[2];
rz(1.4537571) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5428634) q[1];
sx q[1];
rz(-2.3485564) q[1];
sx q[1];
rz(-1.7584156) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1133075) q[3];
sx q[3];
rz(-2.3445498) q[3];
sx q[3];
rz(-0.68238168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.557495) q[2];
sx q[2];
rz(-2.2258874) q[2];
sx q[2];
rz(-0.49475691) q[2];
rz(-0.90302145) q[3];
sx q[3];
rz(-1.5477864) q[3];
sx q[3];
rz(1.2899227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6506127) q[0];
sx q[0];
rz(-2.3968093) q[0];
sx q[0];
rz(-1.0850798) q[0];
rz(0.96013534) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(0.18403149) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7407496) q[0];
sx q[0];
rz(-1.4662192) q[0];
sx q[0];
rz(-2.8507289) q[0];
x q[1];
rz(0.36501546) q[2];
sx q[2];
rz(-1.8301617) q[2];
sx q[2];
rz(0.58155453) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4689456) q[1];
sx q[1];
rz(-0.99860672) q[1];
sx q[1];
rz(-2.1395626) q[1];
rz(0.72599756) q[3];
sx q[3];
rz(-1.7222002) q[3];
sx q[3];
rz(-2.9263673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2400143) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(-1.1408172) q[2];
rz(0.42282894) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(1.1269349) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7891156) q[0];
sx q[0];
rz(-2.0334091) q[0];
sx q[0];
rz(0.88678962) q[0];
rz(-0.24041644) q[1];
sx q[1];
rz(-1.0188894) q[1];
sx q[1];
rz(2.9930847) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69420544) q[0];
sx q[0];
rz(-2.2930657) q[0];
sx q[0];
rz(1.7859875) q[0];
x q[1];
rz(1.4609023) q[2];
sx q[2];
rz(-1.3434778) q[2];
sx q[2];
rz(-1.1306888) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0420694) q[1];
sx q[1];
rz(-1.6798786) q[1];
sx q[1];
rz(1.8314349) q[1];
rz(-2.4466483) q[3];
sx q[3];
rz(-2.5616025) q[3];
sx q[3];
rz(0.65141962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.285816) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(-0.4883858) q[2];
rz(0.48940247) q[3];
sx q[3];
rz(-2.2198052) q[3];
sx q[3];
rz(1.874812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4483036) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(-2.3235902) q[0];
rz(1.2524293) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(2.0163527) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97911191) q[0];
sx q[0];
rz(-1.9444124) q[0];
sx q[0];
rz(-0.11066779) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0411885) q[2];
sx q[2];
rz(-0.33249582) q[2];
sx q[2];
rz(1.5645129) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4459671) q[1];
sx q[1];
rz(-1.274316) q[1];
sx q[1];
rz(0.28709025) q[1];
x q[2];
rz(1.0057698) q[3];
sx q[3];
rz(-2.1566448) q[3];
sx q[3];
rz(2.1573531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0685048) q[2];
sx q[2];
rz(-0.98098522) q[2];
sx q[2];
rz(-0.64458624) q[2];
rz(-1.6623496) q[3];
sx q[3];
rz(-0.95696604) q[3];
sx q[3];
rz(-1.4054327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1709764) q[0];
sx q[0];
rz(-3.0727486) q[0];
sx q[0];
rz(-1.6059426) q[0];
rz(1.2212785) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(1.01064) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37234026) q[0];
sx q[0];
rz(-1.7307889) q[0];
sx q[0];
rz(-2.3183547) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0755195) q[2];
sx q[2];
rz(-1.8349378) q[2];
sx q[2];
rz(-1.2656982) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.752754) q[1];
sx q[1];
rz(-2.8201582) q[1];
sx q[1];
rz(-2.3366117) q[1];
x q[2];
rz(2.741022) q[3];
sx q[3];
rz(-1.6555602) q[3];
sx q[3];
rz(2.8469507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5236139) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(-0.69407216) q[2];
rz(-0.56898919) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(-0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0555608) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(-0.49474299) q[0];
rz(-2.5231979) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(-3.0659952) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8537124) q[0];
sx q[0];
rz(-1.259385) q[0];
sx q[0];
rz(-2.1561554) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7463023) q[2];
sx q[2];
rz(-1.6298461) q[2];
sx q[2];
rz(-2.4621071) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9675933) q[1];
sx q[1];
rz(-0.6951957) q[1];
sx q[1];
rz(-0.16144067) q[1];
rz(-pi) q[2];
rz(2.5209386) q[3];
sx q[3];
rz(-1.422158) q[3];
sx q[3];
rz(-1.312048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8081234) q[2];
sx q[2];
rz(-1.418891) q[2];
sx q[2];
rz(-1.3405651) q[2];
rz(0.30424413) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8463523) q[0];
sx q[0];
rz(-1.8363991) q[0];
sx q[0];
rz(2.9647968) q[0];
rz(-1.2416174) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(0.44100824) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84687418) q[0];
sx q[0];
rz(-2.8537769) q[0];
sx q[0];
rz(2.3527282) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15372865) q[2];
sx q[2];
rz(-1.4970461) q[2];
sx q[2];
rz(-1.7388625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16377201) q[1];
sx q[1];
rz(-1.2308916) q[1];
sx q[1];
rz(1.1141067) q[1];
rz(-pi) q[2];
rz(-0.40481683) q[3];
sx q[3];
rz(-0.66473085) q[3];
sx q[3];
rz(-0.80617314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56197721) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(-2.8397172) q[2];
rz(-2.2484696) q[3];
sx q[3];
rz(-1.2877269) q[3];
sx q[3];
rz(-0.20475234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.4176035) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(3.1148615) q[1];
sx q[1];
rz(-1.6865128) q[1];
sx q[1];
rz(-1.7105688) q[1];
rz(2.2970207) q[2];
sx q[2];
rz(-1.209134) q[2];
sx q[2];
rz(-0.67187885) q[2];
rz(2.4206352) q[3];
sx q[3];
rz(-0.69098916) q[3];
sx q[3];
rz(0.50222764) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];