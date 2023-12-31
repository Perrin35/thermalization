OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.964736) q[0];
sx q[0];
rz(-2.2778947) q[0];
sx q[0];
rz(2.9404844) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(4.9464524) q[1];
sx q[1];
rz(10.051605) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1683567) q[0];
sx q[0];
rz(-1.7587979) q[0];
sx q[0];
rz(0.12113916) q[0];
rz(-0.5958545) q[2];
sx q[2];
rz(-2.7239954) q[2];
sx q[2];
rz(-0.28996224) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.99170291) q[1];
sx q[1];
rz(-2.457379) q[1];
sx q[1];
rz(-2.0725155) q[1];
x q[2];
rz(2.9653373) q[3];
sx q[3];
rz(-1.8415383) q[3];
sx q[3];
rz(-2.3034629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(1.2935151) q[2];
rz(2.0375997) q[3];
sx q[3];
rz(-2.2938426) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936546) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(-2.1564116) q[0];
rz(-2.7424116) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(-0.10736297) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3283008) q[0];
sx q[0];
rz(-1.0257226) q[0];
sx q[0];
rz(-2.6585447) q[0];
rz(-pi) q[1];
rz(-0.58789247) q[2];
sx q[2];
rz(-1.3304454) q[2];
sx q[2];
rz(-0.27871486) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0641891) q[1];
sx q[1];
rz(-1.2014376) q[1];
sx q[1];
rz(0.87537745) q[1];
x q[2];
rz(2.7483838) q[3];
sx q[3];
rz(-2.7244096) q[3];
sx q[3];
rz(1.259491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5169516) q[2];
sx q[2];
rz(-1.3890283) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(-0.055756904) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(0.24584809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645638) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(-0.04034986) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-0.89825392) q[1];
sx q[1];
rz(-1.0823762) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6425991) q[0];
sx q[0];
rz(-1.9369619) q[0];
sx q[0];
rz(2.5901592) q[0];
x q[1];
rz(2.9897887) q[2];
sx q[2];
rz(-1.9215343) q[2];
sx q[2];
rz(-2.0098067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3861474) q[1];
sx q[1];
rz(-1.8090994) q[1];
sx q[1];
rz(0.2281245) q[1];
rz(-pi) q[2];
rz(-0.97614395) q[3];
sx q[3];
rz(-2.757132) q[3];
sx q[3];
rz(0.36851766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(-2.5877) q[2];
rz(1.6484377) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(0.55251399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4937113) q[0];
sx q[0];
rz(-2.557502) q[0];
sx q[0];
rz(0.035382263) q[0];
rz(2.1454051) q[1];
sx q[1];
rz(-1.6085947) q[1];
sx q[1];
rz(-2.6534973) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1195388) q[0];
sx q[0];
rz(-1.4931803) q[0];
sx q[0];
rz(-2.0366497) q[0];
rz(-2.343802) q[2];
sx q[2];
rz(-2.3893223) q[2];
sx q[2];
rz(-0.70300245) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.19255895) q[1];
sx q[1];
rz(-0.20201905) q[1];
sx q[1];
rz(2.2662524) q[1];
rz(2.048269) q[3];
sx q[3];
rz(-2.100088) q[3];
sx q[3];
rz(2.570591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7867243) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(0.50746894) q[2];
rz(1.9594225) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(-1.2549843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(0.59947157) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(2.7713293) q[0];
rz(-2.8778991) q[1];
sx q[1];
rz(-1.5636684) q[1];
sx q[1];
rz(-2.945074) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5614583) q[0];
sx q[0];
rz(-0.41813865) q[0];
sx q[0];
rz(2.3925376) q[0];
rz(-pi) q[1];
rz(-1.6351661) q[2];
sx q[2];
rz(-0.46191051) q[2];
sx q[2];
rz(2.7342352) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6477485) q[1];
sx q[1];
rz(-1.7462247) q[1];
sx q[1];
rz(-0.26248787) q[1];
x q[2];
rz(-2.5971562) q[3];
sx q[3];
rz(-2.8215373) q[3];
sx q[3];
rz(0.65359945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.013441) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(2.2244942) q[2];
rz(-2.3069978) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.82814) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(-0.20203461) q[0];
rz(1.0728041) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(1.3938168) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48702792) q[0];
sx q[0];
rz(-1.2396887) q[0];
sx q[0];
rz(-1.4902671) q[0];
rz(-pi) q[1];
rz(-1.1962842) q[2];
sx q[2];
rz(-1.4669344) q[2];
sx q[2];
rz(-1.9749157) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2653633) q[1];
sx q[1];
rz(-0.47788844) q[1];
sx q[1];
rz(-0.96546121) q[1];
rz(-pi) q[2];
rz(-0.73563852) q[3];
sx q[3];
rz(-1.2218352) q[3];
sx q[3];
rz(2.5603106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94971925) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(0.8423155) q[2];
rz(2.0541644) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(-1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52952805) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(-1.0213617) q[0];
rz(2.0102603) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(-2.9313415) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70049858) q[0];
sx q[0];
rz(-0.14963089) q[0];
sx q[0];
rz(-2.0307226) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5049997) q[2];
sx q[2];
rz(-1.734126) q[2];
sx q[2];
rz(2.3492299) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5308295) q[1];
sx q[1];
rz(-1.3794583) q[1];
sx q[1];
rz(0.32236871) q[1];
rz(-pi) q[2];
rz(-0.85150163) q[3];
sx q[3];
rz(-2.0264894) q[3];
sx q[3];
rz(1.15629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3147099) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(1.38114) q[2];
rz(-2.3794877) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(-0.41346082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522266) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(2.5543509) q[0];
rz(-0.44627055) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(2.1648724) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644972) q[0];
sx q[0];
rz(-1.7194887) q[0];
sx q[0];
rz(-2.68481) q[0];
rz(-pi) q[1];
rz(1.4090528) q[2];
sx q[2];
rz(-1.6502893) q[2];
sx q[2];
rz(-0.0083991945) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8969438) q[1];
sx q[1];
rz(-1.2243425) q[1];
sx q[1];
rz(-2.0379373) q[1];
rz(-pi) q[2];
rz(0.55240734) q[3];
sx q[3];
rz(-2.2033785) q[3];
sx q[3];
rz(1.6747024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4010767) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-0.55244279) q[2];
rz(-2.6756514) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(-2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733646) q[0];
sx q[0];
rz(-2.3478274) q[0];
sx q[0];
rz(0.92064944) q[0];
rz(-0.051041516) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(-0.89909536) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5024922) q[0];
sx q[0];
rz(-2.149625) q[0];
sx q[0];
rz(-2.6274908) q[0];
rz(-pi) q[1];
rz(2.0752418) q[2];
sx q[2];
rz(-2.5214508) q[2];
sx q[2];
rz(-1.3401741) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8583262) q[1];
sx q[1];
rz(-0.35052931) q[1];
sx q[1];
rz(1.2298898) q[1];
rz(-pi) q[2];
rz(-1.5713463) q[3];
sx q[3];
rz(-2.7714202) q[3];
sx q[3];
rz(2.7412424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94545323) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(-0.83834046) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8024682) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(-1.1178281) q[0];
rz(-0.71240187) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(0.39696473) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3544281) q[0];
sx q[0];
rz(-2.0485326) q[0];
sx q[0];
rz(-0.54425311) q[0];
rz(-pi) q[1];
rz(2.8357387) q[2];
sx q[2];
rz(-1.2199739) q[2];
sx q[2];
rz(-0.84301126) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46170235) q[1];
sx q[1];
rz(-1.2352518) q[1];
sx q[1];
rz(0.093613503) q[1];
rz(0.40542094) q[3];
sx q[3];
rz(-1.1076895) q[3];
sx q[3];
rz(0.4511569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1931856) q[2];
sx q[2];
rz(-0.90017515) q[2];
sx q[2];
rz(-0.096654264) q[2];
rz(-1.4139253) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(-1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5554572) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.2697521) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(0.66432129) q[2];
sx q[2];
rz(-0.13673377) q[2];
sx q[2];
rz(-2.7665334) q[2];
rz(1.9962911) q[3];
sx q[3];
rz(-2.8588061) q[3];
sx q[3];
rz(-2.2929946) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
