OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71896267) q[0];
sx q[0];
rz(2.842272) q[0];
sx q[0];
rz(8.9301803) q[0];
rz(-5.1410723) q[1];
sx q[1];
rz(2.1358868) q[1];
sx q[1];
rz(11.436643) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8324234) q[0];
sx q[0];
rz(-2.4465843) q[0];
sx q[0];
rz(-2.8929936) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84096618) q[2];
sx q[2];
rz(-0.6305002) q[2];
sx q[2];
rz(1.3497242) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3898824) q[1];
sx q[1];
rz(-2.2851508) q[1];
sx q[1];
rz(2.8044279) q[1];
rz(-1.0315597) q[3];
sx q[3];
rz(-2.4162216) q[3];
sx q[3];
rz(-1.6417208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.31892458) q[2];
sx q[2];
rz(-1.8225887) q[2];
sx q[2];
rz(-2.2887716) q[2];
rz(-1.3301814) q[3];
sx q[3];
rz(-2.4460402) q[3];
sx q[3];
rz(3.0947963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3332719) q[0];
sx q[0];
rz(-0.051017314) q[0];
sx q[0];
rz(-1.6774696) q[0];
rz(1.6630215) q[1];
sx q[1];
rz(-1.9824948) q[1];
sx q[1];
rz(-1.2082072) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0883597) q[0];
sx q[0];
rz(-2.3908983) q[0];
sx q[0];
rz(0.73221598) q[0];
rz(-pi) q[1];
rz(2.0457532) q[2];
sx q[2];
rz(-1.37687) q[2];
sx q[2];
rz(2.3467807) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8025) q[1];
sx q[1];
rz(-0.29401699) q[1];
sx q[1];
rz(-1.1400998) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3974819) q[3];
sx q[3];
rz(-1.1220166) q[3];
sx q[3];
rz(0.059343222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6055484) q[2];
sx q[2];
rz(-0.41993419) q[2];
sx q[2];
rz(1.6571244) q[2];
rz(-3.1260955) q[3];
sx q[3];
rz(-1.213538) q[3];
sx q[3];
rz(1.9626455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14979664) q[0];
sx q[0];
rz(-1.3301671) q[0];
sx q[0];
rz(-2.8699744) q[0];
rz(-0.89667165) q[1];
sx q[1];
rz(-2.6397557) q[1];
sx q[1];
rz(1.2976049) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.04682) q[0];
sx q[0];
rz(-1.5808788) q[0];
sx q[0];
rz(1.1077513) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6372725) q[2];
sx q[2];
rz(-1.9078622) q[2];
sx q[2];
rz(0.45318174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.43186346) q[1];
sx q[1];
rz(-1.6846906) q[1];
sx q[1];
rz(1.3224056) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1362757) q[3];
sx q[3];
rz(-0.76197366) q[3];
sx q[3];
rz(0.055082037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4012332) q[2];
sx q[2];
rz(-2.3637502) q[2];
sx q[2];
rz(-3.0878301) q[2];
rz(-1.2939804) q[3];
sx q[3];
rz(-0.79837489) q[3];
sx q[3];
rz(-2.9062041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(0.2006705) q[0];
sx q[0];
rz(-2.5735452) q[0];
sx q[0];
rz(1.6773552) q[0];
rz(-1.1306521) q[1];
sx q[1];
rz(-2.407275) q[1];
sx q[1];
rz(-0.11437036) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0443886) q[0];
sx q[0];
rz(-2.1425793) q[0];
sx q[0];
rz(-1.0685789) q[0];
x q[1];
rz(-2.4855544) q[2];
sx q[2];
rz(-1.1636574) q[2];
sx q[2];
rz(-2.7380057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3072171) q[1];
sx q[1];
rz(-0.67612069) q[1];
sx q[1];
rz(-1.6130877) q[1];
rz(-pi) q[2];
rz(0.032764445) q[3];
sx q[3];
rz(-0.19812852) q[3];
sx q[3];
rz(-1.9447382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4734681) q[2];
sx q[2];
rz(-1.5359842) q[2];
sx q[2];
rz(-3.065897) q[2];
rz(0.43478742) q[3];
sx q[3];
rz(-1.1554759) q[3];
sx q[3];
rz(-0.079631478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39744034) q[0];
sx q[0];
rz(-1.3120774) q[0];
sx q[0];
rz(0.42006668) q[0];
rz(-2.9282667) q[1];
sx q[1];
rz(-1.7330287) q[1];
sx q[1];
rz(1.8992281) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5516226) q[0];
sx q[0];
rz(-1.9931355) q[0];
sx q[0];
rz(0.7606272) q[0];
rz(-pi) q[1];
rz(2.7105646) q[2];
sx q[2];
rz(-1.3044453) q[2];
sx q[2];
rz(-2.9961627) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.865001) q[1];
sx q[1];
rz(-1.5861771) q[1];
sx q[1];
rz(-1.4702142) q[1];
rz(0.20988864) q[3];
sx q[3];
rz(-2.915463) q[3];
sx q[3];
rz(1.9551639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3199978) q[2];
sx q[2];
rz(-0.91511202) q[2];
sx q[2];
rz(-0.37459174) q[2];
rz(-1.0125259) q[3];
sx q[3];
rz(-2.7495224) q[3];
sx q[3];
rz(-2.4969126) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11151611) q[0];
sx q[0];
rz(-2.6426297) q[0];
sx q[0];
rz(2.7110355) q[0];
rz(-2.997609) q[1];
sx q[1];
rz(-1.0824243) q[1];
sx q[1];
rz(-2.5103501) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4743606) q[0];
sx q[0];
rz(-2.1788414) q[0];
sx q[0];
rz(1.505018) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14603931) q[2];
sx q[2];
rz(-1.9873053) q[2];
sx q[2];
rz(3.0768968) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5938607) q[1];
sx q[1];
rz(-1.0775078) q[1];
sx q[1];
rz(-0.6823205) q[1];
rz(-pi) q[2];
rz(2.8064734) q[3];
sx q[3];
rz(-2.0010173) q[3];
sx q[3];
rz(2.8443401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97079078) q[2];
sx q[2];
rz(-2.0826714) q[2];
sx q[2];
rz(-1.4405174) q[2];
rz(-1.3614281) q[3];
sx q[3];
rz(-2.0378588) q[3];
sx q[3];
rz(-2.6543999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68706566) q[0];
sx q[0];
rz(-0.28959689) q[0];
sx q[0];
rz(-0.92426306) q[0];
rz(2.848792) q[1];
sx q[1];
rz(-1.2288789) q[1];
sx q[1];
rz(0.80088314) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7445114) q[0];
sx q[0];
rz(-1.1985221) q[0];
sx q[0];
rz(-2.8651994) q[0];
x q[1];
rz(-2.0037193) q[2];
sx q[2];
rz(-0.70919631) q[2];
sx q[2];
rz(-1.8877754) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0139069) q[1];
sx q[1];
rz(-3.032445) q[1];
sx q[1];
rz(2.8692895) q[1];
x q[2];
rz(-1.5436567) q[3];
sx q[3];
rz(-1.5487897) q[3];
sx q[3];
rz(0.95809551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.31876365) q[2];
sx q[2];
rz(-1.2206565) q[2];
sx q[2];
rz(1.8611543) q[2];
rz(2.473623) q[3];
sx q[3];
rz(-2.6536055) q[3];
sx q[3];
rz(-0.41332301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8845344) q[0];
sx q[0];
rz(-1.9030544) q[0];
sx q[0];
rz(-1.9588233) q[0];
rz(-0.23712748) q[1];
sx q[1];
rz(-3.0393937) q[1];
sx q[1];
rz(3.0822486) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099456121) q[0];
sx q[0];
rz(-1.8435394) q[0];
sx q[0];
rz(-0.12423652) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58299033) q[2];
sx q[2];
rz(-0.81524476) q[2];
sx q[2];
rz(1.5423403) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3415449) q[1];
sx q[1];
rz(-0.28452415) q[1];
sx q[1];
rz(1.9519898) q[1];
x q[2];
rz(-3.0640934) q[3];
sx q[3];
rz(-1.8935673) q[3];
sx q[3];
rz(0.68068824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.55988971) q[2];
sx q[2];
rz(-0.54055944) q[2];
sx q[2];
rz(-1.3524559) q[2];
rz(1.3672359) q[3];
sx q[3];
rz(-1.4227899) q[3];
sx q[3];
rz(2.2860315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.2818114) q[0];
sx q[0];
rz(-2.3563522) q[0];
sx q[0];
rz(2.6819041) q[0];
rz(-0.028060878) q[1];
sx q[1];
rz(-1.1496081) q[1];
sx q[1];
rz(-1.185816) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81521124) q[0];
sx q[0];
rz(-1.5259169) q[0];
sx q[0];
rz(-2.8759967) q[0];
rz(3.0621959) q[2];
sx q[2];
rz(-1.4697452) q[2];
sx q[2];
rz(-3.0011645) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.016077135) q[1];
sx q[1];
rz(-1.0927769) q[1];
sx q[1];
rz(-2.1090871) q[1];
rz(-pi) q[2];
rz(-2.7544153) q[3];
sx q[3];
rz(-2.4868591) q[3];
sx q[3];
rz(0.93965215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6431553) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(-0.15677491) q[2];
rz(0.59569851) q[3];
sx q[3];
rz(-2.7955293) q[3];
sx q[3];
rz(-0.025207635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6237685) q[0];
sx q[0];
rz(-2.1110004) q[0];
sx q[0];
rz(-0.18381707) q[0];
rz(0.078016438) q[1];
sx q[1];
rz(-1.0573496) q[1];
sx q[1];
rz(2.877291) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8562132) q[0];
sx q[0];
rz(-2.1341679) q[0];
sx q[0];
rz(1.8895288) q[0];
rz(1.043522) q[2];
sx q[2];
rz(-2.1855178) q[2];
sx q[2];
rz(1.6312903) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0288913) q[1];
sx q[1];
rz(-1.3910146) q[1];
sx q[1];
rz(-1.7305018) q[1];
rz(-1.6666344) q[3];
sx q[3];
rz(-1.0776099) q[3];
sx q[3];
rz(0.020190369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8514303) q[2];
sx q[2];
rz(-0.94238472) q[2];
sx q[2];
rz(0.22872049) q[2];
rz(-1.7372355) q[3];
sx q[3];
rz(-2.2018933) q[3];
sx q[3];
rz(3.0586045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1514773) q[0];
sx q[0];
rz(-0.84828068) q[0];
sx q[0];
rz(-1.620851) q[0];
rz(-2.4304541) q[1];
sx q[1];
rz(-1.8322721) q[1];
sx q[1];
rz(-2.4087404) q[1];
rz(-2.0696832) q[2];
sx q[2];
rz(-1.6098534) q[2];
sx q[2];
rz(0.44865566) q[2];
rz(2.5855999) q[3];
sx q[3];
rz(-1.0992194) q[3];
sx q[3];
rz(1.9627375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
