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
rz(-0.29932061) q[0];
sx q[0];
rz(-2.646995) q[0];
rz(1.142113) q[1];
sx q[1];
rz(-1.0057058) q[1];
sx q[1];
rz(-2.0118654) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8324234) q[0];
sx q[0];
rz(-0.69500837) q[0];
sx q[0];
rz(-0.2485991) q[0];
rz(-1.0725934) q[2];
sx q[2];
rz(-1.1668201) q[2];
sx q[2];
rz(-0.4046658) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8985624) q[1];
sx q[1];
rz(-2.3645325) q[1];
sx q[1];
rz(-1.9352566) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7143728) q[3];
sx q[3];
rz(-0.96517206) q[3];
sx q[3];
rz(0.96715121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.31892458) q[2];
sx q[2];
rz(-1.8225887) q[2];
sx q[2];
rz(2.2887716) q[2];
rz(-1.8114113) q[3];
sx q[3];
rz(-2.4460402) q[3];
sx q[3];
rz(0.046796355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80832076) q[0];
sx q[0];
rz(-0.051017314) q[0];
sx q[0];
rz(-1.6774696) q[0];
rz(-1.4785712) q[1];
sx q[1];
rz(-1.9824948) q[1];
sx q[1];
rz(1.9333855) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0639275) q[0];
sx q[0];
rz(-2.0443235) q[0];
sx q[0];
rz(-2.5350476) q[0];
rz(-pi) q[1];
rz(-1.0958395) q[2];
sx q[2];
rz(-1.7647226) q[2];
sx q[2];
rz(-2.3467807) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.959001) q[1];
sx q[1];
rz(-1.6920857) q[1];
sx q[1];
rz(-1.3022997) q[1];
rz(-pi) q[2];
rz(1.3974819) q[3];
sx q[3];
rz(-1.1220166) q[3];
sx q[3];
rz(-3.0822494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6055484) q[2];
sx q[2];
rz(-2.7216585) q[2];
sx q[2];
rz(-1.6571244) q[2];
rz(-3.1260955) q[3];
sx q[3];
rz(-1.213538) q[3];
sx q[3];
rz(-1.1789471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(0.14979664) q[0];
sx q[0];
rz(-1.8114256) q[0];
sx q[0];
rz(-0.27161828) q[0];
rz(-0.89667165) q[1];
sx q[1];
rz(-2.6397557) q[1];
sx q[1];
rz(-1.8439878) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51894278) q[0];
sx q[0];
rz(-1.1077767) q[0];
sx q[0];
rz(-3.1303236) q[0];
x q[1];
rz(0.62745749) q[2];
sx q[2];
rz(-2.5431923) q[2];
sx q[2];
rz(0.57777571) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5601968) q[1];
sx q[1];
rz(-2.8688258) q[1];
sx q[1];
rz(2.0062937) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5657229) q[3];
sx q[3];
rz(-2.3327565) q[3];
sx q[3];
rz(0.047732959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4012332) q[2];
sx q[2];
rz(-2.3637502) q[2];
sx q[2];
rz(-0.05376251) q[2];
rz(1.8476123) q[3];
sx q[3];
rz(-0.79837489) q[3];
sx q[3];
rz(0.23538858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9409222) q[0];
sx q[0];
rz(-0.5680474) q[0];
sx q[0];
rz(-1.6773552) q[0];
rz(1.1306521) q[1];
sx q[1];
rz(-0.73431763) q[1];
sx q[1];
rz(3.0272223) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0972041) q[0];
sx q[0];
rz(-2.1425793) q[0];
sx q[0];
rz(1.0685789) q[0];
rz(-pi) q[1];
rz(-1.0724154) q[2];
sx q[2];
rz(-0.97626462) q[2];
sx q[2];
rz(-1.6785113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3072171) q[1];
sx q[1];
rz(-2.465472) q[1];
sx q[1];
rz(-1.5285049) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5773729) q[3];
sx q[3];
rz(-1.7688171) q[3];
sx q[3];
rz(-1.1634367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4734681) q[2];
sx q[2];
rz(-1.6056085) q[2];
sx q[2];
rz(-3.065897) q[2];
rz(2.7068052) q[3];
sx q[3];
rz(-1.1554759) q[3];
sx q[3];
rz(0.079631478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7441523) q[0];
sx q[0];
rz(-1.8295153) q[0];
sx q[0];
rz(-2.721526) q[0];
rz(2.9282667) q[1];
sx q[1];
rz(-1.408564) q[1];
sx q[1];
rz(-1.2423645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7540589) q[0];
sx q[0];
rz(-2.2926169) q[0];
sx q[0];
rz(0.57768627) q[0];
x q[1];
rz(1.2790643) q[2];
sx q[2];
rz(-1.9856678) q[2];
sx q[2];
rz(1.8366829) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.99862) q[1];
sx q[1];
rz(-0.10174739) q[1];
sx q[1];
rz(1.7228026) q[1];
rz(-0.20988864) q[3];
sx q[3];
rz(-2.915463) q[3];
sx q[3];
rz(1.1864288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8215948) q[2];
sx q[2];
rz(-2.2264806) q[2];
sx q[2];
rz(-0.37459174) q[2];
rz(1.0125259) q[3];
sx q[3];
rz(-2.7495224) q[3];
sx q[3];
rz(-0.64468002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0300765) q[0];
sx q[0];
rz(-0.498963) q[0];
sx q[0];
rz(-0.43055713) q[0];
rz(-2.997609) q[1];
sx q[1];
rz(-2.0591683) q[1];
sx q[1];
rz(-0.63124257) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66723204) q[0];
sx q[0];
rz(-2.1788414) q[0];
sx q[0];
rz(-1.505018) q[0];
rz(-pi) q[1];
rz(-1.9912791) q[2];
sx q[2];
rz(-1.4373206) q[2];
sx q[2];
rz(1.4466637) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4858125) q[1];
sx q[1];
rz(-0.98193278) q[1];
sx q[1];
rz(-0.9649802) q[1];
x q[2];
rz(2.1926375) q[3];
sx q[3];
rz(-0.53882155) q[3];
sx q[3];
rz(-2.1486189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.97079078) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(-1.7010752) q[2];
rz(1.3614281) q[3];
sx q[3];
rz(-1.1037339) q[3];
sx q[3];
rz(0.48719278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.454527) q[0];
sx q[0];
rz(-0.28959689) q[0];
sx q[0];
rz(0.92426306) q[0];
rz(0.29280064) q[1];
sx q[1];
rz(-1.2288789) q[1];
sx q[1];
rz(2.3407095) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3970812) q[0];
sx q[0];
rz(-1.1985221) q[0];
sx q[0];
rz(-2.8651994) q[0];
rz(-2.7960294) q[2];
sx q[2];
rz(-2.2032732) q[2];
sx q[2];
rz(1.3407624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1276857) q[1];
sx q[1];
rz(-0.10914762) q[1];
sx q[1];
rz(-0.27230316) q[1];
rz(-pi) q[2];
rz(1.5436567) q[3];
sx q[3];
rz(-1.592803) q[3];
sx q[3];
rz(0.95809551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31876365) q[2];
sx q[2];
rz(-1.9209361) q[2];
sx q[2];
rz(1.2804383) q[2];
rz(-0.66796962) q[3];
sx q[3];
rz(-0.48798713) q[3];
sx q[3];
rz(-2.7282696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25705826) q[0];
sx q[0];
rz(-1.9030544) q[0];
sx q[0];
rz(-1.9588233) q[0];
rz(-0.23712748) q[1];
sx q[1];
rz(-0.10219899) q[1];
sx q[1];
rz(0.059344083) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5049669) q[0];
sx q[0];
rz(-1.4511746) q[0];
sx q[0];
rz(-1.2960394) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5586023) q[2];
sx q[2];
rz(-2.3263479) q[2];
sx q[2];
rz(-1.5423403) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0036111) q[1];
sx q[1];
rz(-1.6754158) q[1];
sx q[1];
rz(-1.8358747) q[1];
rz(-1.7982676) q[3];
sx q[3];
rz(-0.331628) q[3];
sx q[3];
rz(-0.92078269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.55988971) q[2];
sx q[2];
rz(-0.54055944) q[2];
sx q[2];
rz(-1.7891368) q[2];
rz(1.7743568) q[3];
sx q[3];
rz(-1.4227899) q[3];
sx q[3];
rz(0.8555612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8597813) q[0];
sx q[0];
rz(-0.78524041) q[0];
sx q[0];
rz(2.6819041) q[0];
rz(3.1135318) q[1];
sx q[1];
rz(-1.1496081) q[1];
sx q[1];
rz(1.9557767) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59212771) q[0];
sx q[0];
rz(-0.26927265) q[0];
sx q[0];
rz(-2.9721391) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90699754) q[2];
sx q[2];
rz(-0.12842783) q[2];
sx q[2];
rz(-0.80824404) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2433127) q[1];
sx q[1];
rz(-2.4377258) q[1];
sx q[1];
rz(2.3614592) q[1];
rz(-1.8529296) q[3];
sx q[3];
rz(-2.1698639) q[3];
sx q[3];
rz(1.4144858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6431553) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(-2.9848177) q[2];
rz(-2.5458941) q[3];
sx q[3];
rz(-0.34606338) q[3];
sx q[3];
rz(0.025207635) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51782411) q[0];
sx q[0];
rz(-1.0305923) q[0];
sx q[0];
rz(-0.18381707) q[0];
rz(3.0635762) q[1];
sx q[1];
rz(-1.0573496) q[1];
sx q[1];
rz(-2.877291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28537946) q[0];
sx q[0];
rz(-1.0074248) q[0];
sx q[0];
rz(1.8895288) q[0];
x q[1];
rz(-1.043522) q[2];
sx q[2];
rz(-0.95607483) q[2];
sx q[2];
rz(-1.5103024) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1127013) q[1];
sx q[1];
rz(-1.750578) q[1];
sx q[1];
rz(1.7305018) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6666344) q[3];
sx q[3];
rz(-2.0639827) q[3];
sx q[3];
rz(-0.020190369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2901624) q[2];
sx q[2];
rz(-0.94238472) q[2];
sx q[2];
rz(0.22872049) q[2];
rz(1.7372355) q[3];
sx q[3];
rz(-2.2018933) q[3];
sx q[3];
rz(0.082988113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9901154) q[0];
sx q[0];
rz(-0.84828068) q[0];
sx q[0];
rz(-1.620851) q[0];
rz(-0.71113853) q[1];
sx q[1];
rz(-1.3093206) q[1];
sx q[1];
rz(0.73285229) q[1];
rz(0.044471519) q[2];
sx q[2];
rz(-2.0692678) q[2];
sx q[2];
rz(-1.1434126) q[2];
rz(0.76821297) q[3];
sx q[3];
rz(-2.429001) q[3];
sx q[3];
rz(2.902239) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
