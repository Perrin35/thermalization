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
rz(0.49459767) q[0];
rz(-5.1410723) q[1];
sx q[1];
rz(2.1358868) q[1];
sx q[1];
rz(11.436643) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.513213) q[0];
sx q[0];
rz(-2.2404379) q[0];
sx q[0];
rz(-1.3684526) q[0];
x q[1];
rz(-0.84096618) q[2];
sx q[2];
rz(-0.6305002) q[2];
sx q[2];
rz(-1.3497242) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5482481) q[1];
sx q[1];
rz(-1.3181837) q[1];
sx q[1];
rz(2.3139364) q[1];
x q[2];
rz(-1.0315597) q[3];
sx q[3];
rz(-2.4162216) q[3];
sx q[3];
rz(-1.6417208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8226681) q[2];
sx q[2];
rz(-1.319004) q[2];
sx q[2];
rz(-0.85282105) q[2];
rz(1.8114113) q[3];
sx q[3];
rz(-0.69555247) q[3];
sx q[3];
rz(-3.0947963) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3332719) q[0];
sx q[0];
rz(-0.051017314) q[0];
sx q[0];
rz(1.464123) q[0];
rz(1.4785712) q[1];
sx q[1];
rz(-1.9824948) q[1];
sx q[1];
rz(-1.9333855) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0639275) q[0];
sx q[0];
rz(-2.0443235) q[0];
sx q[0];
rz(-0.60654503) q[0];
rz(-1.1651498) q[2];
sx q[2];
rz(-2.6313836) q[2];
sx q[2];
rz(2.0070397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.959001) q[1];
sx q[1];
rz(-1.6920857) q[1];
sx q[1];
rz(1.3022997) q[1];
rz(1.7441107) q[3];
sx q[3];
rz(-2.0195761) q[3];
sx q[3];
rz(0.059343222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.53604424) q[2];
sx q[2];
rz(-0.41993419) q[2];
sx q[2];
rz(1.6571244) q[2];
rz(-3.1260955) q[3];
sx q[3];
rz(-1.213538) q[3];
sx q[3];
rz(-1.1789471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14979664) q[0];
sx q[0];
rz(-1.8114256) q[0];
sx q[0];
rz(-2.8699744) q[0];
rz(0.89667165) q[1];
sx q[1];
rz(-2.6397557) q[1];
sx q[1];
rz(1.8439878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54416849) q[0];
sx q[0];
rz(-2.6784458) q[0];
sx q[0];
rz(-1.5482272) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9515368) q[2];
sx q[2];
rz(-2.0443161) q[2];
sx q[2];
rz(-1.2981594) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1101164) q[1];
sx q[1];
rz(-1.3240485) q[1];
sx q[1];
rz(-3.0241248) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3796261) q[3];
sx q[3];
rz(-1.5671258) q[3];
sx q[3];
rz(1.519561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4012332) q[2];
sx q[2];
rz(-2.3637502) q[2];
sx q[2];
rz(3.0878301) q[2];
rz(1.2939804) q[3];
sx q[3];
rz(-2.3432178) q[3];
sx q[3];
rz(0.23538858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9409222) q[0];
sx q[0];
rz(-2.5735452) q[0];
sx q[0];
rz(1.4642375) q[0];
rz(-2.0109406) q[1];
sx q[1];
rz(-2.407275) q[1];
sx q[1];
rz(-3.0272223) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30442023) q[0];
sx q[0];
rz(-2.3995598) q[0];
sx q[0];
rz(2.4993308) q[0];
rz(-pi) q[1];
rz(2.4855544) q[2];
sx q[2];
rz(-1.9779352) q[2];
sx q[2];
rz(-2.7380057) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83437551) q[1];
sx q[1];
rz(-0.67612069) q[1];
sx q[1];
rz(-1.6130877) q[1];
rz(-1.5642197) q[3];
sx q[3];
rz(-1.3727756) q[3];
sx q[3];
rz(-1.1634367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4734681) q[2];
sx q[2];
rz(-1.5359842) q[2];
sx q[2];
rz(3.065897) q[2];
rz(-0.43478742) q[3];
sx q[3];
rz(-1.9861168) q[3];
sx q[3];
rz(-0.079631478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39744034) q[0];
sx q[0];
rz(-1.3120774) q[0];
sx q[0];
rz(0.42006668) q[0];
rz(2.9282667) q[1];
sx q[1];
rz(-1.7330287) q[1];
sx q[1];
rz(1.2423645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5899701) q[0];
sx q[0];
rz(-1.9931355) q[0];
sx q[0];
rz(-0.7606272) q[0];
rz(0.57849291) q[2];
sx q[2];
rz(-2.6393386) q[2];
sx q[2];
rz(1.1963716) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2957569) q[1];
sx q[1];
rz(-1.4702262) q[1];
sx q[1];
rz(-0.015458903) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5228988) q[3];
sx q[3];
rz(-1.7918799) q[3];
sx q[3];
rz(-0.97122279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3199978) q[2];
sx q[2];
rz(-0.91511202) q[2];
sx q[2];
rz(-2.7670009) q[2];
rz(-1.0125259) q[3];
sx q[3];
rz(-0.39207021) q[3];
sx q[3];
rz(2.4969126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(-2.0591683) q[1];
sx q[1];
rz(2.5103501) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4743606) q[0];
sx q[0];
rz(-0.96275126) q[0];
sx q[0];
rz(-1.6365746) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1503136) q[2];
sx q[2];
rz(-1.704272) q[2];
sx q[2];
rz(-1.4466637) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5938607) q[1];
sx q[1];
rz(-1.0775078) q[1];
sx q[1];
rz(2.4592722) q[1];
x q[2];
rz(-0.9489551) q[3];
sx q[3];
rz(-0.53882155) q[3];
sx q[3];
rz(-2.1486189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1708019) q[2];
sx q[2];
rz(-2.0826714) q[2];
sx q[2];
rz(1.4405174) q[2];
rz(1.3614281) q[3];
sx q[3];
rz(-2.0378588) q[3];
sx q[3];
rz(-0.48719278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68706566) q[0];
sx q[0];
rz(-2.8519958) q[0];
sx q[0];
rz(-2.2173296) q[0];
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
rz(1.0821874) q[0];
sx q[0];
rz(-0.45978433) q[0];
sx q[0];
rz(-2.1807488) q[0];
rz(-pi) q[1];
rz(-2.0037193) q[2];
sx q[2];
rz(-2.4323963) q[2];
sx q[2];
rz(-1.2538172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.427721) q[1];
sx q[1];
rz(-1.5414951) q[1];
sx q[1];
rz(-0.10515611) q[1];
rz(-pi) q[2];
x q[2];
rz(1.597936) q[3];
sx q[3];
rz(-1.592803) q[3];
sx q[3];
rz(2.1834971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.31876365) q[2];
sx q[2];
rz(-1.2206565) q[2];
sx q[2];
rz(1.8611543) q[2];
rz(0.66796962) q[3];
sx q[3];
rz(-0.48798713) q[3];
sx q[3];
rz(-0.41332301) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25705826) q[0];
sx q[0];
rz(-1.9030544) q[0];
sx q[0];
rz(1.1827693) q[0];
rz(2.9044652) q[1];
sx q[1];
rz(-3.0393937) q[1];
sx q[1];
rz(-0.059344083) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0421365) q[0];
sx q[0];
rz(-1.8435394) q[0];
sx q[0];
rz(-0.12423652) q[0];
rz(2.4164532) q[2];
sx q[2];
rz(-1.1584917) q[2];
sx q[2];
rz(-2.7453842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3415449) q[1];
sx q[1];
rz(-0.28452415) q[1];
sx q[1];
rz(-1.9519898) q[1];
rz(1.8944727) q[3];
sx q[3];
rz(-1.4973065) q[3];
sx q[3];
rz(2.2761114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5817029) q[2];
sx q[2];
rz(-2.6010332) q[2];
sx q[2];
rz(-1.3524559) q[2];
rz(-1.7743568) q[3];
sx q[3];
rz(-1.7188027) q[3];
sx q[3];
rz(-2.2860315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8597813) q[0];
sx q[0];
rz(-0.78524041) q[0];
sx q[0];
rz(-2.6819041) q[0];
rz(-0.028060878) q[1];
sx q[1];
rz(-1.9919845) q[1];
sx q[1];
rz(1.185816) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59212771) q[0];
sx q[0];
rz(-2.87232) q[0];
sx q[0];
rz(-0.16945355) q[0];
rz(-pi) q[1];
rz(-1.469428) q[2];
sx q[2];
rz(-1.4918054) q[2];
sx q[2];
rz(1.7192507) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.016077135) q[1];
sx q[1];
rz(-1.0927769) q[1];
sx q[1];
rz(-1.0325055) q[1];
rz(-1.8529296) q[3];
sx q[3];
rz(-2.1698639) q[3];
sx q[3];
rz(-1.7271068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.49843732) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(-0.15677491) q[2];
rz(0.59569851) q[3];
sx q[3];
rz(-0.34606338) q[3];
sx q[3];
rz(0.025207635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6237685) q[0];
sx q[0];
rz(-2.1110004) q[0];
sx q[0];
rz(-0.18381707) q[0];
rz(3.0635762) q[1];
sx q[1];
rz(-1.0573496) q[1];
sx q[1];
rz(0.26430166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3027356) q[0];
sx q[0];
rz(-2.502901) q[0];
sx q[0];
rz(-2.6810718) q[0];
rz(2.0980706) q[2];
sx q[2];
rz(-2.1855178) q[2];
sx q[2];
rz(1.5103024) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8459611) q[1];
sx q[1];
rz(-0.23989883) q[1];
sx q[1];
rz(0.71883454) q[1];
rz(-pi) q[2];
rz(0.49510689) q[3];
sx q[3];
rz(-1.6551842) q[3];
sx q[3];
rz(1.5051248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2901624) q[2];
sx q[2];
rz(-0.94238472) q[2];
sx q[2];
rz(-2.9128722) q[2];
rz(1.7372355) q[3];
sx q[3];
rz(-0.93969932) q[3];
sx q[3];
rz(-0.082988113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1514773) q[0];
sx q[0];
rz(-0.84828068) q[0];
sx q[0];
rz(-1.620851) q[0];
rz(-0.71113853) q[1];
sx q[1];
rz(-1.3093206) q[1];
sx q[1];
rz(0.73285229) q[1];
rz(-1.489303) q[2];
sx q[2];
rz(-0.50028481) q[2];
sx q[2];
rz(-1.0505983) q[2];
rz(0.55599273) q[3];
sx q[3];
rz(-2.0423732) q[3];
sx q[3];
rz(-1.1788551) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
