OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5390227) q[0];
sx q[0];
rz(-2.5780926) q[0];
sx q[0];
rz(-0.45698693) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(-1.8656123) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7514483) q[0];
sx q[0];
rz(-1.4353308) q[0];
sx q[0];
rz(2.9209903) q[0];
rz(-pi) q[1];
rz(-3.0511191) q[2];
sx q[2];
rz(-1.6455368) q[2];
sx q[2];
rz(-0.59404101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4133271) q[1];
sx q[1];
rz(-1.3682377) q[1];
sx q[1];
rz(-3.0800746) q[1];
rz(-0.39235093) q[3];
sx q[3];
rz(-2.1808743) q[3];
sx q[3];
rz(-0.073079212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0322545) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(2.9764552) q[2];
rz(-1.5103546) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(-0.74222773) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6363268) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(-0.045036137) q[0];
rz(-2.8024407) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(0.80274686) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4883603) q[0];
sx q[0];
rz(-2.687504) q[0];
sx q[0];
rz(-1.0351719) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3834125) q[2];
sx q[2];
rz(-0.51088453) q[2];
sx q[2];
rz(-2.5201706) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3992577) q[1];
sx q[1];
rz(-2.5207673) q[1];
sx q[1];
rz(1.4820815) q[1];
x q[2];
rz(-1.734415) q[3];
sx q[3];
rz(-0.67577261) q[3];
sx q[3];
rz(2.8676652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37614432) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(-1.4260028) q[2];
rz(2.5358893) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894492) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(-0.28513518) q[0];
rz(0.51672283) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(1.8018988) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8782608) q[0];
sx q[0];
rz(-0.86255951) q[0];
sx q[0];
rz(1.0017298) q[0];
rz(1.2013024) q[2];
sx q[2];
rz(-0.72417799) q[2];
sx q[2];
rz(-0.61831123) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36741751) q[1];
sx q[1];
rz(-1.7304825) q[1];
sx q[1];
rz(1.3819441) q[1];
x q[2];
rz(-2.6684746) q[3];
sx q[3];
rz(-0.74672532) q[3];
sx q[3];
rz(2.8603922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67119917) q[2];
sx q[2];
rz(-0.931804) q[2];
sx q[2];
rz(-1.408067) q[2];
rz(-2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(-0.071921913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6614439) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(-2.9273422) q[0];
rz(2.7125773) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(2.4750211) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4144856) q[0];
sx q[0];
rz(-1.7761199) q[0];
sx q[0];
rz(1.7617102) q[0];
x q[1];
rz(1.2810983) q[2];
sx q[2];
rz(-0.87756598) q[2];
sx q[2];
rz(-1.2115657) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2506634) q[1];
sx q[1];
rz(-0.86859413) q[1];
sx q[1];
rz(2.3198747) q[1];
x q[2];
rz(0.14931071) q[3];
sx q[3];
rz(-1.4243323) q[3];
sx q[3];
rz(-3.1130476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.2864705) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(0.66037035) q[2];
rz(1.9541698) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(-2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289571) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(2.3045325) q[0];
rz(0.95056668) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(-1.9821092) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75393049) q[0];
sx q[0];
rz(-1.2571063) q[0];
sx q[0];
rz(2.4051106) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1690364) q[2];
sx q[2];
rz(-1.6998569) q[2];
sx q[2];
rz(1.3364524) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0598462) q[1];
sx q[1];
rz(-2.1085599) q[1];
sx q[1];
rz(-2.6408225) q[1];
rz(-pi) q[2];
rz(-0.67622185) q[3];
sx q[3];
rz(-2.3822504) q[3];
sx q[3];
rz(-0.29952213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7065457) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(1.8699899) q[2];
rz(2.0909677) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(-3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5746675) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(-1.1151474) q[0];
rz(-3.0976345) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(-0.10087068) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6504165) q[0];
sx q[0];
rz(-2.4917779) q[0];
sx q[0];
rz(2.8894436) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51668824) q[2];
sx q[2];
rz(-0.54737216) q[2];
sx q[2];
rz(1.237243) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.069177376) q[1];
sx q[1];
rz(-0.25670708) q[1];
sx q[1];
rz(-0.91255811) q[1];
x q[2];
rz(0.50562596) q[3];
sx q[3];
rz(-0.45590948) q[3];
sx q[3];
rz(-2.0924007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4042525) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(2.2163088) q[2];
rz(-2.960079) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7589384) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(2.1784901) q[0];
rz(-1.5513647) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(0.68774736) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15935414) q[0];
sx q[0];
rz(-1.237545) q[0];
sx q[0];
rz(-2.9100145) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1805448) q[2];
sx q[2];
rz(-0.38182048) q[2];
sx q[2];
rz(-2.080999) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2076599) q[1];
sx q[1];
rz(-1.7001517) q[1];
sx q[1];
rz(-1.5473066) q[1];
rz(-pi) q[2];
rz(-1.5958105) q[3];
sx q[3];
rz(-2.5731312) q[3];
sx q[3];
rz(-0.75077885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8035651) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(-2.2040099) q[2];
rz(-2.1607416) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6136318) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(-3.0623073) q[0];
rz(-1.1212564) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(-1.942873) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51533651) q[0];
sx q[0];
rz(-1.6394098) q[0];
sx q[0];
rz(-2.7068044) q[0];
rz(3.0555658) q[2];
sx q[2];
rz(-1.2089529) q[2];
sx q[2];
rz(-1.8334243) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.63003507) q[1];
sx q[1];
rz(-0.50543284) q[1];
sx q[1];
rz(-0.10314421) q[1];
x q[2];
rz(2.6506181) q[3];
sx q[3];
rz(-2.6262865) q[3];
sx q[3];
rz(1.5332424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85598677) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(-2.3525227) q[2];
rz(0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(-1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.566074) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(-2.5532706) q[0];
rz(3.1148124) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(2.2081597) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.371884) q[0];
sx q[0];
rz(-3.0196307) q[0];
sx q[0];
rz(2.6321649) q[0];
rz(1.9503715) q[2];
sx q[2];
rz(-1.2588725) q[2];
sx q[2];
rz(2.1385857) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.55798462) q[1];
sx q[1];
rz(-0.68478157) q[1];
sx q[1];
rz(0.058137356) q[1];
x q[2];
rz(2.4890355) q[3];
sx q[3];
rz(-1.2823294) q[3];
sx q[3];
rz(-1.0126589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0788706) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(-0.36994568) q[2];
rz(-2.2423819) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9545492) q[0];
sx q[0];
rz(-1.3267013) q[0];
sx q[0];
rz(-0.6533587) q[0];
rz(-1.6409142) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(-0.18383372) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18098022) q[0];
sx q[0];
rz(-1.4502589) q[0];
sx q[0];
rz(1.4535849) q[0];
x q[1];
rz(-2.3887563) q[2];
sx q[2];
rz(-0.92217731) q[2];
sx q[2];
rz(2.2019049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.93563423) q[1];
sx q[1];
rz(-1.7365672) q[1];
sx q[1];
rz(-0.88840719) q[1];
rz(-pi) q[2];
rz(-3.1179908) q[3];
sx q[3];
rz(-2.5581103) q[3];
sx q[3];
rz(-2.6445146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9188149) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(-2.9392021) q[2];
rz(-2.0683794) q[3];
sx q[3];
rz(-1.2349293) q[3];
sx q[3];
rz(-1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4595173) q[0];
sx q[0];
rz(-2.2820602) q[0];
sx q[0];
rz(2.5862502) q[0];
rz(0.36322414) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(-1.1584875) q[2];
sx q[2];
rz(-1.9956797) q[2];
sx q[2];
rz(3.0394625) q[2];
rz(0.98060676) q[3];
sx q[3];
rz(-1.7173613) q[3];
sx q[3];
rz(2.2955256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
