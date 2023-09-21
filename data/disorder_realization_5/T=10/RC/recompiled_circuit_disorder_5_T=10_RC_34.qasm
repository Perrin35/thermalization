OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6025699) q[0];
sx q[0];
rz(-0.56350001) q[0];
sx q[0];
rz(-2.6846057) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(1.2759804) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2109283) q[0];
sx q[0];
rz(-1.7893447) q[0];
sx q[0];
rz(1.7095837) q[0];
rz(2.4495507) q[2];
sx q[2];
rz(-0.11728742) q[2];
sx q[2];
rz(0.28809822) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7282655) q[1];
sx q[1];
rz(-1.3682377) q[1];
sx q[1];
rz(0.061518016) q[1];
rz(-pi) q[2];
rz(2.7492417) q[3];
sx q[3];
rz(-2.1808743) q[3];
sx q[3];
rz(3.0685134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0322545) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(-0.16513744) q[2];
rz(1.6312381) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(2.3993649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363268) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(0.045036137) q[0];
rz(-0.33915195) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(-0.80274686) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4883603) q[0];
sx q[0];
rz(-0.45408861) q[0];
sx q[0];
rz(1.0351719) q[0];
rz(-3.0375508) q[2];
sx q[2];
rz(-1.0696971) q[2];
sx q[2];
rz(0.83545557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74233494) q[1];
sx q[1];
rz(-2.5207673) q[1];
sx q[1];
rz(1.4820815) q[1];
rz(-3.0117399) q[3];
sx q[3];
rz(-2.2359071) q[3];
sx q[3];
rz(0.065404281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7654483) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(1.7155898) q[2];
rz(0.60570335) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894492) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(0.28513518) q[0];
rz(-0.51672283) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(1.3396938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8782608) q[0];
sx q[0];
rz(-0.86255951) q[0];
sx q[0];
rz(1.0017298) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88111934) q[2];
sx q[2];
rz(-1.3291877) q[2];
sx q[2];
rz(1.9067665) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.50946402) q[1];
sx q[1];
rz(-0.2466991) q[1];
sx q[1];
rz(0.86172608) q[1];
rz(-pi) q[2];
rz(-1.1717103) q[3];
sx q[3];
rz(-0.92150021) q[3];
sx q[3];
rz(-2.2513575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67119917) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(-1.408067) q[2];
rz(2.9584598) q[3];
sx q[3];
rz(-1.0147084) q[3];
sx q[3];
rz(3.0696707) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6614439) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(-2.9273422) q[0];
rz(-0.42901531) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(-2.4750211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19569451) q[0];
sx q[0];
rz(-1.7576522) q[0];
sx q[0];
rz(0.20901434) q[0];
x q[1];
rz(2.8104086) q[2];
sx q[2];
rz(-2.3996155) q[2];
sx q[2];
rz(-0.77510288) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2506634) q[1];
sx q[1];
rz(-0.86859413) q[1];
sx q[1];
rz(-2.3198747) q[1];
x q[2];
rz(2.9922819) q[3];
sx q[3];
rz(-1.4243323) q[3];
sx q[3];
rz(-0.028545054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2864705) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(2.4812223) q[2];
rz(-1.9541698) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(-0.56046265) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289571) q[0];
sx q[0];
rz(-1.2197138) q[0];
sx q[0];
rz(0.83706013) q[0];
rz(-0.95056668) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(1.1594835) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1449074) q[0];
sx q[0];
rz(-0.7888182) q[0];
sx q[0];
rz(-2.6916654) q[0];
rz(-1.3443089) q[2];
sx q[2];
rz(-0.61033568) q[2];
sx q[2];
rz(-0.047705334) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9002478) q[1];
sx q[1];
rz(-0.7175788) q[1];
sx q[1];
rz(-2.248583) q[1];
rz(0.63727832) q[3];
sx q[3];
rz(-1.1253469) q[3];
sx q[3];
rz(-2.3973936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.435047) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(-1.2716028) q[2];
rz(-1.050625) q[3];
sx q[3];
rz(-1.2579066) q[3];
sx q[3];
rz(3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5746675) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(2.0264453) q[0];
rz(3.0976345) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(0.10087068) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2642919) q[0];
sx q[0];
rz(-1.4192686) q[0];
sx q[0];
rz(-0.63440462) q[0];
x q[1];
rz(-2.6249044) q[2];
sx q[2];
rz(-2.5942205) q[2];
sx q[2];
rz(-1.9043497) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2821591) q[1];
sx q[1];
rz(-1.4148501) q[1];
sx q[1];
rz(-1.7755416) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6359667) q[3];
sx q[3];
rz(-0.45590948) q[3];
sx q[3];
rz(-2.0924007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7373401) q[2];
sx q[2];
rz(-1.5161783) q[2];
sx q[2];
rz(0.92528382) q[2];
rz(0.18151367) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(-1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826542) q[0];
sx q[0];
rz(-2.721334) q[0];
sx q[0];
rz(-2.1784901) q[0];
rz(1.5902279) q[1];
sx q[1];
rz(-2.1236877) q[1];
sx q[1];
rz(2.4538453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6531658) q[0];
sx q[0];
rz(-1.7894206) q[0];
sx q[0];
rz(1.9125008) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1805448) q[2];
sx q[2];
rz(-0.38182048) q[2];
sx q[2];
rz(-2.080999) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5016985) q[1];
sx q[1];
rz(-1.5940897) q[1];
sx q[1];
rz(-3.0122019) q[1];
rz(-1.0024768) q[3];
sx q[3];
rz(-1.5842614) q[3];
sx q[3];
rz(-0.79893597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3380276) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(-0.93758279) q[2];
rz(2.1607416) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6136318) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(0.079285346) q[0];
rz(-2.0203363) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(1.942873) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51533651) q[0];
sx q[0];
rz(-1.6394098) q[0];
sx q[0];
rz(2.7068044) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2077246) q[2];
sx q[2];
rz(-1.4903526) q[2];
sx q[2];
rz(0.23210873) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3938013) q[1];
sx q[1];
rz(-1.0683021) q[1];
sx q[1];
rz(1.5138813) q[1];
rz(-pi) q[2];
rz(-2.6783887) q[3];
sx q[3];
rz(-1.8052881) q[3];
sx q[3];
rz(-2.6686058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85598677) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(-0.78906995) q[2];
rz(-0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.566074) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(-3.1148124) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(0.9334329) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.371884) q[0];
sx q[0];
rz(-0.121962) q[0];
sx q[0];
rz(2.6321649) q[0];
rz(-pi) q[1];
x q[1];
rz(2.80745) q[2];
sx q[2];
rz(-1.2103989) q[2];
sx q[2];
rz(0.44597086) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.583608) q[1];
sx q[1];
rz(-0.68478157) q[1];
sx q[1];
rz(0.058137356) q[1];
rz(0.45456072) q[3];
sx q[3];
rz(-2.436736) q[3];
sx q[3];
rz(0.2017894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.062722) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(2.771647) q[2];
rz(-2.2423819) q[3];
sx q[3];
rz(-1.6588541) q[3];
sx q[3];
rz(-0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.9545492) q[0];
sx q[0];
rz(-1.3267013) q[0];
sx q[0];
rz(-2.488234) q[0];
rz(1.6409142) q[1];
sx q[1];
rz(-1.6710284) q[1];
sx q[1];
rz(-0.18383372) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7659347) q[0];
sx q[0];
rz(-1.6871534) q[0];
sx q[0];
rz(-3.0202306) q[0];
rz(-pi) q[1];
rz(0.75283639) q[2];
sx q[2];
rz(-2.2194153) q[2];
sx q[2];
rz(0.93968771) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.4348236) q[1];
sx q[1];
rz(-0.69908792) q[1];
sx q[1];
rz(-1.8301151) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5582383) q[3];
sx q[3];
rz(-1.5577941) q[3];
sx q[3];
rz(1.0540203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2227778) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(0.20239057) q[2];
rz(2.0683794) q[3];
sx q[3];
rz(-1.2349293) q[3];
sx q[3];
rz(-1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.4595173) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(0.36322414) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(-2.4167378) q[2];
sx q[2];
rz(-0.58314322) q[2];
sx q[2];
rz(-2.4287139) q[2];
rz(-2.965747) q[3];
sx q[3];
rz(-2.1538215) q[3];
sx q[3];
rz(0.62721696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
