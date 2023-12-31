OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.1383837) q[0];
sx q[0];
rz(-2.9870343) q[0];
sx q[0];
rz(2.4490693) q[0];
rz(1.9321631) q[1];
sx q[1];
rz(-1.2485319) q[1];
sx q[1];
rz(-1.385153) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6801075) q[0];
sx q[0];
rz(-0.84851096) q[0];
sx q[0];
rz(-1.1397584) q[0];
rz(-1.295747) q[2];
sx q[2];
rz(-1.0359456) q[2];
sx q[2];
rz(-1.5799074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17978046) q[1];
sx q[1];
rz(-1.7331859) q[1];
sx q[1];
rz(0.23602545) q[1];
x q[2];
rz(2.5296506) q[3];
sx q[3];
rz(-0.7512593) q[3];
sx q[3];
rz(0.9179759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2549071) q[2];
sx q[2];
rz(-0.79780769) q[2];
sx q[2];
rz(2.936426) q[2];
rz(0.77130476) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(-1.1024968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40760621) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(2.6876887) q[0];
rz(-2.1167963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(-1.9143547) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42392143) q[0];
sx q[0];
rz(-0.14980355) q[0];
sx q[0];
rz(1.040209) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51867698) q[2];
sx q[2];
rz(-1.1653324) q[2];
sx q[2];
rz(2.5395218) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2482359) q[1];
sx q[1];
rz(-1.8289939) q[1];
sx q[1];
rz(2.8398819) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9407528) q[3];
sx q[3];
rz(-1.4174995) q[3];
sx q[3];
rz(-2.4917345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1002905) q[2];
sx q[2];
rz(-1.1854478) q[2];
sx q[2];
rz(-0.56742898) q[2];
rz(-0.36519095) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(-0.96810961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48297468) q[0];
sx q[0];
rz(-0.56476074) q[0];
sx q[0];
rz(0.89865249) q[0];
rz(-2.1458416) q[1];
sx q[1];
rz(-1.5834705) q[1];
sx q[1];
rz(-0.333289) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8329187) q[0];
sx q[0];
rz(-1.1090288) q[0];
sx q[0];
rz(-0.69899107) q[0];
x q[1];
rz(-2.0731508) q[2];
sx q[2];
rz(-0.2012673) q[2];
sx q[2];
rz(-1.4301436) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7283199) q[1];
sx q[1];
rz(-1.4188758) q[1];
sx q[1];
rz(2.1876213) q[1];
rz(-pi) q[2];
rz(2.1214478) q[3];
sx q[3];
rz(-1.9198951) q[3];
sx q[3];
rz(-2.246644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4553392) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(-1.1509482) q[2];
rz(-0.84093705) q[3];
sx q[3];
rz(-1.1154113) q[3];
sx q[3];
rz(1.1545198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6999917) q[0];
sx q[0];
rz(-1.5804407) q[0];
sx q[0];
rz(0.68471318) q[0];
rz(-2.1060064) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(-1.9365786) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9158463) q[0];
sx q[0];
rz(-1.5374743) q[0];
sx q[0];
rz(-2.2621364) q[0];
rz(-pi) q[1];
rz(2.2150546) q[2];
sx q[2];
rz(-1.4069923) q[2];
sx q[2];
rz(2.0870199) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11338621) q[1];
sx q[1];
rz(-0.75062597) q[1];
sx q[1];
rz(-1.9764465) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0699176) q[3];
sx q[3];
rz(-2.274401) q[3];
sx q[3];
rz(-3.049831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8923607) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(-0.37115804) q[2];
rz(-1.4012198) q[3];
sx q[3];
rz(-0.6597844) q[3];
sx q[3];
rz(-2.0223117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.086833) q[0];
sx q[0];
rz(-2.355447) q[0];
sx q[0];
rz(-3.0084685) q[0];
rz(0.99331028) q[1];
sx q[1];
rz(-1.7555833) q[1];
sx q[1];
rz(-2.5865119) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18300444) q[0];
sx q[0];
rz(-0.31561139) q[0];
sx q[0];
rz(-1.6118227) q[0];
rz(-pi) q[1];
x q[1];
rz(1.416989) q[2];
sx q[2];
rz(-2.7222735) q[2];
sx q[2];
rz(-0.16659444) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3302147) q[1];
sx q[1];
rz(-1.1256309) q[1];
sx q[1];
rz(3.1203169) q[1];
x q[2];
rz(-2.5693232) q[3];
sx q[3];
rz(-0.096147691) q[3];
sx q[3];
rz(-0.26086807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.30620265) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(-3.0026657) q[2];
rz(0.94240087) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(-0.55148235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54365629) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(2.561835) q[0];
rz(-0.12750164) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(1.5396083) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67028763) q[0];
sx q[0];
rz(-1.3677214) q[0];
sx q[0];
rz(-0.85104403) q[0];
rz(-pi) q[1];
rz(-0.13055735) q[2];
sx q[2];
rz(-1.6160384) q[2];
sx q[2];
rz(3.0322078) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.786799) q[1];
sx q[1];
rz(-1.8537632) q[1];
sx q[1];
rz(-2.016504) q[1];
rz(2.2555389) q[3];
sx q[3];
rz(-1.7985117) q[3];
sx q[3];
rz(-0.39623228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.55398983) q[2];
sx q[2];
rz(-2.8911399) q[2];
sx q[2];
rz(0.26947752) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-2.6168489) q[3];
sx q[3];
rz(0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6948029) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(2.457298) q[0];
rz(0.11958312) q[1];
sx q[1];
rz(-1.8493098) q[1];
sx q[1];
rz(2.6228242) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25012384) q[0];
sx q[0];
rz(-0.84202535) q[0];
sx q[0];
rz(-2.9617589) q[0];
rz(0.22612818) q[2];
sx q[2];
rz(-0.78352189) q[2];
sx q[2];
rz(2.4353611) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.24592933) q[1];
sx q[1];
rz(-0.83918011) q[1];
sx q[1];
rz(2.500446) q[1];
x q[2];
rz(1.6325083) q[3];
sx q[3];
rz(-0.40137526) q[3];
sx q[3];
rz(2.9369831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8873022) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(-2.5781412) q[2];
rz(3.0900132) q[3];
sx q[3];
rz(-2.0023465) q[3];
sx q[3];
rz(-2.1896867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.96034399) q[0];
sx q[0];
rz(-0.5287756) q[0];
sx q[0];
rz(1.3990336) q[0];
rz(0.78701204) q[1];
sx q[1];
rz(-1.1287289) q[1];
sx q[1];
rz(0.74434892) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3314914) q[0];
sx q[0];
rz(-1.6252675) q[0];
sx q[0];
rz(-2.9936552) q[0];
rz(-pi) q[1];
rz(1.8366351) q[2];
sx q[2];
rz(-2.611534) q[2];
sx q[2];
rz(0.30345464) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62521711) q[1];
sx q[1];
rz(-2.5853734) q[1];
sx q[1];
rz(-0.43362995) q[1];
rz(-pi) q[2];
rz(-2.6870071) q[3];
sx q[3];
rz(-1.5919519) q[3];
sx q[3];
rz(-1.3190312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7156334) q[2];
sx q[2];
rz(-1.3954433) q[2];
sx q[2];
rz(-2.5320833) q[2];
rz(0.65731796) q[3];
sx q[3];
rz(-2.4980563) q[3];
sx q[3];
rz(-2.8801584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9534) q[0];
sx q[0];
rz(-3.0472026) q[0];
sx q[0];
rz(1.5040065) q[0];
rz(1.2414744) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(0.77493587) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0749045) q[0];
sx q[0];
rz(-2.3432891) q[0];
sx q[0];
rz(-1.0134646) q[0];
rz(-1.1275034) q[2];
sx q[2];
rz(-1.4294251) q[2];
sx q[2];
rz(3.1183426) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96156582) q[1];
sx q[1];
rz(-2.0225836) q[1];
sx q[1];
rz(-1.3614484) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1433467) q[3];
sx q[3];
rz(-0.25902723) q[3];
sx q[3];
rz(1.6365901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.187414) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(0.15787086) q[2];
rz(-1.212451) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(-1.8635748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0697486) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(0.21433314) q[0];
rz(0.65746039) q[1];
sx q[1];
rz(-2.9174556) q[1];
sx q[1];
rz(-2.0956031) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25046529) q[0];
sx q[0];
rz(-1.2123101) q[0];
sx q[0];
rz(-1.6842708) q[0];
rz(-pi) q[1];
rz(3.0984512) q[2];
sx q[2];
rz(-0.59141814) q[2];
sx q[2];
rz(2.6864348) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6194832) q[1];
sx q[1];
rz(-1.2088747) q[1];
sx q[1];
rz(-2.1543909) q[1];
rz(-pi) q[2];
rz(-1.5442185) q[3];
sx q[3];
rz(-1.1655032) q[3];
sx q[3];
rz(-1.4881031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7252698) q[2];
sx q[2];
rz(-3.0815093) q[2];
sx q[2];
rz(1.0894758) q[2];
rz(-1.5754835) q[3];
sx q[3];
rz(-1.243306) q[3];
sx q[3];
rz(0.65264788) q[3];
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
rz(1.2789223) q[0];
sx q[0];
rz(-2.537732) q[0];
sx q[0];
rz(-2.296007) q[0];
rz(-1.6090341) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(2.5074742) q[2];
sx q[2];
rz(-2.6478883) q[2];
sx q[2];
rz(-1.4512856) q[2];
rz(-1.5076751) q[3];
sx q[3];
rz(-2.1199385) q[3];
sx q[3];
rz(-1.4169823) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
