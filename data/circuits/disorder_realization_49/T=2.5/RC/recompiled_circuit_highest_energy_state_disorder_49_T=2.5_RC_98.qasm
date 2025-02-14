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
rz(-3.1168923) q[0];
sx q[0];
rz(-1.8827117) q[0];
sx q[0];
rz(-1.2727241) q[0];
rz(-4.89115) q[1];
sx q[1];
rz(1.7311544) q[1];
sx q[1];
rz(16.601736) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0426038) q[0];
sx q[0];
rz(-2.1402142) q[0];
sx q[0];
rz(0.75890394) q[0];
rz(2.4364901) q[2];
sx q[2];
rz(-2.7573708) q[2];
sx q[2];
rz(0.23300685) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.53613816) q[1];
sx q[1];
rz(-2.1011323) q[1];
sx q[1];
rz(1.0592036) q[1];
x q[2];
rz(2.7734766) q[3];
sx q[3];
rz(-1.5541346) q[3];
sx q[3];
rz(-2.1429495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7749403) q[2];
sx q[2];
rz(-1.6795936) q[2];
sx q[2];
rz(-2.6095663) q[2];
rz(-2.4074647) q[3];
sx q[3];
rz(-2.8668154) q[3];
sx q[3];
rz(-1.1067357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4183913) q[0];
sx q[0];
rz(-0.98576236) q[0];
sx q[0];
rz(0.080408737) q[0];
rz(0.53781992) q[1];
sx q[1];
rz(-2.0729013) q[1];
sx q[1];
rz(0.72371662) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0590032) q[0];
sx q[0];
rz(-1.132611) q[0];
sx q[0];
rz(0.74371145) q[0];
rz(-pi) q[1];
rz(-0.56494464) q[2];
sx q[2];
rz(-0.98472847) q[2];
sx q[2];
rz(-0.63462574) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.5840877) q[1];
sx q[1];
rz(-0.98230108) q[1];
sx q[1];
rz(3.1397538) q[1];
x q[2];
rz(0.81775093) q[3];
sx q[3];
rz(-1.1766542) q[3];
sx q[3];
rz(0.5241636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85372743) q[2];
sx q[2];
rz(-2.008805) q[2];
sx q[2];
rz(2.2354324) q[2];
rz(0.3624889) q[3];
sx q[3];
rz(-1.5443085) q[3];
sx q[3];
rz(3.0947065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2694117) q[0];
sx q[0];
rz(-1.325664) q[0];
sx q[0];
rz(-2.0500702) q[0];
rz(0.12256924) q[1];
sx q[1];
rz(-1.4361607) q[1];
sx q[1];
rz(-1.3465808) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5751981) q[0];
sx q[0];
rz(-1.9267531) q[0];
sx q[0];
rz(0.031513799) q[0];
x q[1];
rz(-2.4819751) q[2];
sx q[2];
rz(-1.3085367) q[2];
sx q[2];
rz(-0.66103092) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0071738) q[1];
sx q[1];
rz(-1.6225623) q[1];
sx q[1];
rz(1.365713) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19119672) q[3];
sx q[3];
rz(-0.53732291) q[3];
sx q[3];
rz(0.69512284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4284105) q[2];
sx q[2];
rz(-2.0117663) q[2];
sx q[2];
rz(-0.23207363) q[2];
rz(-1.170916) q[3];
sx q[3];
rz(-2.5210095) q[3];
sx q[3];
rz(2.8569729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21122268) q[0];
sx q[0];
rz(-2.3486597) q[0];
sx q[0];
rz(-1.6388182) q[0];
rz(-0.59208313) q[1];
sx q[1];
rz(-1.9860257) q[1];
sx q[1];
rz(2.2519462) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1680206) q[0];
sx q[0];
rz(-2.4102927) q[0];
sx q[0];
rz(2.2941053) q[0];
rz(2.8977029) q[2];
sx q[2];
rz(-2.0497344) q[2];
sx q[2];
rz(2.6289491) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0977323) q[1];
sx q[1];
rz(-1.3907038) q[1];
sx q[1];
rz(-1.8921683) q[1];
rz(-pi) q[2];
rz(-1.9806421) q[3];
sx q[3];
rz(-1.92627) q[3];
sx q[3];
rz(2.6419169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4952116) q[2];
sx q[2];
rz(-1.6782328) q[2];
sx q[2];
rz(2.4656673) q[2];
rz(-1.4365139) q[3];
sx q[3];
rz(-0.33994514) q[3];
sx q[3];
rz(2.9743312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9306358) q[0];
sx q[0];
rz(-2.8093331) q[0];
sx q[0];
rz(-2.9462872) q[0];
rz(0.12061128) q[1];
sx q[1];
rz(-0.21574012) q[1];
sx q[1];
rz(0.19048555) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99283907) q[0];
sx q[0];
rz(-1.6354113) q[0];
sx q[0];
rz(3.0535327) q[0];
rz(-pi) q[1];
rz(-0.10730524) q[2];
sx q[2];
rz(-2.4621747) q[2];
sx q[2];
rz(0.36919644) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.021564158) q[1];
sx q[1];
rz(-0.11332527) q[1];
sx q[1];
rz(0.47685195) q[1];
rz(-pi) q[2];
rz(-0.56511527) q[3];
sx q[3];
rz(-2.4451974) q[3];
sx q[3];
rz(2.2415862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5087937) q[2];
sx q[2];
rz(-1.6570647) q[2];
sx q[2];
rz(1.9602027) q[2];
rz(1.3440291) q[3];
sx q[3];
rz(-2.400178) q[3];
sx q[3];
rz(-2.369829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8089499) q[0];
sx q[0];
rz(-1.7985666) q[0];
sx q[0];
rz(1.1395662) q[0];
rz(-0.66967669) q[1];
sx q[1];
rz(-2.2256336) q[1];
sx q[1];
rz(-2.0119827) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78167262) q[0];
sx q[0];
rz(-1.8501213) q[0];
sx q[0];
rz(1.8769699) q[0];
rz(-2.4239842) q[2];
sx q[2];
rz(-1.773218) q[2];
sx q[2];
rz(-0.69684579) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.672513) q[1];
sx q[1];
rz(-1.7653078) q[1];
sx q[1];
rz(-2.7698293) q[1];
rz(-2.4242006) q[3];
sx q[3];
rz(-0.23788568) q[3];
sx q[3];
rz(2.5715076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.42664042) q[2];
sx q[2];
rz(-1.4098097) q[2];
sx q[2];
rz(-2.7362774) q[2];
rz(-2.3146368) q[3];
sx q[3];
rz(-2.5155641) q[3];
sx q[3];
rz(0.85957447) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7760794) q[0];
sx q[0];
rz(-3.1190393) q[0];
sx q[0];
rz(0.93589163) q[0];
rz(2.2731958) q[1];
sx q[1];
rz(-0.52803841) q[1];
sx q[1];
rz(-1.7220928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7129684) q[0];
sx q[0];
rz(-1.0626864) q[0];
sx q[0];
rz(-1.1626121) q[0];
x q[1];
rz(-1.9738036) q[2];
sx q[2];
rz(-2.5984077) q[2];
sx q[2];
rz(2.3576971) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8384435) q[1];
sx q[1];
rz(-1.1080964) q[1];
sx q[1];
rz(-0.82973231) q[1];
rz(2.9299632) q[3];
sx q[3];
rz(-3.1397925) q[3];
sx q[3];
rz(1.7626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.821227) q[2];
sx q[2];
rz(-1.7419523) q[2];
sx q[2];
rz(-1.0199176) q[2];
rz(-2.6323281) q[3];
sx q[3];
rz(-1.441322) q[3];
sx q[3];
rz(-3.063108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34198636) q[0];
sx q[0];
rz(-1.5155563) q[0];
sx q[0];
rz(-1.1588143) q[0];
rz(-0.47053567) q[1];
sx q[1];
rz(-1.7262986) q[1];
sx q[1];
rz(-1.6758957) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5138595) q[0];
sx q[0];
rz(-3.0854221) q[0];
sx q[0];
rz(-0.43561952) q[0];
rz(-pi) q[1];
rz(1.2242555) q[2];
sx q[2];
rz(-0.90991352) q[2];
sx q[2];
rz(-0.61049547) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.640721) q[1];
sx q[1];
rz(-1.6476008) q[1];
sx q[1];
rz(2.4810664) q[1];
rz(-0.37181446) q[3];
sx q[3];
rz(-0.68228693) q[3];
sx q[3];
rz(2.2390847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9965685) q[2];
sx q[2];
rz(-1.327876) q[2];
sx q[2];
rz(2.3222951) q[2];
rz(0.79646349) q[3];
sx q[3];
rz(-2.7345246) q[3];
sx q[3];
rz(-1.6527294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2040937) q[0];
sx q[0];
rz(-3.0452073) q[0];
sx q[0];
rz(-0.82164422) q[0];
rz(-2.4687528) q[1];
sx q[1];
rz(-2.5655589) q[1];
sx q[1];
rz(2.0988665) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48548904) q[0];
sx q[0];
rz(-0.21767958) q[0];
sx q[0];
rz(-0.34164048) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94093948) q[2];
sx q[2];
rz(-0.90412882) q[2];
sx q[2];
rz(-0.43064865) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4303226) q[1];
sx q[1];
rz(-1.9414328) q[1];
sx q[1];
rz(2.1881585) q[1];
rz(-0.25298499) q[3];
sx q[3];
rz(-1.7182941) q[3];
sx q[3];
rz(1.1563039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3599856) q[2];
sx q[2];
rz(-1.2722509) q[2];
sx q[2];
rz(-0.84235111) q[2];
rz(0.53884566) q[3];
sx q[3];
rz(-1.6539961) q[3];
sx q[3];
rz(2.6174788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6738324) q[0];
sx q[0];
rz(-2.9957275) q[0];
sx q[0];
rz(1.9061331) q[0];
rz(-2.1927059) q[1];
sx q[1];
rz(-0.83273879) q[1];
sx q[1];
rz(1.6611151) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4119856) q[0];
sx q[0];
rz(-0.74568891) q[0];
sx q[0];
rz(-1.1545769) q[0];
rz(-pi) q[1];
rz(-0.5625426) q[2];
sx q[2];
rz(-1.0357366) q[2];
sx q[2];
rz(-1.8059274) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2044511) q[1];
sx q[1];
rz(-1.6455263) q[1];
sx q[1];
rz(0.16507574) q[1];
rz(-pi) q[2];
rz(-2.5148095) q[3];
sx q[3];
rz(-2.7434182) q[3];
sx q[3];
rz(-2.6836723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.71198717) q[2];
sx q[2];
rz(-1.4236071) q[2];
sx q[2];
rz(-0.31022662) q[2];
rz(2.6817536) q[3];
sx q[3];
rz(-0.55791563) q[3];
sx q[3];
rz(1.7673813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9643758) q[0];
sx q[0];
rz(-1.1342659) q[0];
sx q[0];
rz(-2.076617) q[0];
rz(0.89827697) q[1];
sx q[1];
rz(-2.5986462) q[1];
sx q[1];
rz(2.6959261) q[1];
rz(-1.0679792) q[2];
sx q[2];
rz(-2.4120654) q[2];
sx q[2];
rz(0.21929693) q[2];
rz(1.0887922) q[3];
sx q[3];
rz(-2.4203151) q[3];
sx q[3];
rz(-2.388777) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
