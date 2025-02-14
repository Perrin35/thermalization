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
rz(3.0124445) q[0];
sx q[0];
rz(-1.4718066) q[0];
sx q[0];
rz(2.1762525) q[0];
rz(-3.1211634) q[1];
sx q[1];
rz(-1.6579707) q[1];
sx q[1];
rz(1.3774011) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55104461) q[0];
sx q[0];
rz(-2.2529896) q[0];
sx q[0];
rz(0.57063535) q[0];
x q[1];
rz(-0.92490478) q[2];
sx q[2];
rz(-1.5965624) q[2];
sx q[2];
rz(2.5534867) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1024485) q[1];
sx q[1];
rz(-1.3833725) q[1];
sx q[1];
rz(1.4226807) q[1];
rz(-pi) q[2];
rz(0.80111472) q[3];
sx q[3];
rz(-1.5571888) q[3];
sx q[3];
rz(1.0281365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.91780245) q[2];
sx q[2];
rz(-1.4897572) q[2];
sx q[2];
rz(1.6415143) q[2];
rz(-2.3598059) q[3];
sx q[3];
rz(-0.89038554) q[3];
sx q[3];
rz(0.75209832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016945275) q[0];
sx q[0];
rz(-0.77777672) q[0];
sx q[0];
rz(0.97050226) q[0];
rz(-1.3264725) q[1];
sx q[1];
rz(-1.7992203) q[1];
sx q[1];
rz(-0.91189799) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7528943) q[0];
sx q[0];
rz(-1.7363743) q[0];
sx q[0];
rz(-2.5491211) q[0];
x q[1];
rz(0.17120338) q[2];
sx q[2];
rz(-1.7641626) q[2];
sx q[2];
rz(0.69583508) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7378547) q[1];
sx q[1];
rz(-2.4399839) q[1];
sx q[1];
rz(1.6340294) q[1];
x q[2];
rz(-1.6025869) q[3];
sx q[3];
rz(-0.45511757) q[3];
sx q[3];
rz(-1.8546752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9087002) q[2];
sx q[2];
rz(-1.1129271) q[2];
sx q[2];
rz(0.98928893) q[2];
rz(-0.42012897) q[3];
sx q[3];
rz(-2.1412854) q[3];
sx q[3];
rz(-0.72107983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99783889) q[0];
sx q[0];
rz(-1.1622575) q[0];
sx q[0];
rz(-1.0379399) q[0];
rz(-2.2833917) q[1];
sx q[1];
rz(-0.52640262) q[1];
sx q[1];
rz(0.16955489) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9301382) q[0];
sx q[0];
rz(-1.4272604) q[0];
sx q[0];
rz(-3.1375196) q[0];
x q[1];
rz(-0.900678) q[2];
sx q[2];
rz(-2.6478516) q[2];
sx q[2];
rz(-2.7639219) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.819014) q[1];
sx q[1];
rz(-0.59016229) q[1];
sx q[1];
rz(2.1165089) q[1];
rz(1.1812594) q[3];
sx q[3];
rz(-2.7523605) q[3];
sx q[3];
rz(1.8336982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9618591) q[2];
sx q[2];
rz(-2.6229975) q[2];
sx q[2];
rz(2.7355984) q[2];
rz(2.3640442) q[3];
sx q[3];
rz(-0.33027875) q[3];
sx q[3];
rz(0.8146666) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6430214) q[0];
sx q[0];
rz(-1.2914456) q[0];
sx q[0];
rz(0.92798573) q[0];
rz(3.0827177) q[1];
sx q[1];
rz(-1.6935655) q[1];
sx q[1];
rz(-1.7108797) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74364036) q[0];
sx q[0];
rz(-0.40983202) q[0];
sx q[0];
rz(-1.7663971) q[0];
rz(-pi) q[1];
rz(0.36305444) q[2];
sx q[2];
rz(-1.3842693) q[2];
sx q[2];
rz(-2.8721953) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9560066) q[1];
sx q[1];
rz(-1.7934467) q[1];
sx q[1];
rz(-2.0757448) q[1];
x q[2];
rz(0.9417504) q[3];
sx q[3];
rz(-1.3247196) q[3];
sx q[3];
rz(-1.6501282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6835988) q[2];
sx q[2];
rz(-1.6432089) q[2];
sx q[2];
rz(-1.9963473) q[2];
rz(2.2942719) q[3];
sx q[3];
rz(-1.8073795) q[3];
sx q[3];
rz(1.639521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97750807) q[0];
sx q[0];
rz(-1.9590398) q[0];
sx q[0];
rz(-2.7567647) q[0];
rz(-2.341914) q[1];
sx q[1];
rz(-2.622602) q[1];
sx q[1];
rz(1.5934561) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037796) q[0];
sx q[0];
rz(-1.8124969) q[0];
sx q[0];
rz(1.3194618) q[0];
rz(-1.9673011) q[2];
sx q[2];
rz(-1.7616211) q[2];
sx q[2];
rz(0.79939524) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0210798) q[1];
sx q[1];
rz(-1.2215633) q[1];
sx q[1];
rz(1.9053188) q[1];
rz(2.1206854) q[3];
sx q[3];
rz(-1.8916777) q[3];
sx q[3];
rz(-0.73559092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3758214) q[2];
sx q[2];
rz(-0.69362005) q[2];
sx q[2];
rz(-3.0613464) q[2];
rz(-2.5806228) q[3];
sx q[3];
rz(-1.6601945) q[3];
sx q[3];
rz(-2.5188353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4866667) q[0];
sx q[0];
rz(-0.66513649) q[0];
sx q[0];
rz(1.2605865) q[0];
rz(0.27443019) q[1];
sx q[1];
rz(-1.6703037) q[1];
sx q[1];
rz(2.4226277) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.052375) q[0];
sx q[0];
rz(-1.2219273) q[0];
sx q[0];
rz(0.87651663) q[0];
x q[1];
rz(-2.3173213) q[2];
sx q[2];
rz(-1.5780515) q[2];
sx q[2];
rz(-1.5430918) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8985892) q[1];
sx q[1];
rz(-0.92783992) q[1];
sx q[1];
rz(1.7220294) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0522862) q[3];
sx q[3];
rz(-0.42123367) q[3];
sx q[3];
rz(-0.13430922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4869953) q[2];
sx q[2];
rz(-1.2664814) q[2];
sx q[2];
rz(2.5385762) q[2];
rz(-1.4740137) q[3];
sx q[3];
rz(-1.0242198) q[3];
sx q[3];
rz(-0.41072861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6239768) q[0];
sx q[0];
rz(-1.5331974) q[0];
sx q[0];
rz(0.18727592) q[0];
rz(0.97224832) q[1];
sx q[1];
rz(-2.9863803) q[1];
sx q[1];
rz(0.42047277) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0197821) q[0];
sx q[0];
rz(-2.3276637) q[0];
sx q[0];
rz(2.5518083) q[0];
x q[1];
rz(-2.1537909) q[2];
sx q[2];
rz(-0.62623411) q[2];
sx q[2];
rz(-0.69863182) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5198316) q[1];
sx q[1];
rz(-0.48096195) q[1];
sx q[1];
rz(2.2051808) q[1];
rz(1.2240846) q[3];
sx q[3];
rz(-0.60255749) q[3];
sx q[3];
rz(-1.4977221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.095852701) q[2];
sx q[2];
rz(-2.0407712) q[2];
sx q[2];
rz(-0.25071684) q[2];
rz(-1.2480674) q[3];
sx q[3];
rz(-1.7496795) q[3];
sx q[3];
rz(-0.59984508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8641758) q[0];
sx q[0];
rz(-0.59159788) q[0];
sx q[0];
rz(-2.199882) q[0];
rz(0.038453728) q[1];
sx q[1];
rz(-1.4310623) q[1];
sx q[1];
rz(-3.0252735) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4394035) q[0];
sx q[0];
rz(-0.84745126) q[0];
sx q[0];
rz(-0.96203104) q[0];
rz(-pi) q[1];
rz(1.3384678) q[2];
sx q[2];
rz(-0.69455494) q[2];
sx q[2];
rz(-1.6603927) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1585576) q[1];
sx q[1];
rz(-1.8280892) q[1];
sx q[1];
rz(1.757181) q[1];
x q[2];
rz(2.8781901) q[3];
sx q[3];
rz(-2.7399094) q[3];
sx q[3];
rz(-2.4954737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2945127) q[2];
sx q[2];
rz(-2.3307762) q[2];
sx q[2];
rz(1.026356) q[2];
rz(-1.9319084) q[3];
sx q[3];
rz(-1.0299725) q[3];
sx q[3];
rz(0.9616372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.472979) q[0];
sx q[0];
rz(-1.0740148) q[0];
sx q[0];
rz(2.7630254) q[0];
rz(2.0319132) q[1];
sx q[1];
rz(-0.30866426) q[1];
sx q[1];
rz(3.1139156) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8016114) q[0];
sx q[0];
rz(-1.2072354) q[0];
sx q[0];
rz(-0.63940639) q[0];
rz(-1.7871066) q[2];
sx q[2];
rz(-1.7515108) q[2];
sx q[2];
rz(1.9150776) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82904774) q[1];
sx q[1];
rz(-1.2133737) q[1];
sx q[1];
rz(-0.19695671) q[1];
rz(-pi) q[2];
rz(0.31881551) q[3];
sx q[3];
rz(-1.1622815) q[3];
sx q[3];
rz(-0.2948979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8370342) q[2];
sx q[2];
rz(-0.95418945) q[2];
sx q[2];
rz(-2.7216116) q[2];
rz(-2.3000681) q[3];
sx q[3];
rz(-2.1244815) q[3];
sx q[3];
rz(2.7628472) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751752) q[0];
sx q[0];
rz(-3.0247122) q[0];
sx q[0];
rz(2.8564659) q[0];
rz(2.9410703) q[1];
sx q[1];
rz(-1.2031809) q[1];
sx q[1];
rz(2.3060422) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052976089) q[0];
sx q[0];
rz(-2.0248013) q[0];
sx q[0];
rz(0.80843057) q[0];
x q[1];
rz(-1.5642688) q[2];
sx q[2];
rz(-0.92374974) q[2];
sx q[2];
rz(-1.4359635) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0613868) q[1];
sx q[1];
rz(-2.5701984) q[1];
sx q[1];
rz(-0.6462884) q[1];
x q[2];
rz(0.16160139) q[3];
sx q[3];
rz(-0.23862632) q[3];
sx q[3];
rz(1.0922701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7150813) q[2];
sx q[2];
rz(-1.0914404) q[2];
sx q[2];
rz(2.4439028) q[2];
rz(0.52608025) q[3];
sx q[3];
rz(-0.79280058) q[3];
sx q[3];
rz(1.5066159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0195011) q[0];
sx q[0];
rz(-2.2631336) q[0];
sx q[0];
rz(2.7418131) q[0];
rz(1.9922235) q[1];
sx q[1];
rz(-0.72723564) q[1];
sx q[1];
rz(-0.74692187) q[1];
rz(-1.9952075) q[2];
sx q[2];
rz(-0.72327264) q[2];
sx q[2];
rz(-2.8930152) q[2];
rz(1.2739194) q[3];
sx q[3];
rz(-2.0991201) q[3];
sx q[3];
rz(0.42602628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
