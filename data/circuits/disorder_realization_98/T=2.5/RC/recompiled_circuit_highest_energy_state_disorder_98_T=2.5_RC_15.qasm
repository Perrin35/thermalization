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
rz(-0.9653402) q[0];
rz(0.020429285) q[1];
sx q[1];
rz(-1.483622) q[1];
sx q[1];
rz(-1.3774011) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3453774) q[0];
sx q[0];
rz(-0.85897972) q[0];
sx q[0];
rz(0.9839566) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1093311) q[2];
sx q[2];
rz(-0.92515495) q[2];
sx q[2];
rz(-2.1394859) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1024485) q[1];
sx q[1];
rz(-1.7582201) q[1];
sx q[1];
rz(-1.718912) q[1];
rz(-pi) q[2];
rz(1.5512439) q[3];
sx q[3];
rz(-0.76977714) q[3];
sx q[3];
rz(2.584892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91780245) q[2];
sx q[2];
rz(-1.6518355) q[2];
sx q[2];
rz(-1.6415143) q[2];
rz(-2.3598059) q[3];
sx q[3];
rz(-0.89038554) q[3];
sx q[3];
rz(-2.3894943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1246474) q[0];
sx q[0];
rz(-2.3638159) q[0];
sx q[0];
rz(-0.97050226) q[0];
rz(-1.3264725) q[1];
sx q[1];
rz(-1.3423723) q[1];
sx q[1];
rz(0.91189799) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3886984) q[0];
sx q[0];
rz(-1.4052183) q[0];
sx q[0];
rz(2.5491211) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17120338) q[2];
sx q[2];
rz(-1.3774301) q[2];
sx q[2];
rz(-0.69583508) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9262097) q[1];
sx q[1];
rz(-1.5299986) q[1];
sx q[1];
rz(-2.271419) q[1];
rz(-1.5390057) q[3];
sx q[3];
rz(-2.6864751) q[3];
sx q[3];
rz(1.2869175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9087002) q[2];
sx q[2];
rz(-2.0286655) q[2];
sx q[2];
rz(-0.98928893) q[2];
rz(-0.42012897) q[3];
sx q[3];
rz(-1.0003072) q[3];
sx q[3];
rz(-2.4205128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1437538) q[0];
sx q[0];
rz(-1.1622575) q[0];
sx q[0];
rz(1.0379399) q[0];
rz(2.2833917) q[1];
sx q[1];
rz(-0.52640262) q[1];
sx q[1];
rz(-0.16955489) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9586048) q[0];
sx q[0];
rz(-0.14359328) q[0];
sx q[0];
rz(1.5426226) q[0];
rz(-1.9699668) q[2];
sx q[2];
rz(-1.2720275) q[2];
sx q[2];
rz(-0.58387652) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.21903164) q[1];
sx q[1];
rz(-1.2777849) q[1];
sx q[1];
rz(-1.0508363) q[1];
rz(-1.1812594) q[3];
sx q[3];
rz(-0.3892322) q[3];
sx q[3];
rz(-1.3078944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1797336) q[2];
sx q[2];
rz(-2.6229975) q[2];
sx q[2];
rz(2.7355984) q[2];
rz(2.3640442) q[3];
sx q[3];
rz(-2.8113139) q[3];
sx q[3];
rz(2.3269261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6430214) q[0];
sx q[0];
rz(-1.850147) q[0];
sx q[0];
rz(0.92798573) q[0];
rz(-3.0827177) q[1];
sx q[1];
rz(-1.4480271) q[1];
sx q[1];
rz(-1.7108797) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3979523) q[0];
sx q[0];
rz(-0.40983202) q[0];
sx q[0];
rz(1.3751956) q[0];
x q[1];
rz(2.7785382) q[2];
sx q[2];
rz(-1.3842693) q[2];
sx q[2];
rz(2.8721953) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2637564) q[1];
sx q[1];
rz(-1.0794414) q[1];
sx q[1];
rz(-0.25313913) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9417504) q[3];
sx q[3];
rz(-1.3247196) q[3];
sx q[3];
rz(-1.4914644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4579939) q[2];
sx q[2];
rz(-1.6432089) q[2];
sx q[2];
rz(-1.9963473) q[2];
rz(0.84732071) q[3];
sx q[3];
rz(-1.3342131) q[3];
sx q[3];
rz(-1.5020717) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97750807) q[0];
sx q[0];
rz(-1.1825528) q[0];
sx q[0];
rz(0.384828) q[0];
rz(-2.341914) q[1];
sx q[1];
rz(-0.5189907) q[1];
sx q[1];
rz(-1.5934561) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.528397) q[0];
sx q[0];
rz(-1.3269182) q[0];
sx q[0];
rz(-2.8923678) q[0];
x q[1];
rz(1.9673011) q[2];
sx q[2];
rz(-1.7616211) q[2];
sx q[2];
rz(-0.79939524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1205129) q[1];
sx q[1];
rz(-1.2215633) q[1];
sx q[1];
rz(-1.9053188) q[1];
x q[2];
rz(2.1372651) q[3];
sx q[3];
rz(-2.5133555) q[3];
sx q[3];
rz(0.35997501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3758214) q[2];
sx q[2];
rz(-2.4479726) q[2];
sx q[2];
rz(3.0613464) q[2];
rz(-2.5806228) q[3];
sx q[3];
rz(-1.6601945) q[3];
sx q[3];
rz(0.62275732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65492594) q[0];
sx q[0];
rz(-0.66513649) q[0];
sx q[0];
rz(1.2605865) q[0];
rz(0.27443019) q[1];
sx q[1];
rz(-1.6703037) q[1];
sx q[1];
rz(-0.71896499) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089217606) q[0];
sx q[0];
rz(-1.9196654) q[0];
sx q[0];
rz(2.265076) q[0];
x q[1];
rz(-1.5601129) q[2];
sx q[2];
rz(-2.3950393) q[2];
sx q[2];
rz(-0.03554666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.58108854) q[1];
sx q[1];
rz(-1.6916654) q[1];
sx q[1];
rz(-0.64847364) q[1];
rz(-pi) q[2];
x q[2];
rz(0.218504) q[3];
sx q[3];
rz(-1.9338622) q[3];
sx q[3];
rz(-0.6930815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4869953) q[2];
sx q[2];
rz(-1.8751112) q[2];
sx q[2];
rz(-0.60301644) q[2];
rz(1.6675789) q[3];
sx q[3];
rz(-2.1173729) q[3];
sx q[3];
rz(0.41072861) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5176158) q[0];
sx q[0];
rz(-1.6083953) q[0];
sx q[0];
rz(-2.9543167) q[0];
rz(0.97224832) q[1];
sx q[1];
rz(-2.9863803) q[1];
sx q[1];
rz(0.42047277) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1218106) q[0];
sx q[0];
rz(-0.81392899) q[0];
sx q[0];
rz(-2.5518083) q[0];
x q[1];
rz(-2.114061) q[2];
sx q[2];
rz(-1.2422556) q[2];
sx q[2];
rz(1.3628886) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.17294614) q[1];
sx q[1];
rz(-1.9526281) q[1];
sx q[1];
rz(0.29995014) q[1];
rz(-0.99655788) q[3];
sx q[3];
rz(-1.7645932) q[3];
sx q[3];
rz(-2.7791948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.04574) q[2];
sx q[2];
rz(-1.1008215) q[2];
sx q[2];
rz(0.25071684) q[2];
rz(-1.2480674) q[3];
sx q[3];
rz(-1.7496795) q[3];
sx q[3];
rz(-0.59984508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8641758) q[0];
sx q[0];
rz(-0.59159788) q[0];
sx q[0];
rz(-2.199882) q[0];
rz(3.1031389) q[1];
sx q[1];
rz(-1.7105303) q[1];
sx q[1];
rz(0.11631913) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4424902) q[0];
sx q[0];
rz(-1.1277871) q[0];
sx q[0];
rz(2.3194314) q[0];
rz(-1.3384678) q[2];
sx q[2];
rz(-2.4470377) q[2];
sx q[2];
rz(-1.6603927) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47919905) q[1];
sx q[1];
rz(-2.8251007) q[1];
sx q[1];
rz(2.5280158) q[1];
x q[2];
rz(0.26340254) q[3];
sx q[3];
rz(-0.40168328) q[3];
sx q[3];
rz(0.64611891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8470799) q[2];
sx q[2];
rz(-0.81081644) q[2];
sx q[2];
rz(-1.026356) q[2];
rz(1.2096842) q[3];
sx q[3];
rz(-2.1116202) q[3];
sx q[3];
rz(2.1799555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66861361) q[0];
sx q[0];
rz(-2.0675779) q[0];
sx q[0];
rz(-2.7630254) q[0];
rz(2.0319132) q[1];
sx q[1];
rz(-2.8329284) q[1];
sx q[1];
rz(-3.1139156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33998128) q[0];
sx q[0];
rz(-1.9343573) q[0];
sx q[0];
rz(-2.5021863) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86556025) q[2];
sx q[2];
rz(-2.8606374) q[2];
sx q[2];
rz(1.0298426) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82904774) q[1];
sx q[1];
rz(-1.9282189) q[1];
sx q[1];
rz(0.19695671) q[1];
rz(-pi) q[2];
rz(1.1430986) q[3];
sx q[3];
rz(-1.2790247) q[3];
sx q[3];
rz(1.1455208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3045584) q[2];
sx q[2];
rz(-2.1874032) q[2];
sx q[2];
rz(2.7216116) q[2];
rz(-2.3000681) q[3];
sx q[3];
rz(-1.0171112) q[3];
sx q[3];
rz(-2.7628472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3751752) q[0];
sx q[0];
rz(-3.0247122) q[0];
sx q[0];
rz(-2.8564659) q[0];
rz(-2.9410703) q[1];
sx q[1];
rz(-1.9384117) q[1];
sx q[1];
rz(-0.83555046) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886166) q[0];
sx q[0];
rz(-1.1167913) q[0];
sx q[0];
rz(0.80843057) q[0];
rz(1.5642688) q[2];
sx q[2];
rz(-2.2178429) q[2];
sx q[2];
rz(-1.4359635) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.056014765) q[1];
sx q[1];
rz(-1.2390572) q[1];
sx q[1];
rz(-0.47420926) q[1];
rz(-pi) q[2];
rz(1.6099168) q[3];
sx q[3];
rz(-1.3353383) q[3];
sx q[3];
rz(2.2155516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7150813) q[2];
sx q[2];
rz(-2.0501523) q[2];
sx q[2];
rz(0.69768989) q[2];
rz(-0.52608025) q[3];
sx q[3];
rz(-2.3487921) q[3];
sx q[3];
rz(-1.6349767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0195011) q[0];
sx q[0];
rz(-0.87845907) q[0];
sx q[0];
rz(-0.39977951) q[0];
rz(-1.1493692) q[1];
sx q[1];
rz(-0.72723564) q[1];
sx q[1];
rz(-0.74692187) q[1];
rz(0.34869817) q[2];
sx q[2];
rz(-0.92338466) q[2];
sx q[2];
rz(0.79104214) q[2];
rz(-2.6769842) q[3];
sx q[3];
rz(-2.5425442) q[3];
sx q[3];
rz(3.0221593) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
