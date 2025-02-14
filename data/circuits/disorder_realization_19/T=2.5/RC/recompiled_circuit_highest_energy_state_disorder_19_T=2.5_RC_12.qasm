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
rz(2.1228696) q[0];
sx q[0];
rz(-2.2824204) q[0];
sx q[0];
rz(-2.3265042) q[0];
rz(1.9563142) q[1];
sx q[1];
rz(-1.7307245) q[1];
sx q[1];
rz(-1.0676395) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1045221) q[0];
sx q[0];
rz(-0.43209729) q[0];
sx q[0];
rz(-1.62754) q[0];
rz(-pi) q[1];
rz(-0.63500603) q[2];
sx q[2];
rz(-2.5189812) q[2];
sx q[2];
rz(2.9787763) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8274535) q[1];
sx q[1];
rz(-1.8874536) q[1];
sx q[1];
rz(-3.0400671) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26431636) q[3];
sx q[3];
rz(-2.2393763) q[3];
sx q[3];
rz(0.17449915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4472569) q[2];
sx q[2];
rz(-0.7434291) q[2];
sx q[2];
rz(1.2747964) q[2];
rz(0.46191195) q[3];
sx q[3];
rz(-2.4670944) q[3];
sx q[3];
rz(-1.1326724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3502515) q[0];
sx q[0];
rz(-2.8332062) q[0];
sx q[0];
rz(-1.8885008) q[0];
rz(-0.14532267) q[1];
sx q[1];
rz(-1.3959613) q[1];
sx q[1];
rz(-1.0911509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017830124) q[0];
sx q[0];
rz(-1.1457256) q[0];
sx q[0];
rz(0.2353038) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45812313) q[2];
sx q[2];
rz(-1.2451742) q[2];
sx q[2];
rz(1.2026915) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1205463) q[1];
sx q[1];
rz(-1.7851549) q[1];
sx q[1];
rz(1.3938741) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1717779) q[3];
sx q[3];
rz(-1.1091091) q[3];
sx q[3];
rz(-0.53129133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48065177) q[2];
sx q[2];
rz(-1.7512243) q[2];
sx q[2];
rz(2.6118028) q[2];
rz(-0.79331136) q[3];
sx q[3];
rz(-1.6042234) q[3];
sx q[3];
rz(-0.62354273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48524258) q[0];
sx q[0];
rz(-2.1915477) q[0];
sx q[0];
rz(-2.6128838) q[0];
rz(-0.54620019) q[1];
sx q[1];
rz(-0.9587973) q[1];
sx q[1];
rz(0.34034696) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8204644) q[0];
sx q[0];
rz(-1.8978776) q[0];
sx q[0];
rz(1.8595247) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5916948) q[2];
sx q[2];
rz(-0.27355121) q[2];
sx q[2];
rz(3.105148) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18225154) q[1];
sx q[1];
rz(-1.6653582) q[1];
sx q[1];
rz(1.574081) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3119218) q[3];
sx q[3];
rz(-1.2874914) q[3];
sx q[3];
rz(-1.0782575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9693552) q[2];
sx q[2];
rz(-0.41219553) q[2];
sx q[2];
rz(2.7117512) q[2];
rz(2.3853081) q[3];
sx q[3];
rz(-0.1736621) q[3];
sx q[3];
rz(2.3156796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38465685) q[0];
sx q[0];
rz(-1.5058368) q[0];
sx q[0];
rz(1.4935619) q[0];
rz(2.3736296) q[1];
sx q[1];
rz(-2.6710644) q[1];
sx q[1];
rz(2.5951662) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35432974) q[0];
sx q[0];
rz(-3.0453735) q[0];
sx q[0];
rz(2.3090099) q[0];
rz(2.619952) q[2];
sx q[2];
rz(-1.120479) q[2];
sx q[2];
rz(-0.7618103) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.029303251) q[1];
sx q[1];
rz(-1.7186856) q[1];
sx q[1];
rz(1.6690955) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68164556) q[3];
sx q[3];
rz(-2.1053616) q[3];
sx q[3];
rz(-0.93972423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4297318) q[2];
sx q[2];
rz(-1.6542566) q[2];
sx q[2];
rz(2.8374953) q[2];
rz(-1.7763304) q[3];
sx q[3];
rz(-1.920776) q[3];
sx q[3];
rz(2.064866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.6292608) q[0];
sx q[0];
rz(-0.4442232) q[0];
sx q[0];
rz(-0.00057922676) q[0];
rz(-2.0594788) q[1];
sx q[1];
rz(-0.52434701) q[1];
sx q[1];
rz(0.93926114) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9094641) q[0];
sx q[0];
rz(-0.8117903) q[0];
sx q[0];
rz(-0.62007298) q[0];
rz(-pi) q[1];
rz(1.0910735) q[2];
sx q[2];
rz(-2.4726598) q[2];
sx q[2];
rz(2.875653) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1402005) q[1];
sx q[1];
rz(-1.0085229) q[1];
sx q[1];
rz(-2.8881025) q[1];
rz(-pi) q[2];
rz(2.9143798) q[3];
sx q[3];
rz(-1.876653) q[3];
sx q[3];
rz(1.3671041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.038736343) q[2];
sx q[2];
rz(-1.2612217) q[2];
sx q[2];
rz(2.8803414) q[2];
rz(2.2767565) q[3];
sx q[3];
rz(-1.7295001) q[3];
sx q[3];
rz(-0.29204667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84457266) q[0];
sx q[0];
rz(-1.7140056) q[0];
sx q[0];
rz(0.20508668) q[0];
rz(1.0467485) q[1];
sx q[1];
rz(-1.3958684) q[1];
sx q[1];
rz(-2.7632025) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8240067) q[0];
sx q[0];
rz(-1.7574302) q[0];
sx q[0];
rz(-1.179479) q[0];
x q[1];
rz(-1.7360052) q[2];
sx q[2];
rz(-1.1045611) q[2];
sx q[2];
rz(1.2638484) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4503895) q[1];
sx q[1];
rz(-0.65237633) q[1];
sx q[1];
rz(2.1689586) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50726733) q[3];
sx q[3];
rz(-0.75921042) q[3];
sx q[3];
rz(1.678148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4947074) q[2];
sx q[2];
rz(-1.9834221) q[2];
sx q[2];
rz(1.0542487) q[2];
rz(-2.4833637) q[3];
sx q[3];
rz(-0.64526486) q[3];
sx q[3];
rz(0.48404199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9996416) q[0];
sx q[0];
rz(-1.9331837) q[0];
sx q[0];
rz(-2.5010338) q[0];
rz(-0.41796747) q[1];
sx q[1];
rz(-1.9211946) q[1];
sx q[1];
rz(-2.3366065) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1520839) q[0];
sx q[0];
rz(-1.594127) q[0];
sx q[0];
rz(-2.4341466) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4067602) q[2];
sx q[2];
rz(-0.8443588) q[2];
sx q[2];
rz(-0.11876362) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1678214) q[1];
sx q[1];
rz(-0.034618363) q[1];
sx q[1];
rz(0.95914118) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42607362) q[3];
sx q[3];
rz(-2.5059627) q[3];
sx q[3];
rz(2.9866708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5927222) q[2];
sx q[2];
rz(-2.1465116) q[2];
sx q[2];
rz(0.76118809) q[2];
rz(-0.74782863) q[3];
sx q[3];
rz(-1.3176094) q[3];
sx q[3];
rz(-2.4650011) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51650301) q[0];
sx q[0];
rz(-2.857132) q[0];
sx q[0];
rz(-1.4768584) q[0];
rz(-2.6853216) q[1];
sx q[1];
rz(-1.7218593) q[1];
sx q[1];
rz(-0.87108535) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3696156) q[0];
sx q[0];
rz(-2.5273558) q[0];
sx q[0];
rz(-3.1200039) q[0];
rz(0.34452166) q[2];
sx q[2];
rz(-1.7167175) q[2];
sx q[2];
rz(-0.61071705) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2312517) q[1];
sx q[1];
rz(-1.7897072) q[1];
sx q[1];
rz(-0.33556767) q[1];
rz(3.0814287) q[3];
sx q[3];
rz(-1.9964661) q[3];
sx q[3];
rz(-0.7983467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2386834) q[2];
sx q[2];
rz(-2.6148655) q[2];
sx q[2];
rz(0.61863679) q[2];
rz(3.1213308) q[3];
sx q[3];
rz(-2.198115) q[3];
sx q[3];
rz(2.3616135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0778377) q[0];
sx q[0];
rz(-1.5564593) q[0];
sx q[0];
rz(-0.42386398) q[0];
rz(0.18691143) q[1];
sx q[1];
rz(-2.321545) q[1];
sx q[1];
rz(-1.4580457) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33015051) q[0];
sx q[0];
rz(-1.6887661) q[0];
sx q[0];
rz(-1.5565722) q[0];
rz(-pi) q[1];
rz(-2.0293268) q[2];
sx q[2];
rz(-2.5800309) q[2];
sx q[2];
rz(2.4379345) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.600454) q[1];
sx q[1];
rz(-0.14800528) q[1];
sx q[1];
rz(2.7861094) q[1];
rz(-pi) q[2];
rz(1.4814754) q[3];
sx q[3];
rz(-1.9253988) q[3];
sx q[3];
rz(-1.8338628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75108782) q[2];
sx q[2];
rz(-2.0743399) q[2];
sx q[2];
rz(-1.0941774) q[2];
rz(-1.4607726) q[3];
sx q[3];
rz(-1.7655617) q[3];
sx q[3];
rz(1.8028397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3222892) q[0];
sx q[0];
rz(-1.7604473) q[0];
sx q[0];
rz(-1.5018916) q[0];
rz(-2.8627401) q[1];
sx q[1];
rz(-2.1809705) q[1];
sx q[1];
rz(1.0241114) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8968643) q[0];
sx q[0];
rz(-2.6435268) q[0];
sx q[0];
rz(2.9480272) q[0];
x q[1];
rz(0.050692888) q[2];
sx q[2];
rz(-1.2963352) q[2];
sx q[2];
rz(1.7982184) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2453277) q[1];
sx q[1];
rz(-0.30054856) q[1];
sx q[1];
rz(-0.73722571) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34124438) q[3];
sx q[3];
rz(-0.8559209) q[3];
sx q[3];
rz(1.9174066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9378822) q[2];
sx q[2];
rz(-1.8546162) q[2];
sx q[2];
rz(-0.53696519) q[2];
rz(1.2282061) q[3];
sx q[3];
rz(-2.1824586) q[3];
sx q[3];
rz(2.8502407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1163597) q[0];
sx q[0];
rz(-1.7740842) q[0];
sx q[0];
rz(-2.0728003) q[0];
rz(-2.4346726) q[1];
sx q[1];
rz(-2.0126577) q[1];
sx q[1];
rz(-0.001002034) q[1];
rz(1.3798643) q[2];
sx q[2];
rz(-1.9700865) q[2];
sx q[2];
rz(-1.3843591) q[2];
rz(-0.14866004) q[3];
sx q[3];
rz(-1.3547263) q[3];
sx q[3];
rz(2.3480036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
