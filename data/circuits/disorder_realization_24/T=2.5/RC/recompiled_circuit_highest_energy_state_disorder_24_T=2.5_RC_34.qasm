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
rz(-0.51198045) q[0];
sx q[0];
rz(-2.5706302) q[0];
sx q[0];
rz(-2.9538739) q[0];
rz(0.81387782) q[1];
sx q[1];
rz(4.6133572) q[1];
sx q[1];
rz(7.5724966) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1027689) q[0];
sx q[0];
rz(-1.3631522) q[0];
sx q[0];
rz(2.8994292) q[0];
rz(1.600483) q[2];
sx q[2];
rz(-2.1412537) q[2];
sx q[2];
rz(2.1726959) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.52367585) q[1];
sx q[1];
rz(-0.53046983) q[1];
sx q[1];
rz(-0.50650774) q[1];
rz(-pi) q[2];
rz(1.2175519) q[3];
sx q[3];
rz(-1.3526256) q[3];
sx q[3];
rz(0.45807236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8744371) q[2];
sx q[2];
rz(-1.1031373) q[2];
sx q[2];
rz(-2.3647986) q[2];
rz(-1.8393501) q[3];
sx q[3];
rz(-1.6603419) q[3];
sx q[3];
rz(1.9384025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5341107) q[0];
sx q[0];
rz(-2.6597839) q[0];
sx q[0];
rz(2.1433461) q[0];
rz(-0.36969319) q[1];
sx q[1];
rz(-1.0414711) q[1];
sx q[1];
rz(1.1179771) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3415754) q[0];
sx q[0];
rz(-1.9722) q[0];
sx q[0];
rz(-0.97824162) q[0];
rz(0.68685617) q[2];
sx q[2];
rz(-2.0389028) q[2];
sx q[2];
rz(-1.3324225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4395968) q[1];
sx q[1];
rz(-2.0148104) q[1];
sx q[1];
rz(-0.16298144) q[1];
rz(1.9309773) q[3];
sx q[3];
rz(-0.66347117) q[3];
sx q[3];
rz(2.0455895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.86576858) q[2];
sx q[2];
rz(-1.3752702) q[2];
sx q[2];
rz(-2.7562874) q[2];
rz(3.1028808) q[3];
sx q[3];
rz(-0.35274115) q[3];
sx q[3];
rz(-2.501131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5019219) q[0];
sx q[0];
rz(-0.47698912) q[0];
sx q[0];
rz(1.3013526) q[0];
rz(-3.0838857) q[1];
sx q[1];
rz(-2.6951908) q[1];
sx q[1];
rz(1.8147963) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2952174) q[0];
sx q[0];
rz(-1.0344939) q[0];
sx q[0];
rz(0.55732507) q[0];
x q[1];
rz(-0.71452629) q[2];
sx q[2];
rz(-2.9305612) q[2];
sx q[2];
rz(0.85725609) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51242764) q[1];
sx q[1];
rz(-1.2841793) q[1];
sx q[1];
rz(2.5932339) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7884198) q[3];
sx q[3];
rz(-0.88017094) q[3];
sx q[3];
rz(1.4598802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79746276) q[2];
sx q[2];
rz(-0.8351438) q[2];
sx q[2];
rz(-0.75378913) q[2];
rz(-1.3695184) q[3];
sx q[3];
rz(-1.6794208) q[3];
sx q[3];
rz(-2.4863825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21514431) q[0];
sx q[0];
rz(-1.0855874) q[0];
sx q[0];
rz(-1.3947067) q[0];
rz(-1.9649547) q[1];
sx q[1];
rz(-1.5086915) q[1];
sx q[1];
rz(-0.36697695) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6539388) q[0];
sx q[0];
rz(-1.6375082) q[0];
sx q[0];
rz(0.056746527) q[0];
x q[1];
rz(-1.3064477) q[2];
sx q[2];
rz(-1.8943139) q[2];
sx q[2];
rz(2.8217962) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6527892) q[1];
sx q[1];
rz(-1.7646043) q[1];
sx q[1];
rz(0.36313063) q[1];
rz(0.36716299) q[3];
sx q[3];
rz(-1.490574) q[3];
sx q[3];
rz(-2.1047161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41161141) q[2];
sx q[2];
rz(-1.2282164) q[2];
sx q[2];
rz(-1.270594) q[2];
rz(0.96108428) q[3];
sx q[3];
rz(-1.4189439) q[3];
sx q[3];
rz(0.050475033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56080317) q[0];
sx q[0];
rz(-0.92629782) q[0];
sx q[0];
rz(-1.3128989) q[0];
rz(-1.987223) q[1];
sx q[1];
rz(-1.9998974) q[1];
sx q[1];
rz(-0.55783522) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2656888) q[0];
sx q[0];
rz(-0.77123986) q[0];
sx q[0];
rz(-0.79637595) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9809772) q[2];
sx q[2];
rz(-1.4945507) q[2];
sx q[2];
rz(2.3996224) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.055649672) q[1];
sx q[1];
rz(-0.80971566) q[1];
sx q[1];
rz(1.147406) q[1];
rz(-0.29739012) q[3];
sx q[3];
rz(-0.47978401) q[3];
sx q[3];
rz(-0.23272091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.64342) q[2];
sx q[2];
rz(-1.8489445) q[2];
sx q[2];
rz(2.2922929) q[2];
rz(0.15527209) q[3];
sx q[3];
rz(-2.1657491) q[3];
sx q[3];
rz(-0.37080216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7840541) q[0];
sx q[0];
rz(-2.5548866) q[0];
sx q[0];
rz(-0.46911711) q[0];
rz(2.6672089) q[1];
sx q[1];
rz(-0.7862888) q[1];
sx q[1];
rz(-0.75538409) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2952399) q[0];
sx q[0];
rz(-1.8420514) q[0];
sx q[0];
rz(-3.1370107) q[0];
rz(-0.62894507) q[2];
sx q[2];
rz(-1.2997932) q[2];
sx q[2];
rz(-0.48297255) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0363706) q[1];
sx q[1];
rz(-1.7835938) q[1];
sx q[1];
rz(-0.12963055) q[1];
rz(2.7128025) q[3];
sx q[3];
rz(-2.2615823) q[3];
sx q[3];
rz(-2.5429728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0581806) q[2];
sx q[2];
rz(-1.3231134) q[2];
sx q[2];
rz(-2.9537436) q[2];
rz(-1.2457054) q[3];
sx q[3];
rz(-1.3789504) q[3];
sx q[3];
rz(1.0904306) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.93194) q[0];
sx q[0];
rz(-1.2670452) q[0];
sx q[0];
rz(-2.7506822) q[0];
rz(0.93005013) q[1];
sx q[1];
rz(-1.7101733) q[1];
sx q[1];
rz(-1.3040868) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6041038) q[0];
sx q[0];
rz(-1.2199645) q[0];
sx q[0];
rz(-3.119191) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9043998) q[2];
sx q[2];
rz(-0.35065864) q[2];
sx q[2];
rz(1.2192977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2432989) q[1];
sx q[1];
rz(-2.0212681) q[1];
sx q[1];
rz(2.9551278) q[1];
x q[2];
rz(-1.5897069) q[3];
sx q[3];
rz(-1.2745556) q[3];
sx q[3];
rz(2.4774574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2685252) q[2];
sx q[2];
rz(-2.8748942) q[2];
sx q[2];
rz(0.11338691) q[2];
rz(-0.015965613) q[3];
sx q[3];
rz(-1.6458052) q[3];
sx q[3];
rz(0.16460831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78405821) q[0];
sx q[0];
rz(-1.6679732) q[0];
sx q[0];
rz(-3.0175324) q[0];
rz(-0.1217753) q[1];
sx q[1];
rz(-2.3901794) q[1];
sx q[1];
rz(1.6990936) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7989726) q[0];
sx q[0];
rz(-2.7443287) q[0];
sx q[0];
rz(-1.8333866) q[0];
rz(-pi) q[1];
rz(-0.29192544) q[2];
sx q[2];
rz(-1.6124083) q[2];
sx q[2];
rz(0.49753639) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0066083) q[1];
sx q[1];
rz(-1.4378752) q[1];
sx q[1];
rz(-1.9167109) q[1];
rz(-pi) q[2];
rz(0.27713953) q[3];
sx q[3];
rz(-0.74944121) q[3];
sx q[3];
rz(-0.14271862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.96559912) q[2];
sx q[2];
rz(-1.3600574) q[2];
sx q[2];
rz(2.3967801) q[2];
rz(-0.92721573) q[3];
sx q[3];
rz(-1.6330481) q[3];
sx q[3];
rz(2.0492699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6532779) q[0];
sx q[0];
rz(-1.7019685) q[0];
sx q[0];
rz(-2.9162245) q[0];
rz(-1.1124181) q[1];
sx q[1];
rz(-2.367159) q[1];
sx q[1];
rz(-2.2427799) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32466896) q[0];
sx q[0];
rz(-2.3153439) q[0];
sx q[0];
rz(-0.81314317) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5139047) q[2];
sx q[2];
rz(-1.9348839) q[2];
sx q[2];
rz(0.68825005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0045709) q[1];
sx q[1];
rz(-0.70024509) q[1];
sx q[1];
rz(1.2546468) q[1];
rz(-pi) q[2];
rz(1.2779253) q[3];
sx q[3];
rz(-1.5808269) q[3];
sx q[3];
rz(2.9603279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.044346873) q[2];
sx q[2];
rz(-1.4298507) q[2];
sx q[2];
rz(2.782235) q[2];
rz(2.8906631) q[3];
sx q[3];
rz(-1.072262) q[3];
sx q[3];
rz(0.76771626) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74442416) q[0];
sx q[0];
rz(-0.1305307) q[0];
sx q[0];
rz(-0.8771483) q[0];
rz(-1.5769222) q[1];
sx q[1];
rz(-1.501187) q[1];
sx q[1];
rz(-2.3559779) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3451259) q[0];
sx q[0];
rz(-0.32497999) q[0];
sx q[0];
rz(2.3152405) q[0];
rz(-1.9989817) q[2];
sx q[2];
rz(-2.0642363) q[2];
sx q[2];
rz(-2.7827415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0141416) q[1];
sx q[1];
rz(-1.2268865) q[1];
sx q[1];
rz(1.6456855) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7212058) q[3];
sx q[3];
rz(-2.0724943) q[3];
sx q[3];
rz(-3.0034163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2239573) q[2];
sx q[2];
rz(-0.79006299) q[2];
sx q[2];
rz(0.7381953) q[2];
rz(-0.92292845) q[3];
sx q[3];
rz(-1.3083369) q[3];
sx q[3];
rz(2.8452528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7947163) q[0];
sx q[0];
rz(-1.7807757) q[0];
sx q[0];
rz(-2.1697252) q[0];
rz(-0.90732668) q[1];
sx q[1];
rz(-2.0569888) q[1];
sx q[1];
rz(2.8203698) q[1];
rz(-0.99030607) q[2];
sx q[2];
rz(-1.2377501) q[2];
sx q[2];
rz(-0.35828423) q[2];
rz(-0.76413705) q[3];
sx q[3];
rz(-2.4485179) q[3];
sx q[3];
rz(-2.5701523) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
