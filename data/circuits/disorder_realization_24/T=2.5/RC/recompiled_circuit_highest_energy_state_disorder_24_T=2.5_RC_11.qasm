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
rz(-1.6698281) q[1];
sx q[1];
rz(-1.8522813) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1027689) q[0];
sx q[0];
rz(-1.3631522) q[0];
sx q[0];
rz(2.8994292) q[0];
x q[1];
rz(1.600483) q[2];
sx q[2];
rz(-1.0003389) q[2];
sx q[2];
rz(0.96889673) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.047886176) q[1];
sx q[1];
rz(-1.1125065) q[1];
sx q[1];
rz(-1.2935733) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23203316) q[3];
sx q[3];
rz(-1.9153144) q[3];
sx q[3];
rz(-1.0330878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2671555) q[2];
sx q[2];
rz(-2.0384553) q[2];
sx q[2];
rz(-0.77679408) q[2];
rz(-1.8393501) q[3];
sx q[3];
rz(-1.6603419) q[3];
sx q[3];
rz(1.9384025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5341107) q[0];
sx q[0];
rz(-0.48180875) q[0];
sx q[0];
rz(0.99824655) q[0];
rz(-0.36969319) q[1];
sx q[1];
rz(-1.0414711) q[1];
sx q[1];
rz(-2.0236156) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29660615) q[0];
sx q[0];
rz(-0.70195508) q[0];
sx q[0];
rz(-0.92092307) q[0];
rz(-2.4684899) q[2];
sx q[2];
rz(-0.80922283) q[2];
sx q[2];
rz(0.26462091) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0740023) q[1];
sx q[1];
rz(-2.6704881) q[1];
sx q[1];
rz(1.2420688) q[1];
rz(-pi) q[2];
rz(2.2023674) q[3];
sx q[3];
rz(-1.3520006) q[3];
sx q[3];
rz(-2.95524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86576858) q[2];
sx q[2];
rz(-1.7663225) q[2];
sx q[2];
rz(0.38530525) q[2];
rz(3.1028808) q[3];
sx q[3];
rz(-2.7888515) q[3];
sx q[3];
rz(2.501131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63967079) q[0];
sx q[0];
rz(-2.6646035) q[0];
sx q[0];
rz(1.84024) q[0];
rz(0.057706984) q[1];
sx q[1];
rz(-0.4464018) q[1];
sx q[1];
rz(-1.8147963) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1039057) q[0];
sx q[0];
rz(-0.7531727) q[0];
sx q[0];
rz(0.84363787) q[0];
rz(0.71452629) q[2];
sx q[2];
rz(-2.9305612) q[2];
sx q[2];
rz(-0.85725609) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2542299) q[1];
sx q[1];
rz(-1.0471724) q[1];
sx q[1];
rz(1.2382522) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1744099) q[3];
sx q[3];
rz(-0.76226888) q[3];
sx q[3];
rz(-1.1571194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3441299) q[2];
sx q[2];
rz(-0.8351438) q[2];
sx q[2];
rz(0.75378913) q[2];
rz(1.7720743) q[3];
sx q[3];
rz(-1.4621719) q[3];
sx q[3];
rz(-0.65521017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21514431) q[0];
sx q[0];
rz(-2.0560052) q[0];
sx q[0];
rz(-1.3947067) q[0];
rz(1.1766379) q[1];
sx q[1];
rz(-1.5086915) q[1];
sx q[1];
rz(2.7746157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1934051) q[0];
sx q[0];
rz(-0.087554878) q[0];
sx q[0];
rz(2.2746536) q[0];
rz(-pi) q[1];
rz(-2.4796333) q[2];
sx q[2];
rz(-0.41482224) q[2];
sx q[2];
rz(0.38554672) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1550459) q[1];
sx q[1];
rz(-1.2147702) q[1];
sx q[1];
rz(1.3638391) q[1];
x q[2];
rz(-0.36716299) q[3];
sx q[3];
rz(-1.6510186) q[3];
sx q[3];
rz(1.0368766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.41161141) q[2];
sx q[2];
rz(-1.9133762) q[2];
sx q[2];
rz(-1.270594) q[2];
rz(0.96108428) q[3];
sx q[3];
rz(-1.4189439) q[3];
sx q[3];
rz(-3.0911176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56080317) q[0];
sx q[0];
rz(-0.92629782) q[0];
sx q[0];
rz(1.3128989) q[0];
rz(-1.1543697) q[1];
sx q[1];
rz(-1.1416953) q[1];
sx q[1];
rz(2.5837574) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2656888) q[0];
sx q[0];
rz(-0.77123986) q[0];
sx q[0];
rz(-2.3452167) q[0];
rz(-pi) q[1];
rz(1.4935605) q[2];
sx q[2];
rz(-1.7309411) q[2];
sx q[2];
rz(0.84116615) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8164929) q[1];
sx q[1];
rz(-1.268728) q[1];
sx q[1];
rz(-0.80728553) q[1];
rz(-0.29739012) q[3];
sx q[3];
rz(-0.47978401) q[3];
sx q[3];
rz(-0.23272091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4981726) q[2];
sx q[2];
rz(-1.8489445) q[2];
sx q[2];
rz(-0.84929973) q[2];
rz(-2.9863206) q[3];
sx q[3];
rz(-2.1657491) q[3];
sx q[3];
rz(-0.37080216) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3575386) q[0];
sx q[0];
rz(-0.58670601) q[0];
sx q[0];
rz(-0.46911711) q[0];
rz(2.6672089) q[1];
sx q[1];
rz(-2.3553039) q[1];
sx q[1];
rz(-2.3862086) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82925382) q[0];
sx q[0];
rz(-2.8702998) q[0];
sx q[0];
rz(1.5872699) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5126476) q[2];
sx q[2];
rz(-1.8417995) q[2];
sx q[2];
rz(2.6586201) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6346917) q[1];
sx q[1];
rz(-1.4441057) q[1];
sx q[1];
rz(1.3562528) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3085318) q[3];
sx q[3];
rz(-1.2445881) q[3];
sx q[3];
rz(1.8859832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0834121) q[2];
sx q[2];
rz(-1.3231134) q[2];
sx q[2];
rz(2.9537436) q[2];
rz(1.8958873) q[3];
sx q[3];
rz(-1.3789504) q[3];
sx q[3];
rz(1.0904306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2096527) q[0];
sx q[0];
rz(-1.2670452) q[0];
sx q[0];
rz(2.7506822) q[0];
rz(-2.2115425) q[1];
sx q[1];
rz(-1.7101733) q[1];
sx q[1];
rz(1.8375058) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.539003) q[0];
sx q[0];
rz(-2.7900759) q[0];
sx q[0];
rz(-1.6319265) q[0];
x q[1];
rz(-1.6565336) q[2];
sx q[2];
rz(-1.2303599) q[2];
sx q[2];
rz(-2.1742333) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3871349) q[1];
sx q[1];
rz(-1.4031193) q[1];
sx q[1];
rz(-2.0281591) q[1];
rz(-pi) q[2];
x q[2];
rz(2.845302) q[3];
sx q[3];
rz(-1.5888831) q[3];
sx q[3];
rz(0.91218218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8730674) q[2];
sx q[2];
rz(-0.26669845) q[2];
sx q[2];
rz(-0.11338691) q[2];
rz(-0.015965613) q[3];
sx q[3];
rz(-1.6458052) q[3];
sx q[3];
rz(-2.9769843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78405821) q[0];
sx q[0];
rz(-1.6679732) q[0];
sx q[0];
rz(0.12406021) q[0];
rz(0.1217753) q[1];
sx q[1];
rz(-2.3901794) q[1];
sx q[1];
rz(1.442499) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0589813) q[0];
sx q[0];
rz(-1.9537203) q[0];
sx q[0];
rz(-0.10848606) q[0];
rz(-pi) q[1];
rz(1.5273483) q[2];
sx q[2];
rz(-1.8624616) q[2];
sx q[2];
rz(1.0857605) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1349843) q[1];
sx q[1];
rz(-1.4378752) q[1];
sx q[1];
rz(1.2248817) q[1];
rz(-pi) q[2];
rz(2.4115218) q[3];
sx q[3];
rz(-1.3833117) q[3];
sx q[3];
rz(1.9188251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96559912) q[2];
sx q[2];
rz(-1.3600574) q[2];
sx q[2];
rz(-0.74481258) q[2];
rz(-2.2143769) q[3];
sx q[3];
rz(-1.5085446) q[3];
sx q[3];
rz(-1.0923227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48831478) q[0];
sx q[0];
rz(-1.4396242) q[0];
sx q[0];
rz(-2.9162245) q[0];
rz(1.1124181) q[1];
sx q[1];
rz(-0.77443361) q[1];
sx q[1];
rz(-2.2427799) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5170636) q[0];
sx q[0];
rz(-1.0072021) q[0];
sx q[0];
rz(0.64080142) q[0];
x q[1];
rz(2.4821539) q[2];
sx q[2];
rz(-0.62022479) q[2];
sx q[2];
rz(0.31977113) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.73279335) q[1];
sx q[1];
rz(-0.91161722) q[1];
sx q[1];
rz(-0.2562457) q[1];
rz(-1.2779253) q[3];
sx q[3];
rz(-1.5607657) q[3];
sx q[3];
rz(-0.18126479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.044346873) q[2];
sx q[2];
rz(-1.4298507) q[2];
sx q[2];
rz(-2.782235) q[2];
rz(2.8906631) q[3];
sx q[3];
rz(-2.0693306) q[3];
sx q[3];
rz(2.3738764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.74442416) q[0];
sx q[0];
rz(-3.011062) q[0];
sx q[0];
rz(0.8771483) q[0];
rz(-1.5646704) q[1];
sx q[1];
rz(-1.501187) q[1];
sx q[1];
rz(-0.78561479) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1164074) q[0];
sx q[0];
rz(-1.8078363) q[0];
sx q[0];
rz(-2.9171506) q[0];
rz(-pi) q[1];
x q[1];
rz(1.142611) q[2];
sx q[2];
rz(-1.0773563) q[2];
sx q[2];
rz(-0.35885119) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4686376) q[1];
sx q[1];
rz(-1.6412927) q[1];
sx q[1];
rz(-0.34480178) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7212058) q[3];
sx q[3];
rz(-2.0724943) q[3];
sx q[3];
rz(0.13817638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2239573) q[2];
sx q[2];
rz(-2.3515297) q[2];
sx q[2];
rz(0.7381953) q[2];
rz(-0.92292845) q[3];
sx q[3];
rz(-1.8332558) q[3];
sx q[3];
rz(0.2963399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.7947163) q[0];
sx q[0];
rz(-1.3608169) q[0];
sx q[0];
rz(0.97186744) q[0];
rz(-2.234266) q[1];
sx q[1];
rz(-1.0846039) q[1];
sx q[1];
rz(-0.32122282) q[1];
rz(0.39225929) q[2];
sx q[2];
rz(-2.1156559) q[2];
sx q[2];
rz(1.0012913) q[2];
rz(-1.0492269) q[3];
sx q[3];
rz(-2.0502301) q[3];
sx q[3];
rz(-0.32296317) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
