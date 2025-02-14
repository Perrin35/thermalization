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
rz(2.6296122) q[0];
sx q[0];
rz(-0.57096243) q[0];
sx q[0];
rz(-0.1877187) q[0];
rz(-2.3277148) q[1];
sx q[1];
rz(-1.4717646) q[1];
sx q[1];
rz(-1.2893113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5587414) q[0];
sx q[0];
rz(-1.8076573) q[0];
sx q[0];
rz(1.357097) q[0];
rz(-pi) q[1];
rz(-1.600483) q[2];
sx q[2];
rz(-2.1412537) q[2];
sx q[2];
rz(-2.1726959) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0937065) q[1];
sx q[1];
rz(-1.1125065) q[1];
sx q[1];
rz(1.8480193) q[1];
x q[2];
rz(-1.0008775) q[3];
sx q[3];
rz(-2.7288247) q[3];
sx q[3];
rz(1.4980396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8744371) q[2];
sx q[2];
rz(-1.1031373) q[2];
sx q[2];
rz(0.77679408) q[2];
rz(-1.8393501) q[3];
sx q[3];
rz(-1.6603419) q[3];
sx q[3];
rz(-1.2031901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60748196) q[0];
sx q[0];
rz(-0.48180875) q[0];
sx q[0];
rz(-2.1433461) q[0];
rz(-0.36969319) q[1];
sx q[1];
rz(-2.1001215) q[1];
sx q[1];
rz(2.0236156) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3415754) q[0];
sx q[0];
rz(-1.9722) q[0];
sx q[0];
rz(-0.97824162) q[0];
rz(-pi) q[1];
rz(2.4547365) q[2];
sx q[2];
rz(-2.0389028) q[2];
sx q[2];
rz(1.3324225) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4395968) q[1];
sx q[1];
rz(-2.0148104) q[1];
sx q[1];
rz(-0.16298144) q[1];
x q[2];
rz(-0.93922521) q[3];
sx q[3];
rz(-1.7895921) q[3];
sx q[3];
rz(-0.18635264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.86576858) q[2];
sx q[2];
rz(-1.3752702) q[2];
sx q[2];
rz(-0.38530525) q[2];
rz(-3.1028808) q[3];
sx q[3];
rz(-0.35274115) q[3];
sx q[3];
rz(2.501131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5019219) q[0];
sx q[0];
rz(-2.6646035) q[0];
sx q[0];
rz(1.84024) q[0];
rz(-3.0838857) q[1];
sx q[1];
rz(-0.4464018) q[1];
sx q[1];
rz(-1.8147963) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1039057) q[0];
sx q[0];
rz(-0.7531727) q[0];
sx q[0];
rz(0.84363787) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9811601) q[2];
sx q[2];
rz(-1.7084873) q[2];
sx q[2];
rz(1.4169803) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2542299) q[1];
sx q[1];
rz(-1.0471724) q[1];
sx q[1];
rz(-1.2382522) q[1];
x q[2];
rz(2.7884198) q[3];
sx q[3];
rz(-0.88017094) q[3];
sx q[3];
rz(1.4598802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79746276) q[2];
sx q[2];
rz(-2.3064488) q[2];
sx q[2];
rz(2.3878035) q[2];
rz(-1.3695184) q[3];
sx q[3];
rz(-1.6794208) q[3];
sx q[3];
rz(-2.4863825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-0.21514431) q[0];
sx q[0];
rz(-2.0560052) q[0];
sx q[0];
rz(1.7468859) q[0];
rz(1.1766379) q[1];
sx q[1];
rz(-1.5086915) q[1];
sx q[1];
rz(-0.36697695) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6539388) q[0];
sx q[0];
rz(-1.6375082) q[0];
sx q[0];
rz(-3.0848461) q[0];
x q[1];
rz(-0.66195935) q[2];
sx q[2];
rz(-2.7267704) q[2];
sx q[2];
rz(0.38554672) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9865468) q[1];
sx q[1];
rz(-1.9268225) q[1];
sx q[1];
rz(-1.7777535) q[1];
rz(-pi) q[2];
rz(-1.4848726) q[3];
sx q[3];
rz(-1.20487) q[3];
sx q[3];
rz(-2.576862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7299812) q[2];
sx q[2];
rz(-1.9133762) q[2];
sx q[2];
rz(-1.270594) q[2];
rz(-0.96108428) q[3];
sx q[3];
rz(-1.7226487) q[3];
sx q[3];
rz(0.050475033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5807895) q[0];
sx q[0];
rz(-2.2152948) q[0];
sx q[0];
rz(-1.3128989) q[0];
rz(1.1543697) q[1];
sx q[1];
rz(-1.9998974) q[1];
sx q[1];
rz(2.5837574) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2247705) q[0];
sx q[0];
rz(-1.0616515) q[0];
sx q[0];
rz(2.1780685) q[0];
x q[1];
rz(0.16061546) q[2];
sx q[2];
rz(-1.4945507) q[2];
sx q[2];
rz(0.74197021) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5231042) q[1];
sx q[1];
rz(-2.2918211) q[1];
sx q[1];
rz(-2.7343661) q[1];
x q[2];
rz(-1.4194896) q[3];
sx q[3];
rz(-1.1137205) q[3];
sx q[3];
rz(0.099927038) q[3];
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
rz(0.15527209) q[3];
sx q[3];
rz(-2.1657491) q[3];
sx q[3];
rz(-0.37080216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3575386) q[0];
sx q[0];
rz(-0.58670601) q[0];
sx q[0];
rz(-2.6724755) q[0];
rz(-2.6672089) q[1];
sx q[1];
rz(-2.3553039) q[1];
sx q[1];
rz(2.3862086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84635272) q[0];
sx q[0];
rz(-1.2995412) q[0];
sx q[0];
rz(-3.1370107) q[0];
rz(-0.44122656) q[2];
sx q[2];
rz(-0.6775113) q[2];
sx q[2];
rz(0.73523075) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1052221) q[1];
sx q[1];
rz(-1.7835938) q[1];
sx q[1];
rz(3.0119621) q[1];
x q[2];
rz(-0.83306082) q[3];
sx q[3];
rz(-1.2445881) q[3];
sx q[3];
rz(-1.8859832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0581806) q[2];
sx q[2];
rz(-1.3231134) q[2];
sx q[2];
rz(0.18784909) q[2];
rz(1.2457054) q[3];
sx q[3];
rz(-1.7626423) q[3];
sx q[3];
rz(-2.0511621) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.93194) q[0];
sx q[0];
rz(-1.8745475) q[0];
sx q[0];
rz(-0.39091045) q[0];
rz(-0.93005013) q[1];
sx q[1];
rz(-1.7101733) q[1];
sx q[1];
rz(1.3040868) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6041038) q[0];
sx q[0];
rz(-1.9216282) q[0];
sx q[0];
rz(-0.022401698) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34159625) q[2];
sx q[2];
rz(-1.4899906) q[2];
sx q[2];
rz(-2.5668457) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.48941055) q[1];
sx q[1];
rz(-2.6565108) q[1];
sx q[1];
rz(1.9368383) q[1];
rz(-pi) q[2];
x q[2];
rz(0.061873925) q[3];
sx q[3];
rz(-0.29682595) q[3];
sx q[3];
rz(2.5421531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2685252) q[2];
sx q[2];
rz(-2.8748942) q[2];
sx q[2];
rz(3.0282057) q[2];
rz(-0.015965613) q[3];
sx q[3];
rz(-1.6458052) q[3];
sx q[3];
rz(0.16460831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.78405821) q[0];
sx q[0];
rz(-1.6679732) q[0];
sx q[0];
rz(3.0175324) q[0];
rz(3.0198174) q[1];
sx q[1];
rz(-0.75141326) q[1];
sx q[1];
rz(-1.6990936) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0826113) q[0];
sx q[0];
rz(-1.9537203) q[0];
sx q[0];
rz(0.10848606) q[0];
rz(1.6142443) q[2];
sx q[2];
rz(-1.8624616) q[2];
sx q[2];
rz(2.0558321) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5296764) q[1];
sx q[1];
rz(-1.2280591) q[1];
sx q[1];
rz(3.0004098) q[1];
rz(-pi) q[2];
rz(0.7300709) q[3];
sx q[3];
rz(-1.3833117) q[3];
sx q[3];
rz(1.2227675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1759935) q[2];
sx q[2];
rz(-1.3600574) q[2];
sx q[2];
rz(2.3967801) q[2];
rz(2.2143769) q[3];
sx q[3];
rz(-1.5085446) q[3];
sx q[3];
rz(-2.0492699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.6532779) q[0];
sx q[0];
rz(-1.4396242) q[0];
sx q[0];
rz(2.9162245) q[0];
rz(-2.0291746) q[1];
sx q[1];
rz(-2.367159) q[1];
sx q[1];
rz(-0.89881277) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32466896) q[0];
sx q[0];
rz(-2.3153439) q[0];
sx q[0];
rz(2.3284495) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4821539) q[2];
sx q[2];
rz(-2.5213679) q[2];
sx q[2];
rz(-2.8218215) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1370217) q[1];
sx q[1];
rz(-2.4413476) q[1];
sx q[1];
rz(-1.2546468) q[1];
rz(1.6055272) q[3];
sx q[3];
rz(-2.8485549) q[3];
sx q[3];
rz(-1.422783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.044346873) q[2];
sx q[2];
rz(-1.4298507) q[2];
sx q[2];
rz(-0.35935768) q[2];
rz(-2.8906631) q[3];
sx q[3];
rz(-2.0693306) q[3];
sx q[3];
rz(0.76771626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74442416) q[0];
sx q[0];
rz(-0.1305307) q[0];
sx q[0];
rz(-0.8771483) q[0];
rz(1.5769222) q[1];
sx q[1];
rz(-1.501187) q[1];
sx q[1];
rz(2.3559779) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3451259) q[0];
sx q[0];
rz(-2.8166127) q[0];
sx q[0];
rz(-0.82635211) q[0];
x q[1];
rz(0.65746515) q[2];
sx q[2];
rz(-0.6415002) q[2];
sx q[2];
rz(-2.7335707) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12745107) q[1];
sx q[1];
rz(-1.9147062) q[1];
sx q[1];
rz(-1.4959072) q[1];
x q[2];
rz(0.26667554) q[3];
sx q[3];
rz(-2.6196822) q[3];
sx q[3];
rz(-0.44346186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2239573) q[2];
sx q[2];
rz(-2.3515297) q[2];
sx q[2];
rz(2.4033974) q[2];
rz(-2.2186642) q[3];
sx q[3];
rz(-1.3083369) q[3];
sx q[3];
rz(0.2963399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7947163) q[0];
sx q[0];
rz(-1.3608169) q[0];
sx q[0];
rz(0.97186744) q[0];
rz(2.234266) q[1];
sx q[1];
rz(-2.0569888) q[1];
sx q[1];
rz(2.8203698) q[1];
rz(2.1512866) q[2];
sx q[2];
rz(-1.2377501) q[2];
sx q[2];
rz(-0.35828423) q[2];
rz(0.76413705) q[3];
sx q[3];
rz(-0.69307477) q[3];
sx q[3];
rz(0.57144036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
