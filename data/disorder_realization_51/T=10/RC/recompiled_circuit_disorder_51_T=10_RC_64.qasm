OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5470619) q[0];
sx q[0];
rz(4.2630258) q[0];
sx q[0];
rz(6.1017258) q[0];
rz(2.060086) q[1];
sx q[1];
rz(-0.67343155) q[1];
sx q[1];
rz(-1.0531309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15381972) q[0];
sx q[0];
rz(-1.1115371) q[0];
sx q[0];
rz(2.905373) q[0];
rz(-pi) q[1];
rz(1.8671145) q[2];
sx q[2];
rz(-2.1726492) q[2];
sx q[2];
rz(1.5822496) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3095113) q[1];
sx q[1];
rz(-2.32825) q[1];
sx q[1];
rz(1.8413715) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1374723) q[3];
sx q[3];
rz(-2.7702799) q[3];
sx q[3];
rz(-2.5759047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.16333214) q[2];
sx q[2];
rz(-1.1085199) q[2];
sx q[2];
rz(-1.367761) q[2];
rz(2.1286428) q[3];
sx q[3];
rz(-0.84665853) q[3];
sx q[3];
rz(-3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.098009) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(0.59511551) q[0];
rz(1.0455421) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(1.5140623) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1925416) q[0];
sx q[0];
rz(-2.4727614) q[0];
sx q[0];
rz(2.2632954) q[0];
rz(-2.3081231) q[2];
sx q[2];
rz(-1.5916628) q[2];
sx q[2];
rz(1.8881063) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.086386911) q[1];
sx q[1];
rz(-1.7654164) q[1];
sx q[1];
rz(-0.17533949) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1176923) q[3];
sx q[3];
rz(-3.0278904) q[3];
sx q[3];
rz(-2.5642455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.65511584) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(-0.79616037) q[2];
rz(2.1697309) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(2.9698353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3765091) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(-3.0774975) q[0];
rz(0.31072101) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(1.6832738) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9653939) q[0];
sx q[0];
rz(-1.6178693) q[0];
sx q[0];
rz(-1.7182299) q[0];
x q[1];
rz(-2.4887423) q[2];
sx q[2];
rz(-2.0712426) q[2];
sx q[2];
rz(2.1476538) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7160733) q[1];
sx q[1];
rz(-2.7287373) q[1];
sx q[1];
rz(-0.13456657) q[1];
x q[2];
rz(3.0890373) q[3];
sx q[3];
rz(-0.31523809) q[3];
sx q[3];
rz(1.0759575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1594499) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(0.67908755) q[2];
rz(-1.7689765) q[3];
sx q[3];
rz(-1.8434098) q[3];
sx q[3];
rz(2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.8452334) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(-1.1244208) q[0];
rz(1.9212978) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(1.6569998) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.785448) q[0];
sx q[0];
rz(-1.5404285) q[0];
sx q[0];
rz(1.6214451) q[0];
x q[1];
rz(1.053327) q[2];
sx q[2];
rz(-1.7382858) q[2];
sx q[2];
rz(0.53211624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0255819) q[1];
sx q[1];
rz(-0.75041795) q[1];
sx q[1];
rz(-2.2376899) q[1];
rz(-pi) q[2];
rz(-0.72317601) q[3];
sx q[3];
rz(-1.3190862) q[3];
sx q[3];
rz(0.46079208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0174039) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(-1.0162639) q[2];
rz(1.4034363) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6338585) q[0];
sx q[0];
rz(-1.8554747) q[0];
sx q[0];
rz(0.39988363) q[0];
rz(1.9790861) q[1];
sx q[1];
rz(-1.3299273) q[1];
sx q[1];
rz(0.17366017) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7090209) q[0];
sx q[0];
rz(-3.0637494) q[0];
sx q[0];
rz(2.7709333) q[0];
rz(-1.3599672) q[2];
sx q[2];
rz(-1.9569009) q[2];
sx q[2];
rz(-2.2112276) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.793321) q[1];
sx q[1];
rz(-1.6760567) q[1];
sx q[1];
rz(-3.1131016) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33186121) q[3];
sx q[3];
rz(-1.3734986) q[3];
sx q[3];
rz(-2.8062537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6614723) q[2];
sx q[2];
rz(-1.1602594) q[2];
sx q[2];
rz(-0.76812569) q[2];
rz(0.85401946) q[3];
sx q[3];
rz(-1.4203527) q[3];
sx q[3];
rz(2.1693726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0356692) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(1.6712028) q[0];
rz(-0.51180965) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(-1.9981729) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2238732) q[0];
sx q[0];
rz(-1.6544764) q[0];
sx q[0];
rz(2.9907945) q[0];
rz(-0.60322275) q[2];
sx q[2];
rz(-1.2957186) q[2];
sx q[2];
rz(3.0657257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1089576) q[1];
sx q[1];
rz(-2.7320478) q[1];
sx q[1];
rz(-1.1213379) q[1];
rz(-pi) q[2];
rz(2.1120758) q[3];
sx q[3];
rz(-1.8347077) q[3];
sx q[3];
rz(-1.2332066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64289552) q[2];
sx q[2];
rz(-0.84474793) q[2];
sx q[2];
rz(-1.6112304) q[2];
rz(1.4536084) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(-0.24916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11944184) q[0];
sx q[0];
rz(-0.74437737) q[0];
sx q[0];
rz(-1.7013593) q[0];
rz(2.4123736) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(-1.1332606) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6183375) q[0];
sx q[0];
rz(-1.528307) q[0];
sx q[0];
rz(-0.32251127) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.557017) q[2];
sx q[2];
rz(-0.85859495) q[2];
sx q[2];
rz(0.53879246) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.28753528) q[1];
sx q[1];
rz(-0.71811986) q[1];
sx q[1];
rz(1.992804) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8797638) q[3];
sx q[3];
rz(-1.7228408) q[3];
sx q[3];
rz(0.1482299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7363654) q[2];
sx q[2];
rz(-1.1232802) q[2];
sx q[2];
rz(0.48842946) q[2];
rz(1.8296261) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(0.64129889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064780386) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(3.0704165) q[0];
rz(3.1094303) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(-1.2088998) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8457444) q[0];
sx q[0];
rz(-0.83202067) q[0];
sx q[0];
rz(0.31006281) q[0];
x q[1];
rz(-1.3337295) q[2];
sx q[2];
rz(-1.7500688) q[2];
sx q[2];
rz(0.2722424) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80450088) q[1];
sx q[1];
rz(-3.0102804) q[1];
sx q[1];
rz(-1.6930224) q[1];
x q[2];
rz(0.62187059) q[3];
sx q[3];
rz(-1.8727881) q[3];
sx q[3];
rz(-2.5627476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20748392) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(-1.5709546) q[2];
rz(2.2682244) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(-0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97380012) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(-2.8299676) q[0];
rz(2.3198126) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(-1.6315546) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2841543) q[0];
sx q[0];
rz(-0.3785924) q[0];
sx q[0];
rz(2.5749102) q[0];
rz(-pi) q[1];
rz(2.0779583) q[2];
sx q[2];
rz(-0.38804752) q[2];
sx q[2];
rz(1.0840814) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1875302) q[1];
sx q[1];
rz(-2.4937951) q[1];
sx q[1];
rz(2.8026583) q[1];
rz(-1.7543091) q[3];
sx q[3];
rz(-1.3418875) q[3];
sx q[3];
rz(-2.6218417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8682378) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(-2.3256425) q[2];
rz(-2.6319035) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.2508535) q[0];
sx q[0];
rz(-1.7128523) q[0];
sx q[0];
rz(-0.2302641) q[0];
rz(-0.62581217) q[1];
sx q[1];
rz(-0.9451378) q[1];
sx q[1];
rz(-0.65840107) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2021159) q[0];
sx q[0];
rz(-2.3411223) q[0];
sx q[0];
rz(0.70864422) q[0];
rz(0.41215956) q[2];
sx q[2];
rz(-0.71752749) q[2];
sx q[2];
rz(-1.621643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2681594) q[1];
sx q[1];
rz(-0.50926103) q[1];
sx q[1];
rz(-1.9164273) q[1];
rz(-pi) q[2];
rz(-0.53604605) q[3];
sx q[3];
rz(-1.2671766) q[3];
sx q[3];
rz(2.924078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6752424) q[2];
sx q[2];
rz(-0.8224951) q[2];
sx q[2];
rz(2.8038483) q[2];
rz(1.0234458) q[3];
sx q[3];
rz(-1.7815536) q[3];
sx q[3];
rz(-1.9256928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6170549) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(1.9032003) q[1];
sx q[1];
rz(-1.4307784) q[1];
sx q[1];
rz(1.4684114) q[1];
rz(-0.17910437) q[2];
sx q[2];
rz(-2.1894107) q[2];
sx q[2];
rz(-1.4084963) q[2];
rz(-1.5068549) q[3];
sx q[3];
rz(-2.1422374) q[3];
sx q[3];
rz(-0.16850785) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];