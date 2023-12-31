OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.934259) q[0];
sx q[0];
rz(-0.59036314) q[0];
sx q[0];
rz(-2.7705749) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(1.376027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24297548) q[0];
sx q[0];
rz(-0.78865047) q[0];
sx q[0];
rz(0.9289766) q[0];
rz(-0.55144989) q[2];
sx q[2];
rz(-2.3460238) q[2];
sx q[2];
rz(-1.7099107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7264759) q[1];
sx q[1];
rz(-1.8980025) q[1];
sx q[1];
rz(0.16023689) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.052470603) q[3];
sx q[3];
rz(-0.82004181) q[3];
sx q[3];
rz(-2.6657871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(-0.98891813) q[2];
rz(-2.3890498) q[3];
sx q[3];
rz(-1.9957333) q[3];
sx q[3];
rz(0.83077103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(1.9616615) q[0];
rz(-2.143899) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(0.72431272) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2253101) q[0];
sx q[0];
rz(-2.3819469) q[0];
sx q[0];
rz(-2.0347974) q[0];
x q[1];
rz(-1.6547336) q[2];
sx q[2];
rz(-1.7695904) q[2];
sx q[2];
rz(-2.4059911) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95485479) q[1];
sx q[1];
rz(-1.2192982) q[1];
sx q[1];
rz(2.9224706) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3734666) q[3];
sx q[3];
rz(-1.337968) q[3];
sx q[3];
rz(-0.80430921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.90536845) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(-0.36188564) q[2];
rz(0.13606717) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(-0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996465) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(-1.746159) q[0];
rz(0.46229258) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(-1.9225072) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8273979) q[0];
sx q[0];
rz(-0.44566804) q[0];
sx q[0];
rz(-2.218194) q[0];
rz(-2.2723324) q[2];
sx q[2];
rz(-1.3420891) q[2];
sx q[2];
rz(-1.2094091) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8428426) q[1];
sx q[1];
rz(-2.0992273) q[1];
sx q[1];
rz(2.8612086) q[1];
x q[2];
rz(3.1261256) q[3];
sx q[3];
rz(-1.3028212) q[3];
sx q[3];
rz(1.0166575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42157713) q[2];
sx q[2];
rz(-2.4013459) q[2];
sx q[2];
rz(-1.6960309) q[2];
rz(0.56882632) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(-3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.9179984) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(-0.51825994) q[0];
rz(0.7154243) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(2.3148361) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0543538) q[0];
sx q[0];
rz(-2.0728489) q[0];
sx q[0];
rz(-1.0994083) q[0];
rz(-pi) q[1];
rz(0.13917285) q[2];
sx q[2];
rz(-2.3690802) q[2];
sx q[2];
rz(0.13538361) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7552232) q[1];
sx q[1];
rz(-1.9519023) q[1];
sx q[1];
rz(-1.2632881) q[1];
rz(-pi) q[2];
rz(0.72439648) q[3];
sx q[3];
rz(-1.8076234) q[3];
sx q[3];
rz(-0.5326007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9892019) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(-1.3789122) q[2];
rz(-0.072323024) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(-1.6453843) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3146661) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(-1.1258874) q[0];
rz(2.2391438) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(0.28516969) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42442214) q[0];
sx q[0];
rz(-1.6132014) q[0];
sx q[0];
rz(1.0958584) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4816149) q[2];
sx q[2];
rz(-2.4184605) q[2];
sx q[2];
rz(2.8358104) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3030745) q[1];
sx q[1];
rz(-2.8285366) q[1];
sx q[1];
rz(1.5721553) q[1];
x q[2];
rz(0.40163715) q[3];
sx q[3];
rz(-1.65997) q[3];
sx q[3];
rz(3.0107486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29331648) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(-2.6021393) q[2];
rz(-0.30682492) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.25813112) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(0.15701292) q[0];
rz(2.4482588) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(1.7745811) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6002944) q[0];
sx q[0];
rz(-1.6389264) q[0];
sx q[0];
rz(-0.033072254) q[0];
rz(-pi) q[1];
rz(-2.1742646) q[2];
sx q[2];
rz(-2.2852995) q[2];
sx q[2];
rz(-2.8514903) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4821266) q[1];
sx q[1];
rz(-2.8409993) q[1];
sx q[1];
rz(-1.8468922) q[1];
x q[2];
rz(-1.6744162) q[3];
sx q[3];
rz(-1.9220256) q[3];
sx q[3];
rz(-2.2675089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8453025) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(2.7381251) q[2];
rz(-2.6599595) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-2.2677299) q[0];
rz(0.44772398) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(1.1706932) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3160352) q[0];
sx q[0];
rz(-0.023840126) q[0];
sx q[0];
rz(-0.26160474) q[0];
x q[1];
rz(1.621775) q[2];
sx q[2];
rz(-1.4547252) q[2];
sx q[2];
rz(-0.34250868) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9782941) q[1];
sx q[1];
rz(-0.98493176) q[1];
sx q[1];
rz(0.75978029) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6550754) q[3];
sx q[3];
rz(-0.58972893) q[3];
sx q[3];
rz(-0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0447023) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(2.5308385) q[2];
rz(0.47510535) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(-0.92774123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.8996745) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(-0.29712594) q[0];
rz(1.3946474) q[1];
sx q[1];
rz(-1.1509832) q[1];
sx q[1];
rz(0.64613211) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0028249) q[0];
sx q[0];
rz(-1.3647623) q[0];
sx q[0];
rz(0.10319184) q[0];
rz(-1.4023151) q[2];
sx q[2];
rz(-1.104276) q[2];
sx q[2];
rz(2.434935) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7568126) q[1];
sx q[1];
rz(-1.2287178) q[1];
sx q[1];
rz(2.6095819) q[1];
rz(0.93692245) q[3];
sx q[3];
rz(-2.2735032) q[3];
sx q[3];
rz(0.015451775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.064676553) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(-0.45483744) q[2];
rz(-2.440195) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(1.1340244) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4421473) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(-0.79750693) q[0];
rz(-0.51756716) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(-3.033175) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3603044) q[0];
sx q[0];
rz(-2.4753248) q[0];
sx q[0];
rz(1.4772619) q[0];
rz(-pi) q[1];
rz(1.3181395) q[2];
sx q[2];
rz(-2.9571819) q[2];
sx q[2];
rz(0.84469634) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0260967) q[1];
sx q[1];
rz(-2.8671088) q[1];
sx q[1];
rz(2.2309169) q[1];
rz(-pi) q[2];
rz(1.3518203) q[3];
sx q[3];
rz(-2.4253546) q[3];
sx q[3];
rz(-2.4718474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-2.8016395) q[2];
rz(-2.7231976) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(-0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6091992) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(-2.4627731) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(-0.055158786) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77627742) q[0];
sx q[0];
rz(-0.3826097) q[0];
sx q[0];
rz(2.9676653) q[0];
rz(-pi) q[1];
rz(-1.360838) q[2];
sx q[2];
rz(-1.0210438) q[2];
sx q[2];
rz(-1.9090261) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89796472) q[1];
sx q[1];
rz(-0.98476714) q[1];
sx q[1];
rz(2.2378053) q[1];
x q[2];
rz(-0.03667128) q[3];
sx q[3];
rz(-2.026537) q[3];
sx q[3];
rz(-2.7367221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1577592) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(-0.49017635) q[2];
rz(0.13752078) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4162083) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(-2.9329119) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(0.31221496) q[2];
sx q[2];
rz(-1.4785462) q[2];
sx q[2];
rz(-1.2350456) q[2];
rz(2.9604838) q[3];
sx q[3];
rz(-2.8977179) q[3];
sx q[3];
rz(0.44117622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
