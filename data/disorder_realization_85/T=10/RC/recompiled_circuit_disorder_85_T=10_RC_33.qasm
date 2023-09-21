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
rz(2.5420904) q[1];
sx q[1];
rz(11.190344) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8986172) q[0];
sx q[0];
rz(-0.78865047) q[0];
sx q[0];
rz(-0.9289766) q[0];
x q[1];
rz(2.5901428) q[2];
sx q[2];
rz(-0.79556888) q[2];
sx q[2];
rz(-1.431682) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4151167) q[1];
sx q[1];
rz(-1.2435902) q[1];
sx q[1];
rz(0.16023689) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6269496) q[3];
sx q[3];
rz(-2.3893642) q[3];
sx q[3];
rz(-0.55263954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0573037) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(-2.1526745) q[2];
rz(2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(0.83077103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(-1.1799312) q[0];
rz(-0.99769366) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(0.72431272) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82942078) q[0];
sx q[0];
rz(-2.2342626) q[0];
sx q[0];
rz(-2.7396766) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39436491) q[2];
sx q[2];
rz(-0.21557237) q[2];
sx q[2];
rz(-2.8087316) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38072941) q[1];
sx q[1];
rz(-0.41178307) q[1];
sx q[1];
rz(2.1058583) q[1];
x q[2];
rz(-0.69085391) q[3];
sx q[3];
rz(-2.8375531) q[3];
sx q[3];
rz(-0.090304852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1419462) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(1.3954337) q[0];
rz(0.46229258) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(-1.9225072) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.483706) q[0];
sx q[0];
rz(-1.833796) q[0];
sx q[0];
rz(1.9348683) q[0];
rz(-pi) q[1];
rz(-0.29580446) q[2];
sx q[2];
rz(-2.25053) q[2];
sx q[2];
rz(-2.5909397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8176986) q[1];
sx q[1];
rz(-0.59191275) q[1];
sx q[1];
rz(1.1281668) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.015467042) q[3];
sx q[3];
rz(-1.8387715) q[3];
sx q[3];
rz(-1.0166575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.4455618) q[2];
rz(2.5727663) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9179984) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(0.51825994) q[0];
rz(0.7154243) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(2.3148361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8985734) q[0];
sx q[0];
rz(-1.1613701) q[0];
sx q[0];
rz(2.5893674) q[0];
rz(-pi) q[1];
rz(0.76765676) q[2];
sx q[2];
rz(-1.4738238) q[2];
sx q[2];
rz(-1.8061639) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0930867) q[1];
sx q[1];
rz(-2.6566681) q[1];
sx q[1];
rz(0.64694689) q[1];
rz(-pi) q[2];
rz(-0.34927807) q[3];
sx q[3];
rz(-2.3861902) q[3];
sx q[3];
rz(-0.7790156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9892019) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(1.3789122) q[2];
rz(-3.0692696) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(-1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82692659) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(2.0157053) q[0];
rz(0.90244883) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(-0.28516969) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0641159) q[0];
sx q[0];
rz(-0.47668326) q[0];
sx q[0];
rz(1.478273) q[0];
x q[1];
rz(-0.6638078) q[2];
sx q[2];
rz(-1.8823349) q[2];
sx q[2];
rz(0.89154348) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73357108) q[1];
sx q[1];
rz(-1.5703778) q[1];
sx q[1];
rz(1.8838521) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6676335) q[3];
sx q[3];
rz(-1.1708461) q[3];
sx q[3];
rz(-1.4777583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29331648) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(0.53945333) q[2];
rz(-2.8347677) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25813112) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(2.4482588) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(1.7745811) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6002944) q[0];
sx q[0];
rz(-1.6389264) q[0];
sx q[0];
rz(3.1085204) q[0];
rz(2.1742646) q[2];
sx q[2];
rz(-0.8562932) q[2];
sx q[2];
rz(-2.8514903) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4821266) q[1];
sx q[1];
rz(-0.30059338) q[1];
sx q[1];
rz(1.8468922) q[1];
x q[2];
rz(0.27512392) q[3];
sx q[3];
rz(-2.776006) q[3];
sx q[3];
rz(-2.5610353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(-2.7381251) q[2];
rz(0.48163313) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(2.6223555) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24213174) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(0.8738628) q[0];
rz(2.6938687) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(1.9708995) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6252977) q[0];
sx q[0];
rz(-1.5769616) q[0];
sx q[0];
rz(3.1185634) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0253719) q[2];
sx q[2];
rz(-1.5201609) q[2];
sx q[2];
rz(-1.2223787) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2504187) q[1];
sx q[1];
rz(-0.95953566) q[1];
sx q[1];
rz(-2.3120018) q[1];
x q[2];
rz(-0.4865173) q[3];
sx q[3];
rz(-0.58972893) q[3];
sx q[3];
rz(0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0447023) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(0.61075413) q[2];
rz(2.6664873) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(0.92774123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8996745) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(-0.29712594) q[0];
rz(-1.7469453) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(-0.64613211) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4532115) q[0];
sx q[0];
rz(-1.6717981) q[0];
sx q[0];
rz(-1.7779011) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6693194) q[2];
sx q[2];
rz(-1.7211282) q[2];
sx q[2];
rz(-0.78778247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1505193) q[1];
sx q[1];
rz(-1.0725613) q[1];
sx q[1];
rz(-1.9626161) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93692245) q[3];
sx q[3];
rz(-2.2735032) q[3];
sx q[3];
rz(0.015451775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.064676553) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(-2.440195) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4421473) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(0.79750693) q[0];
rz(-2.6240255) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(-3.033175) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90005504) q[0];
sx q[0];
rz(-2.2336322) q[0];
sx q[0];
rz(3.0682949) q[0];
rz(1.8234532) q[2];
sx q[2];
rz(-2.9571819) q[2];
sx q[2];
rz(-0.84469634) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79418102) q[1];
sx q[1];
rz(-1.7865744) q[1];
sx q[1];
rz(-2.9706035) q[1];
rz(-pi) q[2];
rz(-0.18687825) q[3];
sx q[3];
rz(-2.2664321) q[3];
sx q[3];
rz(2.7587492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0124399) q[2];
sx q[2];
rz(-1.3468578) q[2];
sx q[2];
rz(-0.33995315) q[2];
rz(2.7231976) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(2.4160014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5323935) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(2.4627731) q[0];
rz(-0.36418307) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(-0.055158786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63294166) q[0];
sx q[0];
rz(-1.635449) q[0];
sx q[0];
rz(2.7642194) q[0];
rz(-pi) q[1];
rz(-2.5818985) q[2];
sx q[2];
rz(-1.392138) q[2];
sx q[2];
rz(0.44911227) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.060505796) q[1];
sx q[1];
rz(-0.85716893) q[1];
sx q[1];
rz(-2.39141) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4961365) q[3];
sx q[3];
rz(-2.6844822) q[3];
sx q[3];
rz(-0.32170579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.98383343) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(-0.49017635) q[2];
rz(-0.13752078) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(0.72538439) q[0];
sx q[0];
rz(-1.890247) q[0];
sx q[0];
rz(2.4631137) q[0];
rz(0.2086808) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(2.8293777) q[2];
sx q[2];
rz(-1.6630465) q[2];
sx q[2];
rz(1.9065471) q[2];
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
