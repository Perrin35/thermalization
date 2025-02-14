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
rz(1.8062502) q[0];
sx q[0];
rz(-0.36646068) q[0];
sx q[0];
rz(-2.6824644) q[0];
rz(-0.65161172) q[1];
sx q[1];
rz(-1.7348644) q[1];
sx q[1];
rz(-2.9098517) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7739959) q[0];
sx q[0];
rz(-1.913842) q[0];
sx q[0];
rz(1.497722) q[0];
rz(-pi) q[1];
rz(3.0475869) q[2];
sx q[2];
rz(-2.5931907) q[2];
sx q[2];
rz(-1.3863871) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5197088) q[1];
sx q[1];
rz(-2.1823898) q[1];
sx q[1];
rz(-2.0382463) q[1];
rz(1.5292589) q[3];
sx q[3];
rz(-1.2327063) q[3];
sx q[3];
rz(-1.5134144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.090652466) q[2];
sx q[2];
rz(-0.53477627) q[2];
sx q[2];
rz(-1.2789307) q[2];
rz(-0.88979641) q[3];
sx q[3];
rz(-1.5225478) q[3];
sx q[3];
rz(-0.33862996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58562529) q[0];
sx q[0];
rz(-2.2323759) q[0];
sx q[0];
rz(2.2157748) q[0];
rz(-2.0766808) q[1];
sx q[1];
rz(-1.5771259) q[1];
sx q[1];
rz(3.056114) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2538214) q[0];
sx q[0];
rz(-1.5740875) q[0];
sx q[0];
rz(-2.7517767) q[0];
rz(1.0930834) q[2];
sx q[2];
rz(-0.6685377) q[2];
sx q[2];
rz(-0.82205176) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.036219941) q[1];
sx q[1];
rz(-2.5565127) q[1];
sx q[1];
rz(-1.1452504) q[1];
rz(-0.55172975) q[3];
sx q[3];
rz(-2.3750616) q[3];
sx q[3];
rz(-0.54099247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3257137) q[2];
sx q[2];
rz(-1.4789944) q[2];
sx q[2];
rz(-0.386664) q[2];
rz(1.9246842) q[3];
sx q[3];
rz(-2.7972126) q[3];
sx q[3];
rz(0.2200505) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81095186) q[0];
sx q[0];
rz(-2.7702259) q[0];
sx q[0];
rz(-0.49705848) q[0];
rz(2.0992384) q[1];
sx q[1];
rz(-0.94756871) q[1];
sx q[1];
rz(-2.846948) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7001273) q[0];
sx q[0];
rz(-1.6396806) q[0];
sx q[0];
rz(2.8089351) q[0];
rz(-pi) q[1];
x q[1];
rz(1.398009) q[2];
sx q[2];
rz(-2.2026988) q[2];
sx q[2];
rz(-1.8586707) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5897419) q[1];
sx q[1];
rz(-0.95847469) q[1];
sx q[1];
rz(-0.53770868) q[1];
rz(-pi) q[2];
rz(2.1717242) q[3];
sx q[3];
rz(-1.6898385) q[3];
sx q[3];
rz(0.62059072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3543388) q[2];
sx q[2];
rz(-2.2806809) q[2];
sx q[2];
rz(-1.5798689) q[2];
rz(-1.076237) q[3];
sx q[3];
rz(-1.302224) q[3];
sx q[3];
rz(-2.2126183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.314986) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(-0.22931799) q[0];
rz(1.5062821) q[1];
sx q[1];
rz(-1.6835667) q[1];
sx q[1];
rz(-0.062072676) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4631043) q[0];
sx q[0];
rz(-1.0764499) q[0];
sx q[0];
rz(0.93017471) q[0];
rz(-pi) q[1];
rz(2.1770085) q[2];
sx q[2];
rz(-0.46614376) q[2];
sx q[2];
rz(-0.073748253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8530555) q[1];
sx q[1];
rz(-1.2934577) q[1];
sx q[1];
rz(-2.4278575) q[1];
rz(0.2009521) q[3];
sx q[3];
rz(-1.4794599) q[3];
sx q[3];
rz(-0.82543711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0486003) q[2];
sx q[2];
rz(-0.91602641) q[2];
sx q[2];
rz(-1.830706) q[2];
rz(0.038289573) q[3];
sx q[3];
rz(-1.7891276) q[3];
sx q[3];
rz(-0.37461764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6467317) q[0];
sx q[0];
rz(-2.4202388) q[0];
sx q[0];
rz(-2.2485961) q[0];
rz(-2.112174) q[1];
sx q[1];
rz(-2.4118377) q[1];
sx q[1];
rz(0.095887862) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3678148) q[0];
sx q[0];
rz(-2.4565963) q[0];
sx q[0];
rz(-0.8789341) q[0];
rz(-pi) q[1];
x q[1];
rz(1.65972) q[2];
sx q[2];
rz(-1.6409573) q[2];
sx q[2];
rz(-1.2360177) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.27285114) q[1];
sx q[1];
rz(-2.4563046) q[1];
sx q[1];
rz(-0.59602316) q[1];
rz(-pi) q[2];
rz(1.9434236) q[3];
sx q[3];
rz(-2.0564787) q[3];
sx q[3];
rz(1.612965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15478495) q[2];
sx q[2];
rz(-1.751535) q[2];
sx q[2];
rz(2.7161157) q[2];
rz(0.42243877) q[3];
sx q[3];
rz(-2.0368302) q[3];
sx q[3];
rz(-2.6384242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4130037) q[0];
sx q[0];
rz(-0.63218963) q[0];
sx q[0];
rz(0.11716209) q[0];
rz(-1.6475742) q[1];
sx q[1];
rz(-2.2393176) q[1];
sx q[1];
rz(-2.07043) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9566222) q[0];
sx q[0];
rz(-1.288646) q[0];
sx q[0];
rz(2.6376702) q[0];
x q[1];
rz(2.6980459) q[2];
sx q[2];
rz(-2.5534592) q[2];
sx q[2];
rz(-3.0556222) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92381645) q[1];
sx q[1];
rz(-2.4144396) q[1];
sx q[1];
rz(1.9132861) q[1];
x q[2];
rz(-1.734524) q[3];
sx q[3];
rz(-1.8843972) q[3];
sx q[3];
rz(-2.0255476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5693207) q[2];
sx q[2];
rz(-0.57491493) q[2];
sx q[2];
rz(-3.102109) q[2];
rz(-0.076400541) q[3];
sx q[3];
rz(-1.1698086) q[3];
sx q[3];
rz(0.45342818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4385248) q[0];
sx q[0];
rz(-0.81145966) q[0];
sx q[0];
rz(-0.95034289) q[0];
rz(-1.1514661) q[1];
sx q[1];
rz(-1.5299503) q[1];
sx q[1];
rz(-1.8340402) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3111585) q[0];
sx q[0];
rz(-1.3518999) q[0];
sx q[0];
rz(-1.9040742) q[0];
rz(2.9588863) q[2];
sx q[2];
rz(-0.38190834) q[2];
sx q[2];
rz(-2.4289102) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3486191) q[1];
sx q[1];
rz(-1.7820616) q[1];
sx q[1];
rz(-1.068685) q[1];
x q[2];
rz(0.96486196) q[3];
sx q[3];
rz(-1.6814702) q[3];
sx q[3];
rz(-1.0718653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3211956) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(0.55366984) q[2];
rz(0.6959483) q[3];
sx q[3];
rz(-1.4724933) q[3];
sx q[3];
rz(1.6863916) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0221136) q[0];
sx q[0];
rz(-1.6895634) q[0];
sx q[0];
rz(2.8346862) q[0];
rz(-1.8652929) q[1];
sx q[1];
rz(-2.6752094) q[1];
sx q[1];
rz(1.2409522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0725482) q[0];
sx q[0];
rz(-2.6676151) q[0];
sx q[0];
rz(1.9182253) q[0];
rz(1.9996793) q[2];
sx q[2];
rz(-2.1253573) q[2];
sx q[2];
rz(2.6311324) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.13682374) q[1];
sx q[1];
rz(-1.3826177) q[1];
sx q[1];
rz(0.31463639) q[1];
x q[2];
rz(-0.60997529) q[3];
sx q[3];
rz(-1.9053359) q[3];
sx q[3];
rz(-1.1776678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3500195) q[2];
sx q[2];
rz(-0.80524033) q[2];
sx q[2];
rz(1.8512858) q[2];
rz(0.6811412) q[3];
sx q[3];
rz(-2.2414424) q[3];
sx q[3];
rz(1.6397938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8660368) q[0];
sx q[0];
rz(-1.0365726) q[0];
sx q[0];
rz(1.9679029) q[0];
rz(0.58700079) q[1];
sx q[1];
rz(-0.78158164) q[1];
sx q[1];
rz(0.079708286) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3143334) q[0];
sx q[0];
rz(-1.9406576) q[0];
sx q[0];
rz(-0.24263675) q[0];
rz(-pi) q[1];
rz(-2.6775421) q[2];
sx q[2];
rz(-2.8949605) q[2];
sx q[2];
rz(2.3753948) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6708128) q[1];
sx q[1];
rz(-0.97320405) q[1];
sx q[1];
rz(-0.87009354) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.856273) q[3];
sx q[3];
rz(-2.3148708) q[3];
sx q[3];
rz(-2.7088425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6012663) q[2];
sx q[2];
rz(-1.2286295) q[2];
sx q[2];
rz(-1.3339174) q[2];
rz(1.5909083) q[3];
sx q[3];
rz(-1.0100789) q[3];
sx q[3];
rz(2.9529115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.659336) q[0];
sx q[0];
rz(-1.5631258) q[0];
sx q[0];
rz(-1.2905066) q[0];
rz(-0.036529649) q[1];
sx q[1];
rz(-1.1643658) q[1];
sx q[1];
rz(2.0972924) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8528906) q[0];
sx q[0];
rz(-1.1452617) q[0];
sx q[0];
rz(0.3279) q[0];
rz(-2.2453111) q[2];
sx q[2];
rz(-0.41322979) q[2];
sx q[2];
rz(-0.21597029) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83936939) q[1];
sx q[1];
rz(-1.4166792) q[1];
sx q[1];
rz(-2.3030512) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32560168) q[3];
sx q[3];
rz(-0.34593098) q[3];
sx q[3];
rz(-0.17445645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2269939) q[2];
sx q[2];
rz(-2.1103766) q[2];
sx q[2];
rz(1.8168137) q[2];
rz(-2.3850208) q[3];
sx q[3];
rz(-0.888266) q[3];
sx q[3];
rz(-2.2911086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.666438) q[0];
sx q[0];
rz(-1.7378687) q[0];
sx q[0];
rz(0.73874656) q[0];
rz(-1.4229763) q[1];
sx q[1];
rz(-1.9444793) q[1];
sx q[1];
rz(0.3107298) q[1];
rz(0.46427609) q[2];
sx q[2];
rz(-2.3515679) q[2];
sx q[2];
rz(-0.49754561) q[2];
rz(-2.4800469) q[3];
sx q[3];
rz(-1.1894124) q[3];
sx q[3];
rz(-0.00581707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
