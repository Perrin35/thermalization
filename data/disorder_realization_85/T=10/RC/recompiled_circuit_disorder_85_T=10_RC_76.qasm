OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(-2.5512295) q[0];
sx q[0];
rz(-0.37101775) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(-1.7655656) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3289514) q[0];
sx q[0];
rz(-2.0094123) q[0];
sx q[0];
rz(-2.2493275) q[0];
x q[1];
rz(1.0797834) q[2];
sx q[2];
rz(-0.91677374) q[2];
sx q[2];
rz(-2.430928) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8810597) q[1];
sx q[1];
rz(-0.36306371) q[1];
sx q[1];
rz(-2.0102324) q[1];
x q[2];
rz(1.6269496) q[3];
sx q[3];
rz(-2.3893642) q[3];
sx q[3];
rz(2.5889531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0573037) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(-0.98891813) q[2];
rz(2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(0.83077103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.20724021) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(1.9616615) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(0.72431272) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2253101) q[0];
sx q[0];
rz(-2.3819469) q[0];
sx q[0];
rz(-2.0347974) q[0];
rz(-pi) q[1];
rz(-2.9421147) q[2];
sx q[2];
rz(-1.6530767) q[2];
sx q[2];
rz(-0.85180887) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7608632) q[1];
sx q[1];
rz(-2.7298096) q[1];
sx q[1];
rz(2.1058583) q[1];
rz(-pi) q[2];
rz(-1.768126) q[3];
sx q[3];
rz(-1.8036246) q[3];
sx q[3];
rz(-2.3372834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2362242) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(2.779707) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(-3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996465) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(1.3954337) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(1.9225072) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65788668) q[0];
sx q[0];
rz(-1.3077967) q[0];
sx q[0];
rz(1.9348683) q[0];
rz(-pi) q[1];
rz(-1.2246386) q[2];
sx q[2];
rz(-2.4097754) q[2];
sx q[2];
rz(0.099230448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3238941) q[1];
sx q[1];
rz(-0.59191275) q[1];
sx q[1];
rz(-1.1281668) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5145281) q[3];
sx q[3];
rz(-0.26841044) q[3];
sx q[3];
rz(2.066582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7200155) q[2];
sx q[2];
rz(-2.4013459) q[2];
sx q[2];
rz(1.6960309) q[2];
rz(-0.56882632) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(-0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9179984) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(2.6233327) q[0];
rz(2.4261684) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(0.82675654) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8985734) q[0];
sx q[0];
rz(-1.1613701) q[0];
sx q[0];
rz(-0.55222521) q[0];
rz(-2.3739359) q[2];
sx q[2];
rz(-1.4738238) q[2];
sx q[2];
rz(-1.8061639) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0930867) q[1];
sx q[1];
rz(-2.6566681) q[1];
sx q[1];
rz(2.4946458) q[1];
rz(-0.34927807) q[3];
sx q[3];
rz(-0.75540245) q[3];
sx q[3];
rz(-2.3625771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9892019) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(3.0692696) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3146661) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(-1.1258874) q[0];
rz(-0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(2.856423) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1681686) q[0];
sx q[0];
rz(-2.0452721) q[0];
sx q[0];
rz(-0.04767496) q[0];
rz(2.4777849) q[2];
sx q[2];
rz(-1.2592578) q[2];
sx q[2];
rz(2.2500492) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73357108) q[1];
sx q[1];
rz(-1.5703778) q[1];
sx q[1];
rz(-1.2577406) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6676335) q[3];
sx q[3];
rz(-1.1708461) q[3];
sx q[3];
rz(-1.6638343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29331648) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(2.6021393) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(-2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25813112) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(-0.15701292) q[0];
rz(-0.69333386) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(-1.3670115) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1098423) q[0];
sx q[0];
rz(-1.6037918) q[0];
sx q[0];
rz(-1.6389636) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81142601) q[2];
sx q[2];
rz(-2.0138513) q[2];
sx q[2];
rz(-1.7048938) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4945592) q[1];
sx q[1];
rz(-1.651598) q[1];
sx q[1];
rz(-1.2809491) q[1];
rz(-1.6744162) q[3];
sx q[3];
rz(-1.9220256) q[3];
sx q[3];
rz(0.87408376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(-0.40346754) q[2];
rz(-2.6599595) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(-0.51923716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8994609) q[0];
sx q[0];
rz(-0.88328981) q[0];
sx q[0];
rz(-2.2677299) q[0];
rz(0.44772398) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(1.1706932) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8255575) q[0];
sx q[0];
rz(-0.023840126) q[0];
sx q[0];
rz(-0.26160474) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11622073) q[2];
sx q[2];
rz(-1.6214317) q[2];
sx q[2];
rz(1.2223787) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16329855) q[1];
sx q[1];
rz(-2.1566609) q[1];
sx q[1];
rz(2.3818124) q[1];
rz(-pi) q[2];
rz(-0.53415926) q[3];
sx q[3];
rz(-1.307752) q[3];
sx q[3];
rz(2.0260889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0447023) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(2.5308385) q[2];
rz(2.6664873) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24191813) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(2.8444667) q[0];
rz(1.3946474) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(2.4954605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6883811) q[0];
sx q[0];
rz(-1.4697945) q[0];
sx q[0];
rz(-1.3636916) q[0];
rz(-pi) q[1];
rz(1.4023151) q[2];
sx q[2];
rz(-2.0373166) q[2];
sx q[2];
rz(-0.70665765) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7039973) q[1];
sx q[1];
rz(-2.5181209) q[1];
sx q[1];
rz(-0.61203476) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5313247) q[3];
sx q[3];
rz(-2.2329997) q[3];
sx q[3];
rz(-0.83412795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0769161) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(0.45483744) q[2];
rz(-2.440195) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69944537) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(2.3440857) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(-0.10841766) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71589564) q[0];
sx q[0];
rz(-1.6285537) q[0];
sx q[0];
rz(-0.9066559) q[0];
rz(-pi) q[1];
rz(0.046594521) q[2];
sx q[2];
rz(-1.3923044) q[2];
sx q[2];
rz(2.5537234) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0260967) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(0.91067578) q[1];
rz(2.2750862) q[3];
sx q[3];
rz(-1.7139072) q[3];
sx q[3];
rz(-2.0742311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-2.8016395) q[2];
rz(-0.41839504) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-0.72559124) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6091992) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(0.6788196) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(0.055158786) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63294166) q[0];
sx q[0];
rz(-1.635449) q[0];
sx q[0];
rz(2.7642194) q[0];
x q[1];
rz(1.360838) q[2];
sx q[2];
rz(-1.0210438) q[2];
sx q[2];
rz(1.9090261) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.058051) q[1];
sx q[1];
rz(-1.0293048) q[1];
sx q[1];
rz(0.70152775) q[1];
rz(-pi) q[2];
rz(2.0268029) q[3];
sx q[3];
rz(-1.5378693) q[3];
sx q[3];
rz(1.1820716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(0.49017635) q[2];
rz(-0.13752078) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(2.2035051) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4162083) q[0];
sx q[0];
rz(-1.890247) q[0];
sx q[0];
rz(2.4631137) q[0];
rz(-0.2086808) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(-2.8293777) q[2];
sx q[2];
rz(-1.4785462) q[2];
sx q[2];
rz(-1.2350456) q[2];
rz(1.6155852) q[3];
sx q[3];
rz(-1.8106034) q[3];
sx q[3];
rz(-2.8869224) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
