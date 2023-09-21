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
rz(3.7319558) q[0];
sx q[0];
rz(9.0537602) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(2.5420904) q[1];
sx q[1];
rz(11.190344) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8126412) q[0];
sx q[0];
rz(-1.1321804) q[0];
sx q[0];
rz(-2.2493275) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4260169) q[2];
sx q[2];
rz(-1.1872429) q[2];
sx q[2];
rz(-2.5959612) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.260533) q[1];
sx q[1];
rz(-2.7785289) q[1];
sx q[1];
rz(-1.1313603) q[1];
rz(-pi) q[2];
rz(-1.514643) q[3];
sx q[3];
rz(-2.3893642) q[3];
sx q[3];
rz(-0.55263954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.084289) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(2.1526745) q[2];
rz(2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(0.83077103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(-1.1799312) q[0];
rz(0.99769366) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(-2.4172799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2253101) q[0];
sx q[0];
rz(-0.75964576) q[0];
sx q[0];
rz(-1.1067953) q[0];
rz(-pi) q[1];
rz(-0.39436491) q[2];
sx q[2];
rz(-2.9260203) q[2];
sx q[2];
rz(2.8087316) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.95485479) q[1];
sx q[1];
rz(-1.9222944) q[1];
sx q[1];
rz(-2.9224706) q[1];
x q[2];
rz(-1.768126) q[3];
sx q[3];
rz(-1.8036246) q[3];
sx q[3];
rz(0.80430921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90536845) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(2.779707) q[2];
rz(-0.13606717) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996465) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(1.3954337) q[0];
rz(-2.6793001) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(1.2190855) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0116545) q[0];
sx q[0];
rz(-1.219795) q[0];
sx q[0];
rz(0.28052335) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2246386) q[2];
sx q[2];
rz(-0.73181728) q[2];
sx q[2];
rz(0.099230448) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3238941) q[1];
sx q[1];
rz(-2.5496799) q[1];
sx q[1];
rz(1.1281668) q[1];
x q[2];
rz(-1.6270646) q[3];
sx q[3];
rz(-2.8731822) q[3];
sx q[3];
rz(1.0750107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.6960309) q[2];
rz(-2.5727663) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(0.11051699) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9179984) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(-2.6233327) q[0];
rz(2.4261684) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(0.82675654) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90137705) q[0];
sx q[0];
rz(-2.4670521) q[0];
sx q[0];
rz(-0.69112372) q[0];
x q[1];
rz(1.4364169) q[2];
sx q[2];
rz(-2.3339286) q[2];
sx q[2];
rz(-2.813051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0485059) q[1];
sx q[1];
rz(-2.6566681) q[1];
sx q[1];
rz(-0.64694689) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4171962) q[3];
sx q[3];
rz(-1.3339692) q[3];
sx q[3];
rz(2.608992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15239079) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(-1.3789122) q[2];
rz(-0.072323024) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0641159) q[0];
sx q[0];
rz(-0.47668326) q[0];
sx q[0];
rz(1.478273) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4777849) q[2];
sx q[2];
rz(-1.2592578) q[2];
sx q[2];
rz(-2.2500492) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3045029) q[1];
sx q[1];
rz(-1.883852) q[1];
sx q[1];
rz(-3.1411527) q[1];
x q[2];
rz(0.40163715) q[3];
sx q[3];
rz(-1.65997) q[3];
sx q[3];
rz(-0.13084403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8482762) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(-0.53945333) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(-0.45421281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.8834615) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(0.69333386) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(-1.7745811) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0317504) q[0];
sx q[0];
rz(-1.5378008) q[0];
sx q[0];
rz(1.6389636) q[0];
rz(2.5622257) q[2];
sx q[2];
rz(-0.89951347) q[2];
sx q[2];
rz(2.6210149) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94783084) q[1];
sx q[1];
rz(-1.8596706) q[1];
sx q[1];
rz(-0.08430251) q[1];
x q[2];
rz(1.4671765) q[3];
sx q[3];
rz(-1.9220256) q[3];
sx q[3];
rz(-2.2675089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29629016) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(2.7381251) q[2];
rz(-2.6599595) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(-2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24213174) q[0];
sx q[0];
rz(-0.88328981) q[0];
sx q[0];
rz(-0.8738628) q[0];
rz(2.6938687) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(-1.1706932) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51629492) q[0];
sx q[0];
rz(-1.5769616) q[0];
sx q[0];
rz(0.023029285) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41202338) q[2];
sx q[2];
rz(-3.014866) q[2];
sx q[2];
rz(2.3840981) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16329855) q[1];
sx q[1];
rz(-2.1566609) q[1];
sx q[1];
rz(-2.3818124) q[1];
rz(-pi) q[2];
rz(-2.6074334) q[3];
sx q[3];
rz(-1.307752) q[3];
sx q[3];
rz(-2.0260889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0447023) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(-0.61075413) q[2];
rz(-2.6664873) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(0.92774123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.24191813) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(-2.8444667) q[0];
rz(1.7469453) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(0.64613211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.811759) q[0];
sx q[0];
rz(-2.9114897) q[0];
sx q[0];
rz(-2.0287201) q[0];
x q[1];
rz(1.4023151) q[2];
sx q[2];
rz(-1.104276) q[2];
sx q[2];
rz(0.70665765) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7568126) q[1];
sx q[1];
rz(-1.2287178) q[1];
sx q[1];
rz(2.6095819) q[1];
x q[2];
rz(-2.5313247) q[3];
sx q[3];
rz(-0.90859298) q[3];
sx q[3];
rz(-2.3074647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0769161) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(-0.70139766) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69944537) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(0.79750693) q[0];
rz(-2.6240255) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(-3.033175) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90005504) q[0];
sx q[0];
rz(-0.90796048) q[0];
sx q[0];
rz(0.073297757) q[0];
rz(-pi) q[1];
rz(-1.7494781) q[2];
sx q[2];
rz(-1.5249426) q[2];
sx q[2];
rz(-2.1669441) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0260967) q[1];
sx q[1];
rz(-2.8671088) q[1];
sx q[1];
rz(-2.2309169) q[1];
rz(-pi) q[2];
rz(-2.9547144) q[3];
sx q[3];
rz(-2.2664321) q[3];
sx q[3];
rz(-2.7587492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0124399) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-0.33995315) q[2];
rz(-2.7231976) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-2.4160014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5323935) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(0.6788196) q[0];
rz(-2.7774096) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(-3.0864339) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.508651) q[0];
sx q[0];
rz(-1.5061437) q[0];
sx q[0];
rz(-2.7642194) q[0];
rz(-pi) q[1];
rz(-1.360838) q[2];
sx q[2];
rz(-1.0210438) q[2];
sx q[2];
rz(-1.9090261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2436279) q[1];
sx q[1];
rz(-0.98476714) q[1];
sx q[1];
rz(-2.2378053) q[1];
rz(-pi) q[2];
rz(-1.6454562) q[3];
sx q[3];
rz(-0.45711043) q[3];
sx q[3];
rz(0.32170579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1577592) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(-2.6514163) q[2];
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
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4162083) q[0];
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
rz(-0.18110885) q[3];
sx q[3];
rz(-2.8977179) q[3];
sx q[3];
rz(0.44117622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
