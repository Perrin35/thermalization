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
rz(1.376027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3289514) q[0];
sx q[0];
rz(-2.0094123) q[0];
sx q[0];
rz(-2.2493275) q[0];
rz(-pi) q[1];
rz(0.55144989) q[2];
sx q[2];
rz(-0.79556888) q[2];
sx q[2];
rz(-1.7099107) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8810597) q[1];
sx q[1];
rz(-2.7785289) q[1];
sx q[1];
rz(-2.0102324) q[1];
rz(2.3222378) q[3];
sx q[3];
rz(-1.6091533) q[3];
sx q[3];
rz(-2.0824144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(-0.98891813) q[2];
rz(2.3890498) q[3];
sx q[3];
rz(-1.9957333) q[3];
sx q[3];
rz(2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20724021) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(-1.1799312) q[0];
rz(0.99769366) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(-2.4172799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1441919) q[0];
sx q[0];
rz(-1.8840944) q[0];
sx q[0];
rz(-2.2749167) q[0];
rz(-pi) q[1];
rz(-2.9421147) q[2];
sx q[2];
rz(-1.6530767) q[2];
sx q[2];
rz(-0.85180887) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7608632) q[1];
sx q[1];
rz(-2.7298096) q[1];
sx q[1];
rz(2.1058583) q[1];
rz(-pi) q[2];
rz(-2.4507387) q[3];
sx q[3];
rz(-2.8375531) q[3];
sx q[3];
rz(0.090304852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90536845) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(-2.779707) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996465) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(-1.746159) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(-1.2190855) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0116545) q[0];
sx q[0];
rz(-1.9217976) q[0];
sx q[0];
rz(0.28052335) q[0];
x q[1];
rz(-1.9169541) q[2];
sx q[2];
rz(-0.73181728) q[2];
sx q[2];
rz(-3.0423622) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0137274) q[1];
sx q[1];
rz(-1.3294819) q[1];
sx q[1];
rz(-2.1167386) q[1];
x q[2];
rz(1.838802) q[3];
sx q[3];
rz(-1.5558814) q[3];
sx q[3];
rz(-0.55004317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7200155) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(-1.4455618) q[2];
rz(0.56882632) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22359426) q[0];
sx q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-2.6233327) q[0];
rz(-2.4261684) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(2.3148361) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2430192) q[0];
sx q[0];
rz(-1.9802226) q[0];
sx q[0];
rz(-2.5893674) q[0];
rz(1.7051758) q[2];
sx q[2];
rz(-0.80766404) q[2];
sx q[2];
rz(-2.813051) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.066847853) q[1];
sx q[1];
rz(-1.8555992) q[1];
sx q[1];
rz(-0.39798255) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7923146) q[3];
sx q[3];
rz(-0.75540245) q[3];
sx q[3];
rz(-2.3625771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9892019) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(1.3789122) q[2];
rz(-0.072323024) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(-1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3146661) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(-1.1258874) q[0];
rz(0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(-2.856423) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0774768) q[0];
sx q[0];
rz(-0.47668326) q[0];
sx q[0];
rz(1.6633196) q[0];
rz(2.4777849) q[2];
sx q[2];
rz(-1.8823349) q[2];
sx q[2];
rz(0.89154348) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.8385182) q[1];
sx q[1];
rz(-2.8285366) q[1];
sx q[1];
rz(-1.5721553) q[1];
x q[2];
rz(-0.40163715) q[3];
sx q[3];
rz(-1.4816227) q[3];
sx q[3];
rz(3.0107486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.29331648) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(0.53945333) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(-0.45421281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25813112) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(-2.4482588) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(-1.3670115) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54129823) q[0];
sx q[0];
rz(-1.5026662) q[0];
sx q[0];
rz(-3.1085204) q[0];
rz(2.1742646) q[2];
sx q[2];
rz(-2.2852995) q[2];
sx q[2];
rz(-0.29010233) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.94783084) q[1];
sx q[1];
rz(-1.8596706) q[1];
sx q[1];
rz(3.0572901) q[1];
rz(2.7886224) q[3];
sx q[3];
rz(-1.6680696) q[3];
sx q[3];
rz(0.73247611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(2.7381251) q[2];
rz(2.6599595) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(-2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-0.8738628) q[0];
rz(-0.44772398) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(-1.1706932) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8255575) q[0];
sx q[0];
rz(-3.1177525) q[0];
sx q[0];
rz(0.26160474) q[0];
x q[1];
rz(2.7295693) q[2];
sx q[2];
rz(-0.12672666) q[2];
sx q[2];
rz(2.3840981) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.16329855) q[1];
sx q[1];
rz(-2.1566609) q[1];
sx q[1];
rz(-2.3818124) q[1];
rz(2.6550754) q[3];
sx q[3];
rz(-0.58972893) q[3];
sx q[3];
rz(-3.1004578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0447023) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(-0.61075413) q[2];
rz(-0.47510535) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8996745) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(0.29712594) q[0];
rz(-1.7469453) q[1];
sx q[1];
rz(-1.1509832) q[1];
sx q[1];
rz(0.64613211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.811759) q[0];
sx q[0];
rz(-0.23010294) q[0];
sx q[0];
rz(-1.1128725) q[0];
x q[1];
rz(0.32142873) q[2];
sx q[2];
rz(-2.6476963) q[2];
sx q[2];
rz(2.0733881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.38478002) q[1];
sx q[1];
rz(-1.2287178) q[1];
sx q[1];
rz(-2.6095819) q[1];
rz(-0.93692245) q[3];
sx q[3];
rz(-0.86808944) q[3];
sx q[3];
rz(0.015451775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.064676553) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(-0.45483744) q[2];
rz(0.70139766) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(-2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69944537) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(2.3440857) q[0];
rz(-0.51756716) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(0.10841766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3603044) q[0];
sx q[0];
rz(-2.4753248) q[0];
sx q[0];
rz(1.4772619) q[0];
x q[1];
rz(-0.046594521) q[2];
sx q[2];
rz(-1.7492883) q[2];
sx q[2];
rz(2.5537234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.91067578) q[1];
rz(-pi) q[2];
rz(1.3518203) q[3];
sx q[3];
rz(-2.4253546) q[3];
sx q[3];
rz(0.66974528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-0.33995315) q[2];
rz(0.41839504) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-2.4160014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5323935) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(-2.4627731) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(-3.0864339) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.508651) q[0];
sx q[0];
rz(-1.635449) q[0];
sx q[0];
rz(2.7642194) q[0];
x q[1];
rz(-0.55969413) q[2];
sx q[2];
rz(-1.392138) q[2];
sx q[2];
rz(2.6924804) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0835417) q[1];
sx q[1];
rz(-1.0293048) q[1];
sx q[1];
rz(2.4400649) q[1];
rz(-1.1147898) q[3];
sx q[3];
rz(-1.6037233) q[3];
sx q[3];
rz(1.959521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1577592) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(2.6514163) q[2];
rz(0.13752078) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(0.93808758) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4162083) q[0];
sx q[0];
rz(-1.890247) q[0];
sx q[0];
rz(2.4631137) q[0];
rz(-2.9329119) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(2.8293777) q[2];
sx q[2];
rz(-1.6630465) q[2];
sx q[2];
rz(1.9065471) q[2];
rz(-2.9604838) q[3];
sx q[3];
rz(-0.24387471) q[3];
sx q[3];
rz(-2.7004164) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
