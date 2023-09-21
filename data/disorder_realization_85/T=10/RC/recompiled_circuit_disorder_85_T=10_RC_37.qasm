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
rz(0.37101775) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(2.5420904) q[1];
sx q[1];
rz(11.190344) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5698619) q[0];
sx q[0];
rz(-0.96643448) q[0];
sx q[0];
rz(-0.5423003) q[0];
rz(-pi) q[1];
rz(-2.5901428) q[2];
sx q[2];
rz(-0.79556888) q[2];
sx q[2];
rz(1.431682) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.260533) q[1];
sx q[1];
rz(-0.36306371) q[1];
sx q[1];
rz(1.1313603) q[1];
rz(1.6269496) q[3];
sx q[3];
rz(-0.75222844) q[3];
sx q[3];
rz(-2.5889531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0573037) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(-2.1526745) q[2];
rz(-2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20724021) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(1.9616615) q[0];
rz(0.99769366) q[1];
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
rz(2.9162826) q[0];
sx q[0];
rz(-2.3819469) q[0];
sx q[0];
rz(1.1067953) q[0];
x q[1];
rz(-1.486859) q[2];
sx q[2];
rz(-1.7695904) q[2];
sx q[2];
rz(-0.73560152) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.69246768) q[1];
sx q[1];
rz(-1.3652703) q[1];
sx q[1];
rz(-1.9301901) q[1];
x q[2];
rz(-1.3734666) q[3];
sx q[3];
rz(-1.337968) q[3];
sx q[3];
rz(0.80430921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90536845) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(2.779707) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
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
rz(0.46229258) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(-1.2190855) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65788668) q[0];
sx q[0];
rz(-1.3077967) q[0];
sx q[0];
rz(1.9348683) q[0];
rz(0.29580446) q[2];
sx q[2];
rz(-0.89106262) q[2];
sx q[2];
rz(-2.5909397) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3238941) q[1];
sx q[1];
rz(-2.5496799) q[1];
sx q[1];
rz(-2.0134258) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6270646) q[3];
sx q[3];
rz(-0.26841044) q[3];
sx q[3];
rz(-1.0750107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7200155) q[2];
sx q[2];
rz(-2.4013459) q[2];
sx q[2];
rz(-1.4455618) q[2];
rz(-0.56882632) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(-3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9179984) q[0];
sx q[0];
rz(-2.693394) q[0];
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
rz(-2.2402156) q[0];
sx q[0];
rz(-2.4670521) q[0];
sx q[0];
rz(2.4504689) q[0];
x q[1];
rz(0.13917285) q[2];
sx q[2];
rz(-2.3690802) q[2];
sx q[2];
rz(0.13538361) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0747448) q[1];
sx q[1];
rz(-1.2859935) q[1];
sx q[1];
rz(0.39798255) q[1];
x q[2];
rz(-0.72439648) q[3];
sx q[3];
rz(-1.3339692) q[3];
sx q[3];
rz(-0.5326007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9892019) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(3.0692696) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(-1.4962083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82692659) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(-1.1258874) q[0];
rz(-0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(-0.28516969) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7171705) q[0];
sx q[0];
rz(-1.5283913) q[0];
sx q[0];
rz(-2.0457343) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9589013) q[2];
sx q[2];
rz(-0.94411196) q[2];
sx q[2];
rz(2.2270122) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3030745) q[1];
sx q[1];
rz(-2.8285366) q[1];
sx q[1];
rz(-1.5721553) q[1];
rz(-pi) q[2];
rz(-0.40163715) q[3];
sx q[3];
rz(-1.4816227) q[3];
sx q[3];
rz(3.0107486) q[3];
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
rz(-0.30682492) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(-2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8834615) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(0.15701292) q[0];
rz(-2.4482588) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(-1.7745811) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088928662) q[0];
sx q[0];
rz(-0.0757218) q[0];
sx q[0];
rz(2.0220387) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3301666) q[2];
sx q[2];
rz(-1.1277414) q[2];
sx q[2];
rz(1.4366988) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94783084) q[1];
sx q[1];
rz(-1.281922) q[1];
sx q[1];
rz(3.0572901) q[1];
x q[2];
rz(-2.8664687) q[3];
sx q[3];
rz(-0.36558662) q[3];
sx q[3];
rz(2.5610353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.29629016) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(-2.7381251) q[2];
rz(2.6599595) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(0.51923716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24213174) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(0.8738628) q[0];
rz(-2.6938687) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(1.1706932) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51629492) q[0];
sx q[0];
rz(-1.5769616) q[0];
sx q[0];
rz(0.023029285) q[0];
rz(-pi) q[1];
rz(1.5198176) q[2];
sx q[2];
rz(-1.6868674) q[2];
sx q[2];
rz(2.799084) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9782941) q[1];
sx q[1];
rz(-2.1566609) q[1];
sx q[1];
rz(0.75978029) q[1];
x q[2];
rz(2.6074334) q[3];
sx q[3];
rz(-1.8338406) q[3];
sx q[3];
rz(-2.0260889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0968904) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(2.5308385) q[2];
rz(0.47510535) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24191813) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(-0.29712594) q[0];
rz(1.3946474) q[1];
sx q[1];
rz(-1.1509832) q[1];
sx q[1];
rz(0.64613211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6883811) q[0];
sx q[0];
rz(-1.6717981) q[0];
sx q[0];
rz(-1.3636916) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32142873) q[2];
sx q[2];
rz(-2.6476963) q[2];
sx q[2];
rz(1.0682046) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38478002) q[1];
sx q[1];
rz(-1.9128748) q[1];
sx q[1];
rz(-2.6095819) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2046702) q[3];
sx q[3];
rz(-2.2735032) q[3];
sx q[3];
rz(0.015451775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.064676553) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(0.70139766) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(1.1340244) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4421473) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(-2.3440857) q[0];
rz(-2.6240255) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(3.033175) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.425697) q[0];
sx q[0];
rz(-1.6285537) q[0];
sx q[0];
rz(-2.2349368) q[0];
x q[1];
rz(-1.7494781) q[2];
sx q[2];
rz(-1.6166501) q[2];
sx q[2];
rz(-0.97464857) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.79418102) q[1];
sx q[1];
rz(-1.3550183) q[1];
sx q[1];
rz(-2.9706035) q[1];
x q[2];
rz(1.7897723) q[3];
sx q[3];
rz(-2.4253546) q[3];
sx q[3];
rz(-0.66974528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(0.33995315) q[2];
rz(2.7231976) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(2.4160014) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5323935) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(-0.6788196) q[0];
rz(-2.7774096) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(-3.0864339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63294166) q[0];
sx q[0];
rz(-1.635449) q[0];
sx q[0];
rz(-2.7642194) q[0];
rz(-pi) q[1];
rz(1.360838) q[2];
sx q[2];
rz(-1.0210438) q[2];
sx q[2];
rz(-1.2325665) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2436279) q[1];
sx q[1];
rz(-0.98476714) q[1];
sx q[1];
rz(-0.90378739) q[1];
rz(-pi) q[2];
rz(3.1049214) q[3];
sx q[3];
rz(-1.1150556) q[3];
sx q[3];
rz(-0.40487056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72538439) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(0.2086808) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(-2.8490534) q[2];
sx q[2];
rz(-0.32513041) q[2];
sx q[2];
rz(0.61376094) q[2];
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
