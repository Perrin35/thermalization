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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8810597) q[1];
sx q[1];
rz(-0.36306371) q[1];
sx q[1];
rz(-1.1313603) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81935482) q[3];
sx q[3];
rz(-1.5324394) q[3];
sx q[3];
rz(-2.0824144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0573037) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(-2.1526745) q[2];
rz(2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(-2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.20724021) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(1.9616615) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(-0.72431272) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2253101) q[0];
sx q[0];
rz(-0.75964576) q[0];
sx q[0];
rz(-2.0347974) q[0];
rz(1.6547336) q[2];
sx q[2];
rz(-1.7695904) q[2];
sx q[2];
rz(2.4059911) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.69246768) q[1];
sx q[1];
rz(-1.3652703) q[1];
sx q[1];
rz(-1.9301901) q[1];
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
rz(2.2362242) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(2.779707) q[2];
rz(0.13606717) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(-3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.1419462) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(1.3954337) q[0];
rz(-0.46229258) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(-1.9225072) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31419471) q[0];
sx q[0];
rz(-0.44566804) q[0];
sx q[0];
rz(-2.218194) q[0];
rz(-2.2723324) q[2];
sx q[2];
rz(-1.7995036) q[2];
sx q[2];
rz(1.2094091) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8176986) q[1];
sx q[1];
rz(-0.59191275) q[1];
sx q[1];
rz(-1.1281668) q[1];
rz(-pi) q[2];
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
rz(-2.4013459) q[2];
sx q[2];
rz(-1.6960309) q[2];
rz(-2.5727663) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(-0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22359426) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(2.6233327) q[0];
rz(-0.7154243) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(0.82675654) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90137705) q[0];
sx q[0];
rz(-2.4670521) q[0];
sx q[0];
rz(0.69112372) q[0];
rz(-pi) q[1];
rz(-0.76765676) q[2];
sx q[2];
rz(-1.6677688) q[2];
sx q[2];
rz(-1.8061639) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0747448) q[1];
sx q[1];
rz(-1.8555992) q[1];
sx q[1];
rz(2.7436101) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8825674) q[3];
sx q[3];
rz(-0.87083737) q[3];
sx q[3];
rz(1.8986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.15239079) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(1.3789122) q[2];
rz(0.072323024) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(1.6453843) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82692659) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(1.1258874) q[0];
rz(0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(-2.856423) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0641159) q[0];
sx q[0];
rz(-0.47668326) q[0];
sx q[0];
rz(-1.6633196) q[0];
x q[1];
rz(-1.1826913) q[2];
sx q[2];
rz(-0.94411196) q[2];
sx q[2];
rz(0.91458048) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3030745) q[1];
sx q[1];
rz(-0.31305602) q[1];
sx q[1];
rz(1.5721553) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4739591) q[3];
sx q[3];
rz(-1.9707465) q[3];
sx q[3];
rz(1.4777583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8482762) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(0.53945333) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(0.45421281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.8834615) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-0.15701292) q[0];
rz(-2.4482588) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(-1.7745811) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.052664) q[0];
sx q[0];
rz(-0.0757218) q[0];
sx q[0];
rz(-2.0220387) q[0];
rz(-pi) q[1];
rz(-0.96732803) q[2];
sx q[2];
rz(-2.2852995) q[2];
sx q[2];
rz(2.8514903) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4945592) q[1];
sx q[1];
rz(-1.651598) q[1];
sx q[1];
rz(-1.8606436) q[1];
x q[2];
rz(-1.6744162) q[3];
sx q[3];
rz(-1.2195671) q[3];
sx q[3];
rz(2.2675089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(0.40346754) q[2];
rz(-0.48163313) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(-2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24213174) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(2.2677299) q[0];
rz(0.44772398) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(-1.9708995) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0543594) q[0];
sx q[0];
rz(-1.5938252) q[0];
sx q[0];
rz(-1.5646294) q[0];
rz(-pi) q[1];
rz(2.7295693) q[2];
sx q[2];
rz(-0.12672666) q[2];
sx q[2];
rz(-0.7574946) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9782941) q[1];
sx q[1];
rz(-0.98493176) q[1];
sx q[1];
rz(-0.75978029) q[1];
x q[2];
rz(2.6074334) q[3];
sx q[3];
rz(-1.307752) q[3];
sx q[3];
rz(2.0260889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0447023) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(-2.5308385) q[2];
rz(2.6664873) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(-2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24191813) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(0.29712594) q[0];
rz(1.7469453) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(0.64613211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.811759) q[0];
sx q[0];
rz(-0.23010294) q[0];
sx q[0];
rz(-1.1128725) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32142873) q[2];
sx q[2];
rz(-0.49389631) q[2];
sx q[2];
rz(-1.0682046) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1505193) q[1];
sx q[1];
rz(-1.0725613) q[1];
sx q[1];
rz(-1.9626161) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2046702) q[3];
sx q[3];
rz(-2.2735032) q[3];
sx q[3];
rz(3.1261409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0769161) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(-0.45483744) q[2];
rz(2.440195) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2415376) q[0];
sx q[0];
rz(-0.90796048) q[0];
sx q[0];
rz(0.073297757) q[0];
rz(1.3921146) q[2];
sx q[2];
rz(-1.6166501) q[2];
sx q[2];
rz(2.1669441) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.79418102) q[1];
sx q[1];
rz(-1.3550183) q[1];
sx q[1];
rz(0.17098917) q[1];
rz(-1.3518203) q[3];
sx q[3];
rz(-0.7162381) q[3];
sx q[3];
rz(-2.4718474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1291528) q[2];
sx q[2];
rz(-1.3468578) q[2];
sx q[2];
rz(2.8016395) q[2];
rz(-2.7231976) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(-0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6091992) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(-2.4627731) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(3.0864339) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.178135) q[0];
sx q[0];
rz(-1.9473416) q[0];
sx q[0];
rz(1.5012653) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55969413) q[2];
sx q[2];
rz(-1.392138) q[2];
sx q[2];
rz(-2.6924804) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.89796472) q[1];
sx q[1];
rz(-0.98476714) q[1];
sx q[1];
rz(-0.90378739) q[1];
rz(3.1049214) q[3];
sx q[3];
rz(-1.1150556) q[3];
sx q[3];
rz(2.7367221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1577592) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(0.49017635) q[2];
rz(3.0040719) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(-2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(pi/2) q[0];
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
rz(-1.6155852) q[3];
sx q[3];
rz(-1.3309892) q[3];
sx q[3];
rz(0.25467024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
