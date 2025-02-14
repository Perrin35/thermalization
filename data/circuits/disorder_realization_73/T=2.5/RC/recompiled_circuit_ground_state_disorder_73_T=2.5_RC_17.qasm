OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8483228) q[0];
sx q[0];
rz(-0.25265101) q[0];
sx q[0];
rz(1.9360315) q[0];
rz(-1.128101) q[1];
sx q[1];
rz(4.4681273) q[1];
sx q[1];
rz(11.820643) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6572031) q[0];
sx q[0];
rz(-1.5108539) q[0];
sx q[0];
rz(-1.3482698) q[0];
rz(-1.5640075) q[2];
sx q[2];
rz(-1.7951843) q[2];
sx q[2];
rz(0.16358384) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.89167833) q[1];
sx q[1];
rz(-0.39977705) q[1];
sx q[1];
rz(-1.0067902) q[1];
rz(1.1806335) q[3];
sx q[3];
rz(-0.50658617) q[3];
sx q[3];
rz(-0.49662874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2241609) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(-1.6708299) q[2];
rz(1.5077695) q[3];
sx q[3];
rz(-0.016851146) q[3];
sx q[3];
rz(-0.9233709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.5725752) q[0];
sx q[0];
rz(-1.1976396) q[0];
sx q[0];
rz(1.5684599) q[0];
rz(-0.16954999) q[1];
sx q[1];
rz(-0.1145656) q[1];
sx q[1];
rz(3.0068908) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41551521) q[0];
sx q[0];
rz(-1.2252136) q[0];
sx q[0];
rz(-2.6879361) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7472166) q[2];
sx q[2];
rz(-0.069709502) q[2];
sx q[2];
rz(1.4331872) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33011064) q[1];
sx q[1];
rz(-2.5977784) q[1];
sx q[1];
rz(2.5303929) q[1];
rz(-pi) q[2];
rz(-0.080218519) q[3];
sx q[3];
rz(-1.748198) q[3];
sx q[3];
rz(2.3631848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0160825) q[2];
sx q[2];
rz(-1.6566015) q[2];
sx q[2];
rz(2.9888195) q[2];
rz(-1.3656535) q[3];
sx q[3];
rz(-0.036866166) q[3];
sx q[3];
rz(-0.16807817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.633054) q[0];
sx q[0];
rz(-2.3336053) q[0];
sx q[0];
rz(0.48164865) q[0];
rz(-0.18474361) q[1];
sx q[1];
rz(-1.7748723) q[1];
sx q[1];
rz(0.99536037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99670519) q[0];
sx q[0];
rz(-1.2249682) q[0];
sx q[0];
rz(-0.24390999) q[0];
rz(-pi) q[1];
rz(-3.0929933) q[2];
sx q[2];
rz(-1.6311797) q[2];
sx q[2];
rz(0.27602613) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9650594) q[1];
sx q[1];
rz(-0.73672026) q[1];
sx q[1];
rz(-2.6882369) q[1];
rz(-0.030582436) q[3];
sx q[3];
rz(-0.97506071) q[3];
sx q[3];
rz(-3.093978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2153726) q[2];
sx q[2];
rz(-0.054457713) q[2];
sx q[2];
rz(-0.064662956) q[2];
rz(-1.0489382) q[3];
sx q[3];
rz(-0.026853042) q[3];
sx q[3];
rz(1.2803199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30508405) q[0];
sx q[0];
rz(-0.11798141) q[0];
sx q[0];
rz(0.82103658) q[0];
rz(0.07846421) q[1];
sx q[1];
rz(-1.6946225) q[1];
sx q[1];
rz(-2.1441114) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2729014) q[0];
sx q[0];
rz(-1.8857546) q[0];
sx q[0];
rz(-0.64906831) q[0];
rz(-2.1199662) q[2];
sx q[2];
rz(-0.06937521) q[2];
sx q[2];
rz(0.79909426) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8197937) q[1];
sx q[1];
rz(-0.60056486) q[1];
sx q[1];
rz(1.0933881) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0313052) q[3];
sx q[3];
rz(-0.48088851) q[3];
sx q[3];
rz(0.65910027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0923826) q[2];
sx q[2];
rz(-3.1123078) q[2];
sx q[2];
rz(1.1537665) q[2];
rz(0.18618259) q[3];
sx q[3];
rz(-3.0598873) q[3];
sx q[3];
rz(-2.8103099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1368197) q[0];
sx q[0];
rz(-0.81757075) q[0];
sx q[0];
rz(-1.1432884) q[0];
rz(1.1413057) q[1];
sx q[1];
rz(-2.3316796) q[1];
sx q[1];
rz(-0.51796651) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.98949) q[0];
sx q[0];
rz(-0.94606646) q[0];
sx q[0];
rz(-2.2476907) q[0];
rz(-pi) q[1];
rz(-1.5094366) q[2];
sx q[2];
rz(-0.01505919) q[2];
sx q[2];
rz(1.3850152) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3822381) q[1];
sx q[1];
rz(-1.2975946) q[1];
sx q[1];
rz(1.9135936) q[1];
x q[2];
rz(-0.2980026) q[3];
sx q[3];
rz(-1.7116705) q[3];
sx q[3];
rz(-2.6948787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3190069) q[2];
sx q[2];
rz(-2.1876882) q[2];
sx q[2];
rz(0.32773584) q[2];
rz(1.94708) q[3];
sx q[3];
rz(-0.12735282) q[3];
sx q[3];
rz(-2.2849042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.0026534) q[0];
sx q[0];
rz(-2.8728573) q[0];
sx q[0];
rz(-0.56185454) q[0];
rz(-1.4909164) q[1];
sx q[1];
rz(-1.5166538) q[1];
sx q[1];
rz(3.0432826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54360169) q[0];
sx q[0];
rz(-2.7181135) q[0];
sx q[0];
rz(-1.4844358) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19616429) q[2];
sx q[2];
rz(-0.00065302424) q[2];
sx q[2];
rz(0.12015039) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5206388) q[1];
sx q[1];
rz(-2.6700182) q[1];
sx q[1];
rz(0.10271272) q[1];
rz(-1.1071854) q[3];
sx q[3];
rz(-2.0258198) q[3];
sx q[3];
rz(0.95038271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.085122846) q[2];
sx q[2];
rz(-0.19530185) q[2];
sx q[2];
rz(-1.0657715) q[2];
rz(-2.7417475) q[3];
sx q[3];
rz(-0.53374922) q[3];
sx q[3];
rz(1.8736418) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14168508) q[0];
sx q[0];
rz(-0.081900224) q[0];
sx q[0];
rz(1.4578693) q[0];
rz(-2.0211925) q[1];
sx q[1];
rz(-0.1381865) q[1];
sx q[1];
rz(2.8021326) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0507815) q[0];
sx q[0];
rz(-0.078074038) q[0];
sx q[0];
rz(-0.17103057) q[0];
rz(-pi) q[1];
rz(2.5869114) q[2];
sx q[2];
rz(-1.5590057) q[2];
sx q[2];
rz(1.5654711) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.554018) q[1];
sx q[1];
rz(-1.4851991) q[1];
sx q[1];
rz(-0.0016335131) q[1];
rz(-pi) q[2];
rz(-1.4711597) q[3];
sx q[3];
rz(-1.2335586) q[3];
sx q[3];
rz(1.3226313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7772943) q[2];
sx q[2];
rz(-0.0019145049) q[2];
sx q[2];
rz(-0.36336362) q[2];
rz(1.0904788) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(-2.0232078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61984396) q[0];
sx q[0];
rz(-0.32829568) q[0];
sx q[0];
rz(-0.92292619) q[0];
rz(-1.6587616) q[1];
sx q[1];
rz(-2.5117579) q[1];
sx q[1];
rz(-3.1373851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.370458) q[0];
sx q[0];
rz(-0.97386375) q[0];
sx q[0];
rz(2.202522) q[0];
x q[1];
rz(0.0059295456) q[2];
sx q[2];
rz(-1.9816795) q[2];
sx q[2];
rz(-3.1290999) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9862822) q[1];
sx q[1];
rz(-0.48008568) q[1];
sx q[1];
rz(1.4309149) q[1];
x q[2];
rz(-1.6170623) q[3];
sx q[3];
rz(-0.95887254) q[3];
sx q[3];
rz(0.16976742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5620455) q[2];
sx q[2];
rz(-1.6041218) q[2];
sx q[2];
rz(-1.2008249) q[2];
rz(-1.7480525) q[3];
sx q[3];
rz(-3.1378919) q[3];
sx q[3];
rz(0.70011955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8414108) q[0];
sx q[0];
rz(-0.44898471) q[0];
sx q[0];
rz(-1.2196983) q[0];
rz(-1.3297184) q[1];
sx q[1];
rz(-2.0025573) q[1];
sx q[1];
rz(-0.16389287) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4224098) q[0];
sx q[0];
rz(-1.9253732) q[0];
sx q[0];
rz(2.9585881) q[0];
rz(-pi) q[1];
rz(0.58456011) q[2];
sx q[2];
rz(-1.341408) q[2];
sx q[2];
rz(3.0363415) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5245318) q[1];
sx q[1];
rz(-1.5111898) q[1];
sx q[1];
rz(1.6888794) q[1];
rz(-pi) q[2];
rz(-0.71299841) q[3];
sx q[3];
rz(-0.87942356) q[3];
sx q[3];
rz(2.612243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4280052) q[2];
sx q[2];
rz(-0.089944936) q[2];
sx q[2];
rz(-0.62291992) q[2];
rz(-0.04341393) q[3];
sx q[3];
rz(-0.89689887) q[3];
sx q[3];
rz(-0.77891946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49122214) q[0];
sx q[0];
rz(-0.0064042052) q[0];
sx q[0];
rz(-2.6553335) q[0];
rz(-0.69475118) q[1];
sx q[1];
rz(-2.7943352) q[1];
sx q[1];
rz(2.6659226) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3846426) q[0];
sx q[0];
rz(-1.5371377) q[0];
sx q[0];
rz(1.606694) q[0];
rz(3.1074355) q[2];
sx q[2];
rz(-2.2766621) q[2];
sx q[2];
rz(-1.5868843) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.025534677) q[1];
sx q[1];
rz(-0.60760159) q[1];
sx q[1];
rz(-0.011354253) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6127897) q[3];
sx q[3];
rz(-0.35548726) q[3];
sx q[3];
rz(-0.58099174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4645369) q[2];
sx q[2];
rz(-0.060364351) q[2];
sx q[2];
rz(1.8439058) q[2];
rz(-1.248598) q[3];
sx q[3];
rz(-0.55515754) q[3];
sx q[3];
rz(-0.36778522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1260592) q[0];
sx q[0];
rz(-1.561469) q[0];
sx q[0];
rz(-1.4373686) q[0];
rz(2.2592648) q[1];
sx q[1];
rz(-0.066451646) q[1];
sx q[1];
rz(0.8393504) q[1];
rz(2.7693314) q[2];
sx q[2];
rz(-0.098975565) q[2];
sx q[2];
rz(3.0437058) q[2];
rz(-0.23918693) q[3];
sx q[3];
rz(-1.5659955) q[3];
sx q[3];
rz(-1.5495054) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
