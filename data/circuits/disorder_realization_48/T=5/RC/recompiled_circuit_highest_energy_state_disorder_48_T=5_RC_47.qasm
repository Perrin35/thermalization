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
rz(-2.2263865) q[0];
sx q[0];
rz(-2.6829166) q[0];
sx q[0];
rz(-0.31804481) q[0];
rz(-1.7755427) q[1];
sx q[1];
rz(-0.69861689) q[1];
sx q[1];
rz(3.0787992) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013115766) q[0];
sx q[0];
rz(-2.4139691) q[0];
sx q[0];
rz(-2.3250871) q[0];
rz(0.89102052) q[2];
sx q[2];
rz(-1.2212379) q[2];
sx q[2];
rz(1.0206136) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.56927189) q[1];
sx q[1];
rz(-0.57479984) q[1];
sx q[1];
rz(2.1009836) q[1];
x q[2];
rz(2.1413689) q[3];
sx q[3];
rz(-0.65815845) q[3];
sx q[3];
rz(-2.6702945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41844765) q[2];
sx q[2];
rz(-1.8258839) q[2];
sx q[2];
rz(-1.0203699) q[2];
rz(-2.4127035) q[3];
sx q[3];
rz(-0.43292361) q[3];
sx q[3];
rz(-1.5322022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51434022) q[0];
sx q[0];
rz(-1.9753375) q[0];
sx q[0];
rz(3.1033206) q[0];
rz(-3.017784) q[1];
sx q[1];
rz(-2.6688711) q[1];
sx q[1];
rz(1.5708057) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0012624) q[0];
sx q[0];
rz(-1.5902983) q[0];
sx q[0];
rz(0.01000761) q[0];
x q[1];
rz(-3.0395987) q[2];
sx q[2];
rz(-2.3677459) q[2];
sx q[2];
rz(0.71600658) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5396186) q[1];
sx q[1];
rz(-1.5718196) q[1];
sx q[1];
rz(-3.1411693) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9306049) q[3];
sx q[3];
rz(-2.6158276) q[3];
sx q[3];
rz(2.7069382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.52246419) q[2];
sx q[2];
rz(-1.3294486) q[2];
sx q[2];
rz(2.5627356) q[2];
rz(-2.0959496) q[3];
sx q[3];
rz(-2.3561616) q[3];
sx q[3];
rz(2.3845909) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66913644) q[0];
sx q[0];
rz(-2.0643015) q[0];
sx q[0];
rz(-0.45397595) q[0];
rz(0.16481915) q[1];
sx q[1];
rz(-0.77607981) q[1];
sx q[1];
rz(-1.4705315) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9605345) q[0];
sx q[0];
rz(-1.1832602) q[0];
sx q[0];
rz(-2.9563532) q[0];
rz(2.8118089) q[2];
sx q[2];
rz(-2.0726786) q[2];
sx q[2];
rz(1.9374073) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78782234) q[1];
sx q[1];
rz(-2.2332237) q[1];
sx q[1];
rz(-1.9595998) q[1];
x q[2];
rz(-0.65916797) q[3];
sx q[3];
rz(-1.1416832) q[3];
sx q[3];
rz(1.8076713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37956023) q[2];
sx q[2];
rz(-1.4780059) q[2];
sx q[2];
rz(1.7851768) q[2];
rz(-2.6421269) q[3];
sx q[3];
rz(-1.6127337) q[3];
sx q[3];
rz(-1.0511901) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9143963) q[0];
sx q[0];
rz(-2.2375186) q[0];
sx q[0];
rz(-2.0830578) q[0];
rz(-1.1610441) q[1];
sx q[1];
rz(-2.5807022) q[1];
sx q[1];
rz(3.0243691) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6446021) q[0];
sx q[0];
rz(-1.0287026) q[0];
sx q[0];
rz(1.5250456) q[0];
rz(1.8294705) q[2];
sx q[2];
rz(-0.66159464) q[2];
sx q[2];
rz(-1.0115397) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8414054) q[1];
sx q[1];
rz(-0.8148199) q[1];
sx q[1];
rz(0.38736613) q[1];
rz(-0.63734148) q[3];
sx q[3];
rz(-1.2351324) q[3];
sx q[3];
rz(-2.6791399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.78493541) q[2];
sx q[2];
rz(-1.4484118) q[2];
sx q[2];
rz(1.3680722) q[2];
rz(-1.8562227) q[3];
sx q[3];
rz(-1.1077489) q[3];
sx q[3];
rz(0.72803298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0990937) q[0];
sx q[0];
rz(-2.2137401) q[0];
sx q[0];
rz(-1.7294783) q[0];
rz(-1.0026576) q[1];
sx q[1];
rz(-0.24191562) q[1];
sx q[1];
rz(-1.5826506) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60022814) q[0];
sx q[0];
rz(-1.1523655) q[0];
sx q[0];
rz(0.61516841) q[0];
x q[1];
rz(-0.71835235) q[2];
sx q[2];
rz(-1.6334772) q[2];
sx q[2];
rz(-2.2557392) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.69935) q[1];
sx q[1];
rz(-2.1210056) q[1];
sx q[1];
rz(-0.027632874) q[1];
x q[2];
rz(1.2941957) q[3];
sx q[3];
rz(-2.8081354) q[3];
sx q[3];
rz(-1.5302304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5530508) q[2];
sx q[2];
rz(-1.7091227) q[2];
sx q[2];
rz(2.7480965) q[2];
rz(-0.95997512) q[3];
sx q[3];
rz(-2.4656651) q[3];
sx q[3];
rz(-2.1737607) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2734964) q[0];
sx q[0];
rz(-0.12729004) q[0];
sx q[0];
rz(-0.18336503) q[0];
rz(2.6496437) q[1];
sx q[1];
rz(-2.349647) q[1];
sx q[1];
rz(1.2194182) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8293835) q[0];
sx q[0];
rz(-1.4903194) q[0];
sx q[0];
rz(3.0521554) q[0];
rz(-1.4923138) q[2];
sx q[2];
rz(-1.6011136) q[2];
sx q[2];
rz(-0.91532367) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.26039429) q[1];
sx q[1];
rz(-1.9250084) q[1];
sx q[1];
rz(0.43448914) q[1];
x q[2];
rz(1.0991092) q[3];
sx q[3];
rz(-1.4955284) q[3];
sx q[3];
rz(-1.690563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.1313974) q[2];
sx q[2];
rz(-1.4364028) q[2];
sx q[2];
rz(0.08237002) q[2];
rz(0.92665893) q[3];
sx q[3];
rz(-0.72946531) q[3];
sx q[3];
rz(-1.5026708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53967404) q[0];
sx q[0];
rz(-0.32783666) q[0];
sx q[0];
rz(-2.4251921) q[0];
rz(-2.7377103) q[1];
sx q[1];
rz(-1.5733746) q[1];
sx q[1];
rz(1.4366038) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2879031) q[0];
sx q[0];
rz(-2.7349964) q[0];
sx q[0];
rz(1.9415744) q[0];
rz(-pi) q[1];
rz(1.7030969) q[2];
sx q[2];
rz(-2.0952333) q[2];
sx q[2];
rz(-2.14545) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.77715141) q[1];
sx q[1];
rz(-1.4713411) q[1];
sx q[1];
rz(-2.742139) q[1];
rz(0.72258805) q[3];
sx q[3];
rz(-2.1919985) q[3];
sx q[3];
rz(-0.91739839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7776103) q[2];
sx q[2];
rz(-0.97475514) q[2];
sx q[2];
rz(0.065356143) q[2];
rz(0.86723793) q[3];
sx q[3];
rz(-1.6403653) q[3];
sx q[3];
rz(1.8170099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6548178) q[0];
sx q[0];
rz(-1.7205394) q[0];
sx q[0];
rz(2.4105893) q[0];
rz(2.3664318) q[1];
sx q[1];
rz(-0.33439264) q[1];
sx q[1];
rz(-2.7269272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0867942) q[0];
sx q[0];
rz(-2.2761184) q[0];
sx q[0];
rz(-0.84673015) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.783466) q[2];
sx q[2];
rz(-1.050282) q[2];
sx q[2];
rz(1.6047603) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6952299) q[1];
sx q[1];
rz(-1.2712423) q[1];
sx q[1];
rz(-1.3079554) q[1];
rz(-0.88736262) q[3];
sx q[3];
rz(-2.931224) q[3];
sx q[3];
rz(2.8223245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23184648) q[2];
sx q[2];
rz(-0.48603386) q[2];
sx q[2];
rz(0.092378423) q[2];
rz(0.77629027) q[3];
sx q[3];
rz(-1.4624566) q[3];
sx q[3];
rz(-2.9379454) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48731503) q[0];
sx q[0];
rz(-1.4947816) q[0];
sx q[0];
rz(0.017729433) q[0];
rz(-2.0954633) q[1];
sx q[1];
rz(-1.4132063) q[1];
sx q[1];
rz(-1.0606934) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26955206) q[0];
sx q[0];
rz(-1.6395131) q[0];
sx q[0];
rz(-1.5848573) q[0];
rz(-pi) q[1];
rz(-2.1589375) q[2];
sx q[2];
rz(-1.7870296) q[2];
sx q[2];
rz(-2.4311291) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.880278) q[1];
sx q[1];
rz(-0.78431097) q[1];
sx q[1];
rz(1.8159291) q[1];
rz(3.0396427) q[3];
sx q[3];
rz(-2.6155439) q[3];
sx q[3];
rz(-2.7395484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9866508) q[2];
sx q[2];
rz(-0.81601802) q[2];
sx q[2];
rz(-0.6443392) q[2];
rz(-0.84195697) q[3];
sx q[3];
rz(-1.569845) q[3];
sx q[3];
rz(-0.90698609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7613206) q[0];
sx q[0];
rz(-0.43571061) q[0];
sx q[0];
rz(-0.45183387) q[0];
rz(-2.1322346) q[1];
sx q[1];
rz(-1.6296856) q[1];
sx q[1];
rz(1.570328) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46635755) q[0];
sx q[0];
rz(-1.5769099) q[0];
sx q[0];
rz(2.9714022) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11650916) q[2];
sx q[2];
rz(-1.0986058) q[2];
sx q[2];
rz(-2.0906868) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.81105197) q[1];
sx q[1];
rz(-1.1833541) q[1];
sx q[1];
rz(0.13918332) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5785923) q[3];
sx q[3];
rz(-0.8484133) q[3];
sx q[3];
rz(2.9786199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.62341225) q[2];
sx q[2];
rz(-2.0698915) q[2];
sx q[2];
rz(1.5709467) q[2];
rz(-3.0359641) q[3];
sx q[3];
rz(-0.94625866) q[3];
sx q[3];
rz(-0.51641974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8138206) q[0];
sx q[0];
rz(-1.0186503) q[0];
sx q[0];
rz(-3.0042197) q[0];
rz(2.9227921) q[1];
sx q[1];
rz(-1.7411502) q[1];
sx q[1];
rz(2.2314744) q[1];
rz(1.1280288) q[2];
sx q[2];
rz(-1.3793066) q[2];
sx q[2];
rz(-2.1366091) q[2];
rz(-1.9365334) q[3];
sx q[3];
rz(-2.1954721) q[3];
sx q[3];
rz(1.1385067) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
