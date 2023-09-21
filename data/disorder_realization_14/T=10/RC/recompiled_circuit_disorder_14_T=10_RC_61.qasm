OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4560661) q[0];
sx q[0];
rz(-0.38903061) q[0];
sx q[0];
rz(2.2580137) q[0];
rz(-0.0097302516) q[1];
sx q[1];
rz(-1.4571804) q[1];
sx q[1];
rz(1.943346) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0177901) q[0];
sx q[0];
rz(-1.666981) q[0];
sx q[0];
rz(-2.9642446) q[0];
rz(-0.3281524) q[2];
sx q[2];
rz(-0.80291623) q[2];
sx q[2];
rz(-2.3053055) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40542422) q[1];
sx q[1];
rz(-0.6286469) q[1];
sx q[1];
rz(0.58341649) q[1];
rz(1.7415813) q[3];
sx q[3];
rz(-1.5251953) q[3];
sx q[3];
rz(3.1053229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1951695) q[2];
sx q[2];
rz(-0.98313466) q[2];
sx q[2];
rz(2.9602489) q[2];
rz(-2.8803853) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(-0.75631022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.8392035) q[0];
sx q[0];
rz(-2.8518682) q[0];
sx q[0];
rz(-2.7547577) q[0];
rz(-2.6392) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(1.5418672) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8097336) q[0];
sx q[0];
rz(-0.74685687) q[0];
sx q[0];
rz(-2.5144308) q[0];
rz(-pi) q[1];
rz(0.16427152) q[2];
sx q[2];
rz(-1.1195682) q[2];
sx q[2];
rz(2.3300366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7911885) q[1];
sx q[1];
rz(-1.5807307) q[1];
sx q[1];
rz(-3.105858) q[1];
rz(2.3719671) q[3];
sx q[3];
rz(-2.5105021) q[3];
sx q[3];
rz(-2.1345994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5388422) q[2];
sx q[2];
rz(-2.2637612) q[2];
sx q[2];
rz(1.7269469) q[2];
rz(-0.85033068) q[3];
sx q[3];
rz(-2.7089705) q[3];
sx q[3];
rz(1.6833646) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26281115) q[0];
sx q[0];
rz(-1.5314064) q[0];
sx q[0];
rz(2.5313654) q[0];
rz(1.2894851) q[1];
sx q[1];
rz(-0.97924966) q[1];
sx q[1];
rz(-0.99197018) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6863166) q[0];
sx q[0];
rz(-1.0193829) q[0];
sx q[0];
rz(2.9247012) q[0];
rz(0.66833468) q[2];
sx q[2];
rz(-1.7303581) q[2];
sx q[2];
rz(2.4192686) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.386258) q[1];
sx q[1];
rz(-1.660368) q[1];
sx q[1];
rz(0.44134015) q[1];
rz(0.46996689) q[3];
sx q[3];
rz(-0.50438687) q[3];
sx q[3];
rz(0.72987635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38301864) q[2];
sx q[2];
rz(-2.980361) q[2];
sx q[2];
rz(-2.8733011) q[2];
rz(-0.39408436) q[3];
sx q[3];
rz(-1.9106617) q[3];
sx q[3];
rz(2.9530853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2878993) q[0];
sx q[0];
rz(-0.51369602) q[0];
sx q[0];
rz(-0.40507856) q[0];
rz(2.4515117) q[1];
sx q[1];
rz(-1.9837374) q[1];
sx q[1];
rz(2.4437723) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989379) q[0];
sx q[0];
rz(-1.5942759) q[0];
sx q[0];
rz(1.6606746) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54450808) q[2];
sx q[2];
rz(-1.9015357) q[2];
sx q[2];
rz(2.7059908) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7792369) q[1];
sx q[1];
rz(-0.52473611) q[1];
sx q[1];
rz(1.1974105) q[1];
rz(-2.1711369) q[3];
sx q[3];
rz(-0.49649039) q[3];
sx q[3];
rz(1.1288201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4219389) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(-1.6323803) q[2];
rz(-2.737282) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(1.4782762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25800911) q[0];
sx q[0];
rz(-1.785935) q[0];
sx q[0];
rz(1.0255381) q[0];
rz(0.57199663) q[1];
sx q[1];
rz(-1.0943741) q[1];
sx q[1];
rz(0.62932032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.660491) q[0];
sx q[0];
rz(-1.8697303) q[0];
sx q[0];
rz(-2.7116508) q[0];
rz(2.8867678) q[2];
sx q[2];
rz(-2.2307768) q[2];
sx q[2];
rz(-3.0072336) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.021796062) q[1];
sx q[1];
rz(-0.52638678) q[1];
sx q[1];
rz(0.22967931) q[1];
rz(-pi) q[2];
x q[2];
rz(2.878177) q[3];
sx q[3];
rz(-0.98589555) q[3];
sx q[3];
rz(-1.5358621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5489674) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(-2.7992451) q[2];
rz(1.3458378) q[3];
sx q[3];
rz(-1.6941518) q[3];
sx q[3];
rz(2.9455744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.489007) q[0];
sx q[0];
rz(-1.9056029) q[0];
sx q[0];
rz(-0.18519369) q[0];
rz(-1.7350896) q[1];
sx q[1];
rz(-2.0506737) q[1];
sx q[1];
rz(-1.7746183) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.470984) q[0];
sx q[0];
rz(-2.1489722) q[0];
sx q[0];
rz(1.3060119) q[0];
x q[1];
rz(-2.30079) q[2];
sx q[2];
rz(-1.7992939) q[2];
sx q[2];
rz(-0.92323869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0095014) q[1];
sx q[1];
rz(-0.17772929) q[1];
sx q[1];
rz(-1.0264261) q[1];
rz(-2.5227138) q[3];
sx q[3];
rz(-0.66580171) q[3];
sx q[3];
rz(0.57951365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8967445) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(2.2582167) q[2];
rz(-1.7287438) q[3];
sx q[3];
rz(-0.69245517) q[3];
sx q[3];
rz(-2.3278055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25154034) q[0];
sx q[0];
rz(-0.10245704) q[0];
sx q[0];
rz(-1.2782156) q[0];
rz(3.1037519) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(-1.7657123) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8517075) q[0];
sx q[0];
rz(-2.2635248) q[0];
sx q[0];
rz(-0.38498199) q[0];
x q[1];
rz(3.0817229) q[2];
sx q[2];
rz(-1.0979568) q[2];
sx q[2];
rz(-1.5945895) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.059767698) q[1];
sx q[1];
rz(-1.6762814) q[1];
sx q[1];
rz(-0.11591537) q[1];
rz(-pi) q[2];
rz(-1.901058) q[3];
sx q[3];
rz(-1.6263279) q[3];
sx q[3];
rz(-0.33952573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6841131) q[2];
sx q[2];
rz(-2.423968) q[2];
sx q[2];
rz(-2.4105371) q[2];
rz(3.030792) q[3];
sx q[3];
rz(-1.585107) q[3];
sx q[3];
rz(2.4462162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5450127) q[0];
sx q[0];
rz(-2.2667363) q[0];
sx q[0];
rz(0.73356432) q[0];
rz(0.60797524) q[1];
sx q[1];
rz(-1.9476451) q[1];
sx q[1];
rz(0.2342934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83345862) q[0];
sx q[0];
rz(-1.8292384) q[0];
sx q[0];
rz(2.0636369) q[0];
x q[1];
rz(0.90227622) q[2];
sx q[2];
rz(-0.96792816) q[2];
sx q[2];
rz(0.55842802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81208166) q[1];
sx q[1];
rz(-2.3595464) q[1];
sx q[1];
rz(1.300632) q[1];
rz(1.0844564) q[3];
sx q[3];
rz(-0.64389766) q[3];
sx q[3];
rz(2.7542496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12604788) q[2];
sx q[2];
rz(-1.8323106) q[2];
sx q[2];
rz(-1.7162494) q[2];
rz(1.4632633) q[3];
sx q[3];
rz(-2.3571456) q[3];
sx q[3];
rz(1.6459758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54206806) q[0];
sx q[0];
rz(-2.8064089) q[0];
sx q[0];
rz(-1.9482127) q[0];
rz(-1.2606196) q[1];
sx q[1];
rz(-1.7766989) q[1];
sx q[1];
rz(2.1967922) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.955907) q[0];
sx q[0];
rz(-2.3392896) q[0];
sx q[0];
rz(-2.5923652) q[0];
x q[1];
rz(-2.8464727) q[2];
sx q[2];
rz(-1.9502389) q[2];
sx q[2];
rz(1.5038568) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.33465696) q[1];
sx q[1];
rz(-1.2922704) q[1];
sx q[1];
rz(1.3929699) q[1];
rz(0.19712574) q[3];
sx q[3];
rz(-1.9204431) q[3];
sx q[3];
rz(0.2998578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9251359) q[2];
sx q[2];
rz(-1.0711121) q[2];
sx q[2];
rz(2.8533868) q[2];
rz(-0.47973412) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(1.5079927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0297246) q[0];
sx q[0];
rz(-0.27619633) q[0];
sx q[0];
rz(0.9129886) q[0];
rz(-2.7669725) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(0.8909117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8091781) q[0];
sx q[0];
rz(-0.42254585) q[0];
sx q[0];
rz(-0.64230625) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.183379) q[2];
sx q[2];
rz(-1.3365067) q[2];
sx q[2];
rz(2.2189552) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9526457) q[1];
sx q[1];
rz(-1.2332321) q[1];
sx q[1];
rz(-1.7811716) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18817801) q[3];
sx q[3];
rz(-1.6450226) q[3];
sx q[3];
rz(-2.1256413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23641071) q[2];
sx q[2];
rz(-2.2832401) q[2];
sx q[2];
rz(-2.6399844) q[2];
rz(-1.2891399) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(-0.44617173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5065153) q[0];
sx q[0];
rz(-1.7000533) q[0];
sx q[0];
rz(0.62361367) q[0];
rz(-1.1322017) q[1];
sx q[1];
rz(-2.3846346) q[1];
sx q[1];
rz(0.089288575) q[1];
rz(1.4176462) q[2];
sx q[2];
rz(-1.8344804) q[2];
sx q[2];
rz(-1.9305965) q[2];
rz(-1.489747) q[3];
sx q[3];
rz(-2.5935018) q[3];
sx q[3];
rz(0.58542585) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
