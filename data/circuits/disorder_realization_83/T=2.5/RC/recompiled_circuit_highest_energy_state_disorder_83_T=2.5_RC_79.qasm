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
rz(2.0496378) q[0];
sx q[0];
rz(-2.0593934) q[0];
sx q[0];
rz(-1.1991294) q[0];
rz(0.24425976) q[1];
sx q[1];
rz(-1.3173988) q[1];
sx q[1];
rz(-0.86950818) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8438709) q[0];
sx q[0];
rz(-1.3605357) q[0];
sx q[0];
rz(-2.3460827) q[0];
rz(0.42596415) q[2];
sx q[2];
rz(-1.1168343) q[2];
sx q[2];
rz(0.23655836) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8972733) q[1];
sx q[1];
rz(-0.88855511) q[1];
sx q[1];
rz(-1.8771434) q[1];
rz(-pi) q[2];
rz(2.364089) q[3];
sx q[3];
rz(-0.99626675) q[3];
sx q[3];
rz(3.0878909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96183744) q[2];
sx q[2];
rz(-1.9840252) q[2];
sx q[2];
rz(-1.301514) q[2];
rz(-2.2573722) q[3];
sx q[3];
rz(-1.0853465) q[3];
sx q[3];
rz(1.4475383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.41534153) q[0];
sx q[0];
rz(-1.911442) q[0];
sx q[0];
rz(-3.0905261) q[0];
rz(0.49634936) q[1];
sx q[1];
rz(-1.8615078) q[1];
sx q[1];
rz(-3.0275717) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3481596) q[0];
sx q[0];
rz(-1.3591477) q[0];
sx q[0];
rz(0.52692174) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5548666) q[2];
sx q[2];
rz(-0.6932879) q[2];
sx q[2];
rz(0.42809286) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.025578827) q[1];
sx q[1];
rz(-0.65372983) q[1];
sx q[1];
rz(2.3313001) q[1];
x q[2];
rz(2.5258371) q[3];
sx q[3];
rz(-2.1020707) q[3];
sx q[3];
rz(2.5221276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.30782917) q[2];
sx q[2];
rz(-1.9139863) q[2];
sx q[2];
rz(0.366079) q[2];
rz(-0.45575538) q[3];
sx q[3];
rz(-2.3594806) q[3];
sx q[3];
rz(-0.2987878) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61641055) q[0];
sx q[0];
rz(-2.1822378) q[0];
sx q[0];
rz(0.33732238) q[0];
rz(-1.6711383) q[1];
sx q[1];
rz(-2.2098139) q[1];
sx q[1];
rz(-0.09387389) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84067837) q[0];
sx q[0];
rz(-2.4196559) q[0];
sx q[0];
rz(-3.0283699) q[0];
x q[1];
rz(-2.6930935) q[2];
sx q[2];
rz(-1.2587768) q[2];
sx q[2];
rz(1.2879368) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15726883) q[1];
sx q[1];
rz(-0.4254057) q[1];
sx q[1];
rz(-1.5509863) q[1];
x q[2];
rz(-1.8444421) q[3];
sx q[3];
rz(-0.85262596) q[3];
sx q[3];
rz(0.20665619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1845392) q[2];
sx q[2];
rz(-2.5856555) q[2];
sx q[2];
rz(2.8957193) q[2];
rz(0.18694123) q[3];
sx q[3];
rz(-1.723946) q[3];
sx q[3];
rz(-2.7174182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4251855) q[0];
sx q[0];
rz(-2.462429) q[0];
sx q[0];
rz(-0.19234046) q[0];
rz(-1.1771857) q[1];
sx q[1];
rz(-2.3277551) q[1];
sx q[1];
rz(2.7331533) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9964938) q[0];
sx q[0];
rz(-1.6276169) q[0];
sx q[0];
rz(-0.15444093) q[0];
x q[1];
rz(2.8772565) q[2];
sx q[2];
rz(-1.2792584) q[2];
sx q[2];
rz(-2.4318926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3689942) q[1];
sx q[1];
rz(-2.1647758) q[1];
sx q[1];
rz(-1.0812378) q[1];
rz(-pi) q[2];
rz(1.3173265) q[3];
sx q[3];
rz(-0.70159021) q[3];
sx q[3];
rz(1.4647558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8583777) q[2];
sx q[2];
rz(-2.242531) q[2];
sx q[2];
rz(-0.176972) q[2];
rz(0.58961287) q[3];
sx q[3];
rz(-1.6556581) q[3];
sx q[3];
rz(-2.1069215) q[3];
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
rz(-1.4347587) q[0];
sx q[0];
rz(-2.28573) q[0];
sx q[0];
rz(0.51704299) q[0];
rz(0.21273461) q[1];
sx q[1];
rz(-2.0566302) q[1];
sx q[1];
rz(-2.3092666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.413899) q[0];
sx q[0];
rz(-3.0988375) q[0];
sx q[0];
rz(0.18563272) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2465785) q[2];
sx q[2];
rz(-0.38215206) q[2];
sx q[2];
rz(-2.7016751) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.019781915) q[1];
sx q[1];
rz(-1.2488718) q[1];
sx q[1];
rz(1.1888629) q[1];
rz(-pi) q[2];
rz(0.30908567) q[3];
sx q[3];
rz(-2.6431572) q[3];
sx q[3];
rz(-0.70685742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0011255) q[2];
sx q[2];
rz(-2.3095864) q[2];
sx q[2];
rz(-2.7531085) q[2];
rz(-1.2515986) q[3];
sx q[3];
rz(-1.3585217) q[3];
sx q[3];
rz(-1.0640915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0851704) q[0];
sx q[0];
rz(-0.52260411) q[0];
sx q[0];
rz(1.1473468) q[0];
rz(1.7583678) q[1];
sx q[1];
rz(-1.5937832) q[1];
sx q[1];
rz(0.17471084) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6679113) q[0];
sx q[0];
rz(-1.5690049) q[0];
sx q[0];
rz(-0.67497336) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60947588) q[2];
sx q[2];
rz(-1.5769488) q[2];
sx q[2];
rz(0.61340082) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6683674) q[1];
sx q[1];
rz(-0.36512926) q[1];
sx q[1];
rz(0.58546807) q[1];
rz(-pi) q[2];
rz(-1.5368323) q[3];
sx q[3];
rz(-1.782825) q[3];
sx q[3];
rz(-2.9370523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2061578) q[2];
sx q[2];
rz(-1.2015227) q[2];
sx q[2];
rz(-2.5344482) q[2];
rz(-0.28572765) q[3];
sx q[3];
rz(-2.5305179) q[3];
sx q[3];
rz(0.40209517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5173335) q[0];
sx q[0];
rz(-1.7835971) q[0];
sx q[0];
rz(-0.72878033) q[0];
rz(1.1392611) q[1];
sx q[1];
rz(-1.1092721) q[1];
sx q[1];
rz(1.8702102) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.073153) q[0];
sx q[0];
rz(-1.5827142) q[0];
sx q[0];
rz(3.1392142) q[0];
x q[1];
rz(2.5554496) q[2];
sx q[2];
rz(-1.8019466) q[2];
sx q[2];
rz(1.4755206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.374124) q[1];
sx q[1];
rz(-2.136655) q[1];
sx q[1];
rz(-1.330862) q[1];
rz(-1.3649105) q[3];
sx q[3];
rz(-0.48501462) q[3];
sx q[3];
rz(2.5554394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.83069688) q[2];
sx q[2];
rz(-1.7252012) q[2];
sx q[2];
rz(0.26965109) q[2];
rz(2.0693178) q[3];
sx q[3];
rz(-2.8282073) q[3];
sx q[3];
rz(2.0872033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4521745) q[0];
sx q[0];
rz(-1.9502689) q[0];
sx q[0];
rz(2.7574975) q[0];
rz(1.1098038) q[1];
sx q[1];
rz(-0.88491076) q[1];
sx q[1];
rz(1.5193411) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1467595) q[0];
sx q[0];
rz(-1.681856) q[0];
sx q[0];
rz(-1.5572118) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7468734) q[2];
sx q[2];
rz(-0.92967722) q[2];
sx q[2];
rz(1.8973998) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.010433179) q[1];
sx q[1];
rz(-1.466806) q[1];
sx q[1];
rz(1.7883744) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7447837) q[3];
sx q[3];
rz(-2.054415) q[3];
sx q[3];
rz(2.015851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9263837) q[2];
sx q[2];
rz(-1.0101725) q[2];
sx q[2];
rz(-2.9877648) q[2];
rz(0.93872968) q[3];
sx q[3];
rz(-1.2490844) q[3];
sx q[3];
rz(1.4332829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4140103) q[0];
sx q[0];
rz(-1.3879956) q[0];
sx q[0];
rz(-0.54186064) q[0];
rz(1.2205203) q[1];
sx q[1];
rz(-0.6747171) q[1];
sx q[1];
rz(-0.064124785) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42196938) q[0];
sx q[0];
rz(-3.0753284) q[0];
sx q[0];
rz(-0.55304773) q[0];
x q[1];
rz(1.4188962) q[2];
sx q[2];
rz(-0.70775251) q[2];
sx q[2];
rz(2.6709556) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.227441) q[1];
sx q[1];
rz(-2.4898448) q[1];
sx q[1];
rz(-0.027795271) q[1];
x q[2];
rz(-1.3702003) q[3];
sx q[3];
rz(-1.163639) q[3];
sx q[3];
rz(-0.98288051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8378143) q[2];
sx q[2];
rz(-2.1393675) q[2];
sx q[2];
rz(0.079806002) q[2];
rz(0.59085733) q[3];
sx q[3];
rz(-2.8320524) q[3];
sx q[3];
rz(-3.0965366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87126842) q[0];
sx q[0];
rz(-2.7239983) q[0];
sx q[0];
rz(0.25750345) q[0];
rz(-2.2389257) q[1];
sx q[1];
rz(-0.84044424) q[1];
sx q[1];
rz(-0.88817516) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0125663) q[0];
sx q[0];
rz(-1.9309224) q[0];
sx q[0];
rz(-8/(11*pi)) q[0];
x q[1];
rz(-0.85827152) q[2];
sx q[2];
rz(-2.2383318) q[2];
sx q[2];
rz(2.7076662) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9714289) q[1];
sx q[1];
rz(-2.2696643) q[1];
sx q[1];
rz(-2.4902043) q[1];
x q[2];
rz(2.2769663) q[3];
sx q[3];
rz(-2.1451604) q[3];
sx q[3];
rz(2.7780617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4070134) q[2];
sx q[2];
rz(-0.96757704) q[2];
sx q[2];
rz(-0.61545294) q[2];
rz(0.19221273) q[3];
sx q[3];
rz(-0.8513611) q[3];
sx q[3];
rz(-1.6063469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94083448) q[0];
sx q[0];
rz(-1.8466908) q[0];
sx q[0];
rz(1.2373663) q[0];
rz(1.0279961) q[1];
sx q[1];
rz(-0.68872394) q[1];
sx q[1];
rz(0.56180305) q[1];
rz(0.061464389) q[2];
sx q[2];
rz(-0.93498421) q[2];
sx q[2];
rz(-1.5056654) q[2];
rz(0.56461188) q[3];
sx q[3];
rz(-2.5228291) q[3];
sx q[3];
rz(2.1822486) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
