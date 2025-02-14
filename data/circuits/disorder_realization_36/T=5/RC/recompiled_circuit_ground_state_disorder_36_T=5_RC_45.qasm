OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1898243) q[0];
sx q[0];
rz(-0.70280743) q[0];
sx q[0];
rz(1.2220609) q[0];
rz(3.089978) q[1];
sx q[1];
rz(-1.6366704) q[1];
sx q[1];
rz(1.6610891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9771271) q[0];
sx q[0];
rz(-1.7344622) q[0];
sx q[0];
rz(0.90552355) q[0];
rz(-pi) q[1];
rz(0.10926508) q[2];
sx q[2];
rz(-1.7393629) q[2];
sx q[2];
rz(1.8336465) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6340629) q[1];
sx q[1];
rz(-2.6802217) q[1];
sx q[1];
rz(2.7410236) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8454523) q[3];
sx q[3];
rz(-2.5293455) q[3];
sx q[3];
rz(2.1795764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1808971) q[2];
sx q[2];
rz(-2.3263859) q[2];
sx q[2];
rz(2.1107471) q[2];
rz(-2.893462) q[3];
sx q[3];
rz(-1.7957567) q[3];
sx q[3];
rz(-0.26721755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2613075) q[0];
sx q[0];
rz(-1.945865) q[0];
sx q[0];
rz(-1.1241359) q[0];
rz(1.0719489) q[1];
sx q[1];
rz(-1.5644466) q[1];
sx q[1];
rz(-0.42676485) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2248707) q[0];
sx q[0];
rz(-1.5735177) q[0];
sx q[0];
rz(2.3140074) q[0];
x q[1];
rz(-0.63964455) q[2];
sx q[2];
rz(-0.42176127) q[2];
sx q[2];
rz(-2.1157921) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.427095) q[1];
sx q[1];
rz(-2.2648304) q[1];
sx q[1];
rz(1.4301894) q[1];
x q[2];
rz(0.20189607) q[3];
sx q[3];
rz(-1.1888973) q[3];
sx q[3];
rz(-2.524162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4830604) q[2];
sx q[2];
rz(-0.51653647) q[2];
sx q[2];
rz(0.10022441) q[2];
rz(-1.077486) q[3];
sx q[3];
rz(-1.5651549) q[3];
sx q[3];
rz(1.9061576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60233068) q[0];
sx q[0];
rz(-0.97719181) q[0];
sx q[0];
rz(1.4343028) q[0];
rz(-0.82646787) q[1];
sx q[1];
rz(-1.4974599) q[1];
sx q[1];
rz(-0.43063393) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42303681) q[0];
sx q[0];
rz(-1.406975) q[0];
sx q[0];
rz(2.8680766) q[0];
rz(1.1330092) q[2];
sx q[2];
rz(-2.8418645) q[2];
sx q[2];
rz(2.9309811) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.90884028) q[1];
sx q[1];
rz(-1.1758057) q[1];
sx q[1];
rz(-0.86887143) q[1];
rz(-pi) q[2];
rz(-2.1149733) q[3];
sx q[3];
rz(-1.2261571) q[3];
sx q[3];
rz(-1.0349864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87381252) q[2];
sx q[2];
rz(-1.9766821) q[2];
sx q[2];
rz(-2.5679892) q[2];
rz(0.38392797) q[3];
sx q[3];
rz(-1.826518) q[3];
sx q[3];
rz(0.65281502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7225994) q[0];
sx q[0];
rz(-2.3498131) q[0];
sx q[0];
rz(-1.0571887) q[0];
rz(-0.81622299) q[1];
sx q[1];
rz(-2.1482601) q[1];
sx q[1];
rz(-2.9073471) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5200096) q[0];
sx q[0];
rz(-2.710973) q[0];
sx q[0];
rz(-1.3060116) q[0];
rz(-pi) q[1];
rz(-1.044687) q[2];
sx q[2];
rz(-1.1981694) q[2];
sx q[2];
rz(1.5043196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1409451) q[1];
sx q[1];
rz(-1.06303) q[1];
sx q[1];
rz(2.341048) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0660421) q[3];
sx q[3];
rz(-1.9947933) q[3];
sx q[3];
rz(0.056886176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8546042) q[2];
sx q[2];
rz(-1.333326) q[2];
sx q[2];
rz(0.90281478) q[2];
rz(1.3826987) q[3];
sx q[3];
rz(-2.8185676) q[3];
sx q[3];
rz(2.2949016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46120241) q[0];
sx q[0];
rz(-1.8620551) q[0];
sx q[0];
rz(-2.4317106) q[0];
rz(-0.4711802) q[1];
sx q[1];
rz(-1.9147562) q[1];
sx q[1];
rz(0.064362854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.956139) q[0];
sx q[0];
rz(-2.6366933) q[0];
sx q[0];
rz(1.8527777) q[0];
rz(2.9079403) q[2];
sx q[2];
rz(-1.086768) q[2];
sx q[2];
rz(0.40586995) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.90366615) q[1];
sx q[1];
rz(-1.1519379) q[1];
sx q[1];
rz(-0.55540696) q[1];
rz(-pi) q[2];
rz(0.10897763) q[3];
sx q[3];
rz(-0.087317467) q[3];
sx q[3];
rz(-0.84978896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71089661) q[2];
sx q[2];
rz(-0.5527834) q[2];
sx q[2];
rz(-2.1882679) q[2];
rz(0.55771762) q[3];
sx q[3];
rz(-1.6744813) q[3];
sx q[3];
rz(-0.52887708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.8271269) q[0];
sx q[0];
rz(-2.0966661) q[0];
sx q[0];
rz(-2.1449828) q[0];
rz(-0.059545513) q[1];
sx q[1];
rz(-2.153502) q[1];
sx q[1];
rz(-1.7893808) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9080862) q[0];
sx q[0];
rz(-0.75859374) q[0];
sx q[0];
rz(-0.80484606) q[0];
x q[1];
rz(-1.900104) q[2];
sx q[2];
rz(-0.8267459) q[2];
sx q[2];
rz(-2.0659826) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3660197) q[1];
sx q[1];
rz(-1.8692353) q[1];
sx q[1];
rz(-1.8587023) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7963403) q[3];
sx q[3];
rz(-0.94178978) q[3];
sx q[3];
rz(-2.0382413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.38702866) q[2];
sx q[2];
rz(-1.2754385) q[2];
sx q[2];
rz(-0.55956364) q[2];
rz(2.3061421) q[3];
sx q[3];
rz(-2.7159034) q[3];
sx q[3];
rz(1.3028418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68592042) q[0];
sx q[0];
rz(-2.4113825) q[0];
sx q[0];
rz(0.36668229) q[0];
rz(-2.2600251) q[1];
sx q[1];
rz(-2.5036948) q[1];
sx q[1];
rz(-1.7104023) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0019497) q[0];
sx q[0];
rz(-2.3649923) q[0];
sx q[0];
rz(-0.5712255) q[0];
rz(-pi) q[1];
rz(1.9164247) q[2];
sx q[2];
rz(-0.91424886) q[2];
sx q[2];
rz(0.57166568) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6122653) q[1];
sx q[1];
rz(-3.021486) q[1];
sx q[1];
rz(-0.24397759) q[1];
x q[2];
rz(1.335698) q[3];
sx q[3];
rz(-0.92609105) q[3];
sx q[3];
rz(2.7489782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7876579) q[2];
sx q[2];
rz(-2.6255609) q[2];
sx q[2];
rz(-2.3616974) q[2];
rz(2.6025313) q[3];
sx q[3];
rz(-1.6293679) q[3];
sx q[3];
rz(-0.61162925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014538177) q[0];
sx q[0];
rz(-0.67033613) q[0];
sx q[0];
rz(0.77914733) q[0];
rz(0.51891333) q[1];
sx q[1];
rz(-1.2823558) q[1];
sx q[1];
rz(-0.40294495) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76501827) q[0];
sx q[0];
rz(-0.85581644) q[0];
sx q[0];
rz(0.9299703) q[0];
rz(-pi) q[1];
rz(2.1261931) q[2];
sx q[2];
rz(-2.5882373) q[2];
sx q[2];
rz(2.8560658) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.25499763) q[1];
sx q[1];
rz(-0.95980058) q[1];
sx q[1];
rz(-2.0605036) q[1];
rz(1.1717848) q[3];
sx q[3];
rz(-2.5305989) q[3];
sx q[3];
rz(2.4965167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11401033) q[2];
sx q[2];
rz(-1.9738013) q[2];
sx q[2];
rz(-2.0243417) q[2];
rz(2.4452325) q[3];
sx q[3];
rz(-1.1098692) q[3];
sx q[3];
rz(1.2874359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6663412) q[0];
sx q[0];
rz(-2.9069558) q[0];
sx q[0];
rz(0.41467211) q[0];
rz(-2.0798202) q[1];
sx q[1];
rz(-0.88645187) q[1];
sx q[1];
rz(-2.3019703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3436326) q[0];
sx q[0];
rz(-2.0164444) q[0];
sx q[0];
rz(-0.4526505) q[0];
x q[1];
rz(2.4213441) q[2];
sx q[2];
rz(-2.0145085) q[2];
sx q[2];
rz(0.72138471) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4794107) q[1];
sx q[1];
rz(-2.4644797) q[1];
sx q[1];
rz(0.22775316) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3247588) q[3];
sx q[3];
rz(-1.2241505) q[3];
sx q[3];
rz(-2.3510166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.35905579) q[2];
sx q[2];
rz(-1.7724719) q[2];
sx q[2];
rz(-3.0699406) q[2];
rz(-3.0960633) q[3];
sx q[3];
rz(-1.9436911) q[3];
sx q[3];
rz(0.45904747) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5422106) q[0];
sx q[0];
rz(-2.5686503) q[0];
sx q[0];
rz(-1.0546767) q[0];
rz(3.1090464) q[1];
sx q[1];
rz(-0.97144214) q[1];
sx q[1];
rz(-0.5221101) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5112177) q[0];
sx q[0];
rz(-1.2652093) q[0];
sx q[0];
rz(0.60141984) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1565151) q[2];
sx q[2];
rz(-2.217948) q[2];
sx q[2];
rz(3.0708608) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.594211) q[1];
sx q[1];
rz(-1.1923596) q[1];
sx q[1];
rz(-0.40706472) q[1];
x q[2];
rz(-0.22618146) q[3];
sx q[3];
rz(-2.1070679) q[3];
sx q[3];
rz(-2.8482343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9944428) q[2];
sx q[2];
rz(-2.9394737) q[2];
sx q[2];
rz(-1.408255) q[2];
rz(-1.555892) q[3];
sx q[3];
rz(-0.2318016) q[3];
sx q[3];
rz(2.2304992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1237352) q[0];
sx q[0];
rz(-1.8386848) q[0];
sx q[0];
rz(0.89717502) q[0];
rz(-0.97902117) q[1];
sx q[1];
rz(-1.9985825) q[1];
sx q[1];
rz(-0.40252007) q[1];
rz(-2.945902) q[2];
sx q[2];
rz(-1.2441845) q[2];
sx q[2];
rz(1.2938538) q[2];
rz(-1.3835945) q[3];
sx q[3];
rz(-0.91619195) q[3];
sx q[3];
rz(-0.34294101) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
