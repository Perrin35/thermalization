OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0140822) q[0];
sx q[0];
rz(-0.2427225) q[0];
sx q[0];
rz(2.2665562) q[0];
rz(1.2526644) q[1];
sx q[1];
rz(-1.5046321) q[1];
sx q[1];
rz(0.59706444) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9677154) q[0];
sx q[0];
rz(-1.1861897) q[0];
sx q[0];
rz(3.0498758) q[0];
rz(-pi) q[1];
rz(2.4389285) q[2];
sx q[2];
rz(-1.7345554) q[2];
sx q[2];
rz(2.5644349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.31103926) q[1];
sx q[1];
rz(-0.90021407) q[1];
sx q[1];
rz(1.5064658) q[1];
x q[2];
rz(1.8641169) q[3];
sx q[3];
rz(-2.9676675) q[3];
sx q[3];
rz(-0.21290646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.732932) q[2];
sx q[2];
rz(-3.0568558) q[2];
sx q[2];
rz(1.5366459) q[2];
rz(1.599865) q[3];
sx q[3];
rz(-2.7432975) q[3];
sx q[3];
rz(-2.1275684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23564553) q[0];
sx q[0];
rz(-2.5962317) q[0];
sx q[0];
rz(1.1646618) q[0];
rz(-2.2230478) q[1];
sx q[1];
rz(-0.48857498) q[1];
sx q[1];
rz(-0.69893828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59832788) q[0];
sx q[0];
rz(-0.65552467) q[0];
sx q[0];
rz(-0.40055029) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.318803) q[2];
sx q[2];
rz(-0.70748913) q[2];
sx q[2];
rz(-1.3327662) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0731502) q[1];
sx q[1];
rz(-0.13769503) q[1];
sx q[1];
rz(-0.87611468) q[1];
rz(-pi) q[2];
rz(-2.9869479) q[3];
sx q[3];
rz(-2.9907812) q[3];
sx q[3];
rz(3.1057764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32236448) q[2];
sx q[2];
rz(-1.7027723) q[2];
sx q[2];
rz(1.2054075) q[2];
rz(0.15959218) q[3];
sx q[3];
rz(-1.3887838) q[3];
sx q[3];
rz(3.0157183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50652248) q[0];
sx q[0];
rz(-0.078190088) q[0];
sx q[0];
rz(0.11678326) q[0];
rz(2.0109239) q[1];
sx q[1];
rz(-3.0394381) q[1];
sx q[1];
rz(3.0932313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581586) q[0];
sx q[0];
rz(-1.5003221) q[0];
sx q[0];
rz(1.6638046) q[0];
rz(1.618967) q[2];
sx q[2];
rz(-1.4496231) q[2];
sx q[2];
rz(0.19751422) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0815403) q[1];
sx q[1];
rz(-2.2339666) q[1];
sx q[1];
rz(-0.40388784) q[1];
rz(1.4428407) q[3];
sx q[3];
rz(-1.6098009) q[3];
sx q[3];
rz(1.6774663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9174663) q[2];
sx q[2];
rz(-1.5175061) q[2];
sx q[2];
rz(0.65829128) q[2];
rz(1.0607464) q[3];
sx q[3];
rz(-2.3297533) q[3];
sx q[3];
rz(0.15061024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0105522) q[0];
sx q[0];
rz(-3.0578767) q[0];
sx q[0];
rz(-2.2678099) q[0];
rz(1.0134169) q[1];
sx q[1];
rz(-0.0031009379) q[1];
sx q[1];
rz(-0.82469624) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.254941) q[0];
sx q[0];
rz(-1.5290197) q[0];
sx q[0];
rz(3.1095563) q[0];
x q[1];
rz(1.0295139) q[2];
sx q[2];
rz(-2.3397987) q[2];
sx q[2];
rz(0.0060334671) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4383126) q[1];
sx q[1];
rz(-1.4442128) q[1];
sx q[1];
rz(-1.8476608) q[1];
rz(-2.0788105) q[3];
sx q[3];
rz(-1.0107147) q[3];
sx q[3];
rz(-2.687272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.970845) q[2];
sx q[2];
rz(-2.8837995) q[2];
sx q[2];
rz(-2.8802059) q[2];
rz(1.1970674) q[3];
sx q[3];
rz(-1.3470998) q[3];
sx q[3];
rz(2.4813467) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9561387) q[0];
sx q[0];
rz(-0.11341299) q[0];
sx q[0];
rz(1.0979106) q[0];
rz(2.5838891) q[1];
sx q[1];
rz(-0.028956078) q[1];
sx q[1];
rz(1.8580233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0484491) q[0];
sx q[0];
rz(-1.6225166) q[0];
sx q[0];
rz(0.038695996) q[0];
x q[1];
rz(-1.6465285) q[2];
sx q[2];
rz(-1.8316187) q[2];
sx q[2];
rz(2.7882831) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9030212) q[1];
sx q[1];
rz(-2.0687177) q[1];
sx q[1];
rz(-1.4378888) q[1];
x q[2];
rz(0.8423644) q[3];
sx q[3];
rz(-1.0972365) q[3];
sx q[3];
rz(2.0061988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6534609) q[2];
sx q[2];
rz(-2.2296495) q[2];
sx q[2];
rz(-0.27257356) q[2];
rz(-1.9778947) q[3];
sx q[3];
rz(-1.3397763) q[3];
sx q[3];
rz(2.1986966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1777451) q[0];
sx q[0];
rz(-3.0680532) q[0];
sx q[0];
rz(-0.22200008) q[0];
rz(-1.6674626) q[1];
sx q[1];
rz(-3.1167897) q[1];
sx q[1];
rz(-2.7892392) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9637102) q[0];
sx q[0];
rz(-0.0019638722) q[0];
sx q[0];
rz(-1.4675115) q[0];
x q[1];
rz(-0.9325005) q[2];
sx q[2];
rz(-0.83915448) q[2];
sx q[2];
rz(0.41590235) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7166802) q[1];
sx q[1];
rz(-1.0031478) q[1];
sx q[1];
rz(-0.54089728) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1389987) q[3];
sx q[3];
rz(-0.96240265) q[3];
sx q[3];
rz(-2.268057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12677975) q[2];
sx q[2];
rz(-1.1819906) q[2];
sx q[2];
rz(-0.78753769) q[2];
rz(-0.13134232) q[3];
sx q[3];
rz(-2.1613439) q[3];
sx q[3];
rz(2.378715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.6545559) q[0];
sx q[0];
rz(-3.1058703) q[0];
sx q[0];
rz(-0.63294739) q[0];
rz(2.5542906) q[1];
sx q[1];
rz(-0.018922806) q[1];
sx q[1];
rz(-2.6515554) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5929664) q[0];
sx q[0];
rz(-1.6621309) q[0];
sx q[0];
rz(1.7305602) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3798602) q[2];
sx q[2];
rz(-1.2314738) q[2];
sx q[2];
rz(0.81449984) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.9783352) q[1];
sx q[1];
rz(-0.73860093) q[1];
sx q[1];
rz(0.51335044) q[1];
x q[2];
rz(-2.7825317) q[3];
sx q[3];
rz(-2.0395055) q[3];
sx q[3];
rz(1.7079197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6106762) q[2];
sx q[2];
rz(-2.5848415) q[2];
sx q[2];
rz(-0.69182932) q[2];
rz(2.1107215) q[3];
sx q[3];
rz(-0.94018799) q[3];
sx q[3];
rz(-2.809439) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0523025) q[0];
sx q[0];
rz(-3.0806627) q[0];
sx q[0];
rz(-2.6887509) q[0];
rz(-0.37372681) q[1];
sx q[1];
rz(-3.0174297) q[1];
sx q[1];
rz(2.9492818) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34626353) q[0];
sx q[0];
rz(-1.457231) q[0];
sx q[0];
rz(1.5553484) q[0];
rz(-pi) q[1];
rz(0.65625729) q[2];
sx q[2];
rz(-2.5263737) q[2];
sx q[2];
rz(-0.27042897) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65147644) q[1];
sx q[1];
rz(-1.3785751) q[1];
sx q[1];
rz(1.0438471) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1753331) q[3];
sx q[3];
rz(-0.34869683) q[3];
sx q[3];
rz(-2.6570081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8082491) q[2];
sx q[2];
rz(-0.054017544) q[2];
sx q[2];
rz(0.56600189) q[2];
rz(-2.443215) q[3];
sx q[3];
rz(-1.799182) q[3];
sx q[3];
rz(-3.0799871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9653387) q[0];
sx q[0];
rz(-0.34602556) q[0];
sx q[0];
rz(1.457343) q[0];
rz(2.1894042) q[1];
sx q[1];
rz(-3.1377276) q[1];
sx q[1];
rz(-1.4176518) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7577292) q[0];
sx q[0];
rz(-1.2653906) q[0];
sx q[0];
rz(2.9440434) q[0];
rz(1.046696) q[2];
sx q[2];
rz(-1.5907739) q[2];
sx q[2];
rz(-1.4422356) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9940852) q[1];
sx q[1];
rz(-1.6030585) q[1];
sx q[1];
rz(2.9659531) q[1];
x q[2];
rz(-3.0921008) q[3];
sx q[3];
rz(-1.4386571) q[3];
sx q[3];
rz(1.8232839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5119322) q[2];
sx q[2];
rz(-2.8544482) q[2];
sx q[2];
rz(0.49893898) q[2];
rz(-2.387909) q[3];
sx q[3];
rz(-0.53871173) q[3];
sx q[3];
rz(2.5628569) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27624398) q[0];
sx q[0];
rz(-2.8911599) q[0];
sx q[0];
rz(2.8540386) q[0];
rz(-0.7175256) q[1];
sx q[1];
rz(-3.0375807) q[1];
sx q[1];
rz(1.8431374) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13253015) q[0];
sx q[0];
rz(-0.98041897) q[0];
sx q[0];
rz(2.4590465) q[0];
rz(0.60628225) q[2];
sx q[2];
rz(-2.0673476) q[2];
sx q[2];
rz(0.50722659) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.76261452) q[1];
sx q[1];
rz(-1.839338) q[1];
sx q[1];
rz(1.9767799) q[1];
rz(-1.6631546) q[3];
sx q[3];
rz(-1.8849117) q[3];
sx q[3];
rz(0.11387728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8269877) q[2];
sx q[2];
rz(-1.2890451) q[2];
sx q[2];
rz(-2.2513466) q[2];
rz(-3.1173949) q[3];
sx q[3];
rz(-0.20166339) q[3];
sx q[3];
rz(-2.3307687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496307) q[0];
sx q[0];
rz(-0.55593119) q[0];
sx q[0];
rz(2.2297106) q[0];
rz(2.1610319) q[1];
sx q[1];
rz(-0.33976561) q[1];
sx q[1];
rz(0.38358546) q[1];
rz(3.0684971) q[2];
sx q[2];
rz(-1.2580032) q[2];
sx q[2];
rz(-2.0138665) q[2];
rz(-1.8611363) q[3];
sx q[3];
rz(-1.3383337) q[3];
sx q[3];
rz(0.70575502) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
