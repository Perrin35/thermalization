OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.088012785) q[0];
sx q[0];
rz(3.454257) q[0];
sx q[0];
rz(9.5587048) q[0];
rz(-0.0066095134) q[1];
sx q[1];
rz(-1.5516888) q[1];
sx q[1];
rz(0.42491999) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3021561) q[0];
sx q[0];
rz(-1.470321) q[0];
sx q[0];
rz(-1.6860486) q[0];
x q[1];
rz(2.7578283) q[2];
sx q[2];
rz(-0.94049938) q[2];
sx q[2];
rz(3.0667194) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.16046187) q[1];
sx q[1];
rz(-0.71804249) q[1];
sx q[1];
rz(-2.1001656) q[1];
rz(-pi) q[2];
rz(2.22394) q[3];
sx q[3];
rz(-0.95953959) q[3];
sx q[3];
rz(-0.33456803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46140823) q[2];
sx q[2];
rz(-2.7646061) q[2];
sx q[2];
rz(0.11412966) q[2];
rz(-2.7528609) q[3];
sx q[3];
rz(-0.26624334) q[3];
sx q[3];
rz(-1.3346599) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1330133) q[0];
sx q[0];
rz(-2.7393434) q[0];
sx q[0];
rz(0.026799686) q[0];
rz(-1.1773479) q[1];
sx q[1];
rz(-2.4939996) q[1];
sx q[1];
rz(-3.1039544) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1397495) q[0];
sx q[0];
rz(-1.5573553) q[0];
sx q[0];
rz(2.3670908) q[0];
x q[1];
rz(-2.9178647) q[2];
sx q[2];
rz(-1.3685952) q[2];
sx q[2];
rz(-3.1333609) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7668132) q[1];
sx q[1];
rz(-1.0008924) q[1];
sx q[1];
rz(0.5809231) q[1];
x q[2];
rz(-1.2977799) q[3];
sx q[3];
rz(-2.3721266) q[3];
sx q[3];
rz(-0.96801341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2521952) q[2];
sx q[2];
rz(-0.9844206) q[2];
sx q[2];
rz(-0.21749116) q[2];
rz(2.7018231) q[3];
sx q[3];
rz(-1.35448) q[3];
sx q[3];
rz(-3.0396438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6169823) q[0];
sx q[0];
rz(-3.1341902) q[0];
sx q[0];
rz(2.5819085) q[0];
rz(0.89645487) q[1];
sx q[1];
rz(-0.47210109) q[1];
sx q[1];
rz(-0.40759531) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89002362) q[0];
sx q[0];
rz(-0.59614658) q[0];
sx q[0];
rz(1.7771276) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6569263) q[2];
sx q[2];
rz(-2.4698832) q[2];
sx q[2];
rz(-1.145878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7201336) q[1];
sx q[1];
rz(-1.4602666) q[1];
sx q[1];
rz(-0.46760749) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0225917) q[3];
sx q[3];
rz(-1.8361409) q[3];
sx q[3];
rz(0.010637334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.018230351) q[2];
sx q[2];
rz(-1.9841649) q[2];
sx q[2];
rz(2.1165712) q[2];
rz(-1.7068663) q[3];
sx q[3];
rz(-0.44308174) q[3];
sx q[3];
rz(0.6492492) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0822815) q[0];
sx q[0];
rz(-2.870443) q[0];
sx q[0];
rz(0.49736381) q[0];
rz(1.2936032) q[1];
sx q[1];
rz(-2.135364) q[1];
sx q[1];
rz(-1.9088378) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8393173) q[0];
sx q[0];
rz(-2.1324106) q[0];
sx q[0];
rz(-1.3297434) q[0];
rz(-1.1511939) q[2];
sx q[2];
rz(-2.2299354) q[2];
sx q[2];
rz(1.4105547) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37356578) q[1];
sx q[1];
rz(-1.7737391) q[1];
sx q[1];
rz(0.25217863) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0972872) q[3];
sx q[3];
rz(-2.7409275) q[3];
sx q[3];
rz(1.6687499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5645912) q[2];
sx q[2];
rz(-2.5265054) q[2];
sx q[2];
rz(-0.1423398) q[2];
rz(-1.1764935) q[3];
sx q[3];
rz(-2.1409972) q[3];
sx q[3];
rz(0.4308027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7560526) q[0];
sx q[0];
rz(-0.98243326) q[0];
sx q[0];
rz(-0.2734215) q[0];
rz(0.39971071) q[1];
sx q[1];
rz(-0.81396657) q[1];
sx q[1];
rz(2.8846557) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32460913) q[0];
sx q[0];
rz(-2.5053609) q[0];
sx q[0];
rz(1.0093635) q[0];
rz(2.1295142) q[2];
sx q[2];
rz(-2.693579) q[2];
sx q[2];
rz(-1.5504993) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94019425) q[1];
sx q[1];
rz(-1.640854) q[1];
sx q[1];
rz(0.045658535) q[1];
rz(2.2502348) q[3];
sx q[3];
rz(-1.3472424) q[3];
sx q[3];
rz(2.4628061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5895245) q[2];
sx q[2];
rz(-0.35577154) q[2];
sx q[2];
rz(-2.3815928) q[2];
rz(0.99203569) q[3];
sx q[3];
rz(-1.2090679) q[3];
sx q[3];
rz(1.9973283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8511667) q[0];
sx q[0];
rz(-2.5100584) q[0];
sx q[0];
rz(-2.5354711) q[0];
rz(-2.0947314) q[1];
sx q[1];
rz(-1.4907587) q[1];
sx q[1];
rz(3.0643588) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4299049) q[0];
sx q[0];
rz(-0.042379286) q[0];
sx q[0];
rz(1.872657) q[0];
rz(-pi) q[1];
rz(-1.0995819) q[2];
sx q[2];
rz(-1.2849286) q[2];
sx q[2];
rz(0.8980823) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97925767) q[1];
sx q[1];
rz(-1.1282776) q[1];
sx q[1];
rz(-2.6975432) q[1];
rz(1.3807543) q[3];
sx q[3];
rz(-0.9166719) q[3];
sx q[3];
rz(3.0669756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91654009) q[2];
sx q[2];
rz(-0.70987916) q[2];
sx q[2];
rz(-2.8705961) q[2];
rz(0.58800507) q[3];
sx q[3];
rz(-2.3376412) q[3];
sx q[3];
rz(2.7325381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3410909) q[0];
sx q[0];
rz(-1.611447) q[0];
sx q[0];
rz(-2.8588168) q[0];
rz(-0.68280363) q[1];
sx q[1];
rz(-0.86137259) q[1];
sx q[1];
rz(-2.2883794) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4597854) q[0];
sx q[0];
rz(-2.2553758) q[0];
sx q[0];
rz(0.48779488) q[0];
rz(-pi) q[1];
rz(1.0710083) q[2];
sx q[2];
rz(-0.82497901) q[2];
sx q[2];
rz(1.210809) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2165378) q[1];
sx q[1];
rz(-1.7345456) q[1];
sx q[1];
rz(1.8064966) q[1];
rz(-3.0158132) q[3];
sx q[3];
rz(-1.7656293) q[3];
sx q[3];
rz(-1.8402214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0736488) q[2];
sx q[2];
rz(-2.6239519) q[2];
sx q[2];
rz(-0.29120564) q[2];
rz(-1.2047042) q[3];
sx q[3];
rz(-0.57345814) q[3];
sx q[3];
rz(1.5709491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35029995) q[0];
sx q[0];
rz(-1.5234103) q[0];
sx q[0];
rz(2.5147901) q[0];
rz(2.0794012) q[1];
sx q[1];
rz(-0.79963446) q[1];
sx q[1];
rz(-0.83745426) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27521321) q[0];
sx q[0];
rz(-1.4555706) q[0];
sx q[0];
rz(2.000745) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8026514) q[2];
sx q[2];
rz(-1.3474479) q[2];
sx q[2];
rz(-1.7898343) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4204605) q[1];
sx q[1];
rz(-1.0612121) q[1];
sx q[1];
rz(-2.9886118) q[1];
rz(-pi) q[2];
rz(-1.4222048) q[3];
sx q[3];
rz(-2.0839543) q[3];
sx q[3];
rz(0.10529127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.31967878) q[2];
sx q[2];
rz(-1.2205114) q[2];
sx q[2];
rz(-2.1004045) q[2];
rz(-2.239481) q[3];
sx q[3];
rz(-1.9769042) q[3];
sx q[3];
rz(-0.56727099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5015471) q[0];
sx q[0];
rz(-0.29842672) q[0];
sx q[0];
rz(-0.10064594) q[0];
rz(1.8177265) q[1];
sx q[1];
rz(-0.83815014) q[1];
sx q[1];
rz(-2.5460338) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8190191) q[0];
sx q[0];
rz(-2.3710476) q[0];
sx q[0];
rz(2.3014803) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8331489) q[2];
sx q[2];
rz(-0.67910087) q[2];
sx q[2];
rz(0.51788143) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2145077) q[1];
sx q[1];
rz(-2.2803398) q[1];
sx q[1];
rz(-1.8983215) q[1];
rz(-pi) q[2];
rz(-1.727739) q[3];
sx q[3];
rz(-2.3683386) q[3];
sx q[3];
rz(-2.2048304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9969534) q[2];
sx q[2];
rz(-0.75012952) q[2];
sx q[2];
rz(1.7993125) q[2];
rz(2.4053549) q[3];
sx q[3];
rz(-0.33595556) q[3];
sx q[3];
rz(-0.098522447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017224273) q[0];
sx q[0];
rz(-0.65171826) q[0];
sx q[0];
rz(-0.29618725) q[0];
rz(1.4172957) q[1];
sx q[1];
rz(-2.2139151) q[1];
sx q[1];
rz(-2.9790624) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86315853) q[0];
sx q[0];
rz(-1.5641644) q[0];
sx q[0];
rz(-1.5970105) q[0];
rz(-2.7835566) q[2];
sx q[2];
rz(-2.0581783) q[2];
sx q[2];
rz(-2.0646937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1028984) q[1];
sx q[1];
rz(-1.7084645) q[1];
sx q[1];
rz(-2.9522411) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1381091) q[3];
sx q[3];
rz(-2.0334019) q[3];
sx q[3];
rz(-1.0779276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6235003) q[2];
sx q[2];
rz(-1.2178428) q[2];
sx q[2];
rz(0.45563844) q[2];
rz(1.0738922) q[3];
sx q[3];
rz(-2.5116601) q[3];
sx q[3];
rz(2.5560801) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940654) q[0];
sx q[0];
rz(-1.8210664) q[0];
sx q[0];
rz(-0.6022712) q[0];
rz(-0.31996721) q[1];
sx q[1];
rz(-2.0260369) q[1];
sx q[1];
rz(1.8021348) q[1];
rz(1.6306277) q[2];
sx q[2];
rz(-1.2345805) q[2];
sx q[2];
rz(2.6474093) q[2];
rz(-2.3665041) q[3];
sx q[3];
rz(-1.6476484) q[3];
sx q[3];
rz(-2.021029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
