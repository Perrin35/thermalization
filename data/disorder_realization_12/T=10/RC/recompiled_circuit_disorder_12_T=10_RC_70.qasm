OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2965887) q[0];
sx q[0];
rz(3.8656524) q[0];
sx q[0];
rz(11.081628) q[0];
rz(1.2031263) q[1];
sx q[1];
rz(-0.523518) q[1];
sx q[1];
rz(2.2533921) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22904299) q[0];
sx q[0];
rz(-0.17670822) q[0];
sx q[0];
rz(0.25632174) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16562478) q[2];
sx q[2];
rz(-1.0399482) q[2];
sx q[2];
rz(-2.3373375) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2792714) q[1];
sx q[1];
rz(-0.16695484) q[1];
sx q[1];
rz(3.1382986) q[1];
rz(-pi) q[2];
rz(-1.6039443) q[3];
sx q[3];
rz(-0.87246934) q[3];
sx q[3];
rz(-0.87227277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.858294) q[2];
sx q[2];
rz(-2.7259939) q[2];
sx q[2];
rz(2.0236012) q[2];
rz(-0.14532146) q[3];
sx q[3];
rz(-1.5829007) q[3];
sx q[3];
rz(-3.0956691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685101) q[0];
sx q[0];
rz(-1.2643603) q[0];
sx q[0];
rz(0.65482393) q[0];
rz(-1.9251992) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(-3.0156946) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1830131) q[0];
sx q[0];
rz(-0.61972451) q[0];
sx q[0];
rz(-0.89967863) q[0];
rz(-0.84463859) q[2];
sx q[2];
rz(-0.8234878) q[2];
sx q[2];
rz(0.57074947) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.8996846) q[1];
sx q[1];
rz(-0.13026127) q[1];
sx q[1];
rz(1.3379407) q[1];
x q[2];
rz(-1.5829854) q[3];
sx q[3];
rz(-2.4180531) q[3];
sx q[3];
rz(0.42750588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.048432365) q[2];
sx q[2];
rz(-1.9530714) q[2];
sx q[2];
rz(0.38802567) q[2];
rz(-1.7175425) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(-0.58732906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5557142) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(-3.0622603) q[0];
rz(0.084005984) q[1];
sx q[1];
rz(-0.80291286) q[1];
sx q[1];
rz(1.1598587) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1201598) q[0];
sx q[0];
rz(-2.624357) q[0];
sx q[0];
rz(1.73818) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94295393) q[2];
sx q[2];
rz(-1.9938333) q[2];
sx q[2];
rz(-1.1689651) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68418903) q[1];
sx q[1];
rz(-2.0261129) q[1];
sx q[1];
rz(1.9940358) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0252762) q[3];
sx q[3];
rz(-1.5226411) q[3];
sx q[3];
rz(0.84850509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40538654) q[2];
sx q[2];
rz(-1.8182886) q[2];
sx q[2];
rz(0.1082871) q[2];
rz(2.4978499) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(0.61557499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.9765587) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(-1.4820341) q[0];
rz(-2.4404793) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(0.28465095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2157001) q[0];
sx q[0];
rz(-1.6452351) q[0];
sx q[0];
rz(2.0490993) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95732032) q[2];
sx q[2];
rz(-2.7942604) q[2];
sx q[2];
rz(-2.6383102) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33672562) q[1];
sx q[1];
rz(-1.6972099) q[1];
sx q[1];
rz(0.079041914) q[1];
rz(-0.90162006) q[3];
sx q[3];
rz(-1.9603143) q[3];
sx q[3];
rz(-0.31182409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.231679) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(-0.56751928) q[2];
rz(0.41401687) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(1.1119941) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48859566) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(-1.8150785) q[0];
rz(1.9891706) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(-0.93793905) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91468231) q[0];
sx q[0];
rz(-1.4786647) q[0];
sx q[0];
rz(3.1196306) q[0];
rz(2.6229834) q[2];
sx q[2];
rz(-2.0687639) q[2];
sx q[2];
rz(-1.2155611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.93463072) q[1];
sx q[1];
rz(-0.043255581) q[1];
sx q[1];
rz(0.91911493) q[1];
x q[2];
rz(1.4411079) q[3];
sx q[3];
rz(-1.218759) q[3];
sx q[3];
rz(0.45613939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1753297) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(-2.1002634) q[2];
rz(1.8390309) q[3];
sx q[3];
rz(-1.1338736) q[3];
sx q[3];
rz(-1.8235824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43907169) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(-0.55737108) q[0];
rz(2.5769261) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(0.55647892) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59069955) q[0];
sx q[0];
rz(-1.7771582) q[0];
sx q[0];
rz(-0.38462374) q[0];
rz(2.9678992) q[2];
sx q[2];
rz(-0.86280338) q[2];
sx q[2];
rz(-0.049335418) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4887052) q[1];
sx q[1];
rz(-2.8294551) q[1];
sx q[1];
rz(2.1991792) q[1];
rz(0.77575404) q[3];
sx q[3];
rz(-1.6072825) q[3];
sx q[3];
rz(-0.072582399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6094728) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(1.9167985) q[2];
rz(-1.6312284) q[3];
sx q[3];
rz(-1.3204201) q[3];
sx q[3];
rz(-1.5009376) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0790134) q[0];
sx q[0];
rz(-1.2591079) q[0];
sx q[0];
rz(0.042908948) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-0.62364548) q[1];
sx q[1];
rz(-2.6409805) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6635054) q[0];
sx q[0];
rz(-0.18112077) q[0];
sx q[0];
rz(-1.3785133) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.074965) q[2];
sx q[2];
rz(-2.5410286) q[2];
sx q[2];
rz(1.9753089) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73270479) q[1];
sx q[1];
rz(-1.171247) q[1];
sx q[1];
rz(-1.7016181) q[1];
rz(1.9686437) q[3];
sx q[3];
rz(-0.85507353) q[3];
sx q[3];
rz(1.0637103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6033972) q[2];
sx q[2];
rz(-2.0624702) q[2];
sx q[2];
rz(-2.4659757) q[2];
rz(2.7006941) q[3];
sx q[3];
rz(-1.7023804) q[3];
sx q[3];
rz(0.52136695) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0764517) q[0];
sx q[0];
rz(-0.32442176) q[0];
sx q[0];
rz(2.0741529) q[0];
rz(-0.43287977) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(1.6759466) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6976801) q[0];
sx q[0];
rz(-1.4364103) q[0];
sx q[0];
rz(1.105955) q[0];
rz(1.6797811) q[2];
sx q[2];
rz(-1.3853067) q[2];
sx q[2];
rz(0.89599228) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.17897478) q[1];
sx q[1];
rz(-1.0777506) q[1];
sx q[1];
rz(-0.38260539) q[1];
x q[2];
rz(-0.64126863) q[3];
sx q[3];
rz(-1.0675758) q[3];
sx q[3];
rz(-0.36676952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.33891588) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(-1.6652997) q[2];
rz(0.87351292) q[3];
sx q[3];
rz(-2.2647808) q[3];
sx q[3];
rz(0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8840238) q[0];
sx q[0];
rz(-0.57618657) q[0];
sx q[0];
rz(2.4043758) q[0];
rz(-0.018741477) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(0.92528701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7172456) q[0];
sx q[0];
rz(-0.92110094) q[0];
sx q[0];
rz(2.2109277) q[0];
x q[1];
rz(0.66321744) q[2];
sx q[2];
rz(-1.6753917) q[2];
sx q[2];
rz(-1.6246206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9463897) q[1];
sx q[1];
rz(-2.9159947) q[1];
sx q[1];
rz(-0.619508) q[1];
x q[2];
rz(0.85261811) q[3];
sx q[3];
rz(-1.8947621) q[3];
sx q[3];
rz(-2.1924803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76901889) q[2];
sx q[2];
rz(-1.5344658) q[2];
sx q[2];
rz(-2.1378689) q[2];
rz(3.051493) q[3];
sx q[3];
rz(-0.026244791) q[3];
sx q[3];
rz(-1.1130921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039624778) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(1.2596624) q[0];
rz(-2.8758077) q[1];
sx q[1];
rz(-2.5140285) q[1];
sx q[1];
rz(2.3840747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3525317) q[0];
sx q[0];
rz(-1.4333945) q[0];
sx q[0];
rz(-1.3760516) q[0];
x q[1];
rz(0.9805571) q[2];
sx q[2];
rz(-1.7014116) q[2];
sx q[2];
rz(-0.57934258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47093686) q[1];
sx q[1];
rz(-0.48479143) q[1];
sx q[1];
rz(-2.0623341) q[1];
x q[2];
rz(2.431589) q[3];
sx q[3];
rz(-1.5922976) q[3];
sx q[3];
rz(-1.8332675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1511128) q[2];
sx q[2];
rz(-1.8258784) q[2];
sx q[2];
rz(2.347836) q[2];
rz(0.63888597) q[3];
sx q[3];
rz(-0.34404889) q[3];
sx q[3];
rz(-1.0206153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6486075) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(-1.5007301) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(-2.5095148) q[2];
sx q[2];
rz(-0.52543228) q[2];
sx q[2];
rz(0.57731522) q[2];
rz(1.100148) q[3];
sx q[3];
rz(-0.63334076) q[3];
sx q[3];
rz(0.28950194) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
