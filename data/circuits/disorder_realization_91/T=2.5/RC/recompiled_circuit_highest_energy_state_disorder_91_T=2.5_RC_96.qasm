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
rz(0.36747992) q[0];
sx q[0];
rz(4.9541192) q[0];
sx q[0];
rz(10.911565) q[0];
rz(-1.6410671) q[1];
sx q[1];
rz(2.5379116) q[1];
sx q[1];
rz(11.707468) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14013962) q[0];
sx q[0];
rz(-1.9705904) q[0];
sx q[0];
rz(0.6116914) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98318451) q[2];
sx q[2];
rz(-1.1074083) q[2];
sx q[2];
rz(-2.7962229) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4564701) q[1];
sx q[1];
rz(-2.2500012) q[1];
sx q[1];
rz(0.14107708) q[1];
rz(-2.1726491) q[3];
sx q[3];
rz(-2.7730983) q[3];
sx q[3];
rz(-0.70337765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.19771244) q[2];
sx q[2];
rz(-0.79215017) q[2];
sx q[2];
rz(-2.0749178) q[2];
rz(-3.0471622) q[3];
sx q[3];
rz(-0.61757278) q[3];
sx q[3];
rz(2.7134231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18148947) q[0];
sx q[0];
rz(-2.8747989) q[0];
sx q[0];
rz(1.9285404) q[0];
rz(2.6890697) q[1];
sx q[1];
rz(-0.25233832) q[1];
sx q[1];
rz(2.4620893) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5951341) q[0];
sx q[0];
rz(-1.6041557) q[0];
sx q[0];
rz(-2.8921158) q[0];
x q[1];
rz(1.7361197) q[2];
sx q[2];
rz(-1.2836873) q[2];
sx q[2];
rz(-1.6197255) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7391239) q[1];
sx q[1];
rz(-1.1796412) q[1];
sx q[1];
rz(0.7372589) q[1];
rz(1.3704471) q[3];
sx q[3];
rz(-0.89290038) q[3];
sx q[3];
rz(1.6163449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54000336) q[2];
sx q[2];
rz(-0.82401472) q[2];
sx q[2];
rz(-2.4448815) q[2];
rz(1.2434897) q[3];
sx q[3];
rz(-1.1949298) q[3];
sx q[3];
rz(2.5534326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531279) q[0];
sx q[0];
rz(-1.0364113) q[0];
sx q[0];
rz(-0.82365197) q[0];
rz(2.9846687) q[1];
sx q[1];
rz(-1.4361959) q[1];
sx q[1];
rz(-0.40208152) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0500129) q[0];
sx q[0];
rz(-0.41917831) q[0];
sx q[0];
rz(-1.8989977) q[0];
rz(-pi) q[1];
rz(0.04871647) q[2];
sx q[2];
rz(-1.8557544) q[2];
sx q[2];
rz(1.6866956) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0515251) q[1];
sx q[1];
rz(-2.0323477) q[1];
sx q[1];
rz(-0.36580795) q[1];
x q[2];
rz(-1.0707466) q[3];
sx q[3];
rz(-2.7038136) q[3];
sx q[3];
rz(2.9851346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0154401) q[2];
sx q[2];
rz(-0.93706477) q[2];
sx q[2];
rz(3.0103969) q[2];
rz(-2.0845856) q[3];
sx q[3];
rz(-1.1934692) q[3];
sx q[3];
rz(-1.2098275) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7159202) q[0];
sx q[0];
rz(-1.965006) q[0];
sx q[0];
rz(1.3878938) q[0];
rz(1.7790986) q[1];
sx q[1];
rz(-1.8905996) q[1];
sx q[1];
rz(-1.9349792) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1183035) q[0];
sx q[0];
rz(-2.1739066) q[0];
sx q[0];
rz(1.009737) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64995857) q[2];
sx q[2];
rz(-1.1372677) q[2];
sx q[2];
rz(-0.94160801) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6304508) q[1];
sx q[1];
rz(-1.2690374) q[1];
sx q[1];
rz(-1.6514227) q[1];
rz(-pi) q[2];
rz(-0.73829262) q[3];
sx q[3];
rz(-1.2311934) q[3];
sx q[3];
rz(0.58635724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44986192) q[2];
sx q[2];
rz(-1.8064156) q[2];
sx q[2];
rz(2.195669) q[2];
rz(0.45390421) q[3];
sx q[3];
rz(-0.64954058) q[3];
sx q[3];
rz(0.18979931) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7396963) q[0];
sx q[0];
rz(-1.3119768) q[0];
sx q[0];
rz(0.78347462) q[0];
rz(2.2354194) q[1];
sx q[1];
rz(-2.2515191) q[1];
sx q[1];
rz(2.5772212) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5106192) q[0];
sx q[0];
rz(-1.6897908) q[0];
sx q[0];
rz(-1.4257159) q[0];
rz(-pi) q[1];
rz(0.059877574) q[2];
sx q[2];
rz(-1.1650601) q[2];
sx q[2];
rz(2.6889888) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7184938) q[1];
sx q[1];
rz(-2.1717487) q[1];
sx q[1];
rz(2.3312631) q[1];
rz(-pi) q[2];
rz(-0.18419097) q[3];
sx q[3];
rz(-0.72316611) q[3];
sx q[3];
rz(1.7890872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2169317) q[2];
sx q[2];
rz(-1.4035808) q[2];
sx q[2];
rz(-3.0972287) q[2];
rz(-1.4102604) q[3];
sx q[3];
rz(-0.20874617) q[3];
sx q[3];
rz(-1.9467719) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5656972) q[0];
sx q[0];
rz(-0.79224753) q[0];
sx q[0];
rz(0.76882291) q[0];
rz(-2.5104751) q[1];
sx q[1];
rz(-1.9335856) q[1];
sx q[1];
rz(-0.15359503) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7241868) q[0];
sx q[0];
rz(-2.119078) q[0];
sx q[0];
rz(-2.7566657) q[0];
x q[1];
rz(1.8472438) q[2];
sx q[2];
rz(-0.9381367) q[2];
sx q[2];
rz(-2.0273493) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2163917) q[1];
sx q[1];
rz(-1.4633302) q[1];
sx q[1];
rz(-2.8587384) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77135857) q[3];
sx q[3];
rz(-1.5017209) q[3];
sx q[3];
rz(-1.815064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.067387335) q[2];
sx q[2];
rz(-1.5437443) q[2];
sx q[2];
rz(0.81573168) q[2];
rz(0.054281209) q[3];
sx q[3];
rz(-1.3868325) q[3];
sx q[3];
rz(-1.774971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62189198) q[0];
sx q[0];
rz(-2.0760355) q[0];
sx q[0];
rz(-0.52072293) q[0];
rz(-2.0479274) q[1];
sx q[1];
rz(-0.92898527) q[1];
sx q[1];
rz(1.7589794) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0313697) q[0];
sx q[0];
rz(-2.5532604) q[0];
sx q[0];
rz(-1.8331241) q[0];
rz(2.0215291) q[2];
sx q[2];
rz(-1.0872957) q[2];
sx q[2];
rz(-2.5978501) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3471038) q[1];
sx q[1];
rz(-0.96318775) q[1];
sx q[1];
rz(0.5083063) q[1];
rz(-pi) q[2];
rz(-2.767278) q[3];
sx q[3];
rz(-2.6427794) q[3];
sx q[3];
rz(2.5588644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9426721) q[2];
sx q[2];
rz(-0.3519381) q[2];
sx q[2];
rz(1.675763) q[2];
rz(-0.28907019) q[3];
sx q[3];
rz(-1.3705658) q[3];
sx q[3];
rz(-1.2808965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11956231) q[0];
sx q[0];
rz(-2.1896095) q[0];
sx q[0];
rz(-0.49222487) q[0];
rz(-2.4960663) q[1];
sx q[1];
rz(-1.5461092) q[1];
sx q[1];
rz(-0.92799497) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6027234) q[0];
sx q[0];
rz(-1.4333748) q[0];
sx q[0];
rz(1.2862872) q[0];
x q[1];
rz(-1.1815588) q[2];
sx q[2];
rz(-1.2368918) q[2];
sx q[2];
rz(-2.2452698) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7240844) q[1];
sx q[1];
rz(-0.88005304) q[1];
sx q[1];
rz(-0.35079415) q[1];
x q[2];
rz(0.31854872) q[3];
sx q[3];
rz(-1.9236698) q[3];
sx q[3];
rz(3.1016289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8202028) q[2];
sx q[2];
rz(-1.5265744) q[2];
sx q[2];
rz(0.75902933) q[2];
rz(1.8697033) q[3];
sx q[3];
rz(-1.3262409) q[3];
sx q[3];
rz(2.429764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-1.3439381) q[0];
sx q[0];
rz(-0.54867083) q[0];
sx q[0];
rz(-0.77242533) q[0];
rz(-1.7732357) q[1];
sx q[1];
rz(-1.1027579) q[1];
sx q[1];
rz(-0.30430421) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1436263) q[0];
sx q[0];
rz(-1.7425101) q[0];
sx q[0];
rz(0.15526662) q[0];
rz(-1.5866942) q[2];
sx q[2];
rz(-0.22581183) q[2];
sx q[2];
rz(-0.91561958) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9615677) q[1];
sx q[1];
rz(-1.3912545) q[1];
sx q[1];
rz(0.7431598) q[1];
x q[2];
rz(-1.4449045) q[3];
sx q[3];
rz(-1.9625809) q[3];
sx q[3];
rz(-2.7916933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7336537) q[2];
sx q[2];
rz(-2.759178) q[2];
sx q[2];
rz(-1.1747053) q[2];
rz(2.5041653) q[3];
sx q[3];
rz(-1.8791684) q[3];
sx q[3];
rz(-1.2229961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8680854) q[0];
sx q[0];
rz(-0.33184505) q[0];
sx q[0];
rz(-0.77274957) q[0];
rz(2.6840774) q[1];
sx q[1];
rz(-1.7465778) q[1];
sx q[1];
rz(1.4332019) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49559418) q[0];
sx q[0];
rz(-2.3353205) q[0];
sx q[0];
rz(-0.92732556) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7501306) q[2];
sx q[2];
rz(-2.7731189) q[2];
sx q[2];
rz(0.76900208) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47327504) q[1];
sx q[1];
rz(-0.5328446) q[1];
sx q[1];
rz(-1.024763) q[1];
x q[2];
rz(2.1674311) q[3];
sx q[3];
rz(-1.0785558) q[3];
sx q[3];
rz(-0.31271586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.92051238) q[2];
sx q[2];
rz(-2.1996193) q[2];
sx q[2];
rz(2.4648049) q[2];
rz(1.6109198) q[3];
sx q[3];
rz(-0.30083209) q[3];
sx q[3];
rz(2.0570741) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9250225) q[0];
sx q[0];
rz(-1.8393479) q[0];
sx q[0];
rz(1.7787697) q[0];
rz(-0.88381797) q[1];
sx q[1];
rz(-0.73278905) q[1];
sx q[1];
rz(-0.97272452) q[1];
rz(2.1402609) q[2];
sx q[2];
rz(-1.9591667) q[2];
sx q[2];
rz(-2.9137076) q[2];
rz(-0.45358446) q[3];
sx q[3];
rz(-0.50752487) q[3];
sx q[3];
rz(-2.8625776) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
