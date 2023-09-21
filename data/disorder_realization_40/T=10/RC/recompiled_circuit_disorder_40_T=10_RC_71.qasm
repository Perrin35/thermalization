OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5326795) q[0];
sx q[0];
rz(-2.764954) q[0];
sx q[0];
rz(-0.11178804) q[0];
rz(1.6821661) q[1];
sx q[1];
rz(-1.4844126) q[1];
sx q[1];
rz(2.9878374) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.114404) q[0];
sx q[0];
rz(-1.9986885) q[0];
sx q[0];
rz(-1.2501636) q[0];
x q[1];
rz(-0.076924952) q[2];
sx q[2];
rz(-1.3300606) q[2];
sx q[2];
rz(-1.7488637) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0220713) q[1];
sx q[1];
rz(-1.540456) q[1];
sx q[1];
rz(-1.8896709) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8044756) q[3];
sx q[3];
rz(-1.9883336) q[3];
sx q[3];
rz(1.9330213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.443632) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(1.4367746) q[2];
rz(0.73389655) q[3];
sx q[3];
rz(-1.5926444) q[3];
sx q[3];
rz(-0.51600391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2519418) q[0];
sx q[0];
rz(-1.9152859) q[0];
sx q[0];
rz(-0.92457986) q[0];
rz(-2.1444767) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(1.8181713) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2307777) q[0];
sx q[0];
rz(-1.1417023) q[0];
sx q[0];
rz(0.50165117) q[0];
rz(-0.1214059) q[2];
sx q[2];
rz(-2.7896023) q[2];
sx q[2];
rz(-2.6999203) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70222774) q[1];
sx q[1];
rz(-0.80184466) q[1];
sx q[1];
rz(0.37978362) q[1];
rz(-pi) q[2];
rz(-0.99625846) q[3];
sx q[3];
rz(-1.2401476) q[3];
sx q[3];
rz(1.4853256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9374342) q[2];
sx q[2];
rz(-1.5519451) q[2];
sx q[2];
rz(-2.3834174) q[2];
rz(0.6289064) q[3];
sx q[3];
rz(-0.40142504) q[3];
sx q[3];
rz(-1.1531856) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4085098) q[0];
sx q[0];
rz(-1.2671616) q[0];
sx q[0];
rz(0.68840233) q[0];
rz(-3.0738661) q[1];
sx q[1];
rz(-1.7522782) q[1];
sx q[1];
rz(-2.6115131) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2358658) q[0];
sx q[0];
rz(-1.6775963) q[0];
sx q[0];
rz(-1.8624767) q[0];
x q[1];
rz(-2.2377551) q[2];
sx q[2];
rz(-2.505216) q[2];
sx q[2];
rz(-2.9449376) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.053777545) q[1];
sx q[1];
rz(-2.3936845) q[1];
sx q[1];
rz(2.0762216) q[1];
rz(-pi) q[2];
rz(-1.3313053) q[3];
sx q[3];
rz(-1.4756087) q[3];
sx q[3];
rz(0.59323192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3337341) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(-2.2401436) q[2];
rz(0.83550134) q[3];
sx q[3];
rz(-1.6136026) q[3];
sx q[3];
rz(-1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9505342) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(-2.3572671) q[0];
rz(3.0803608) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(-0.13664666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.559583) q[0];
sx q[0];
rz(-2.0499381) q[0];
sx q[0];
rz(2.9942306) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6776671) q[2];
sx q[2];
rz(-1.8564965) q[2];
sx q[2];
rz(-2.6505016) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0285443) q[1];
sx q[1];
rz(-1.5946769) q[1];
sx q[1];
rz(-1.5024115) q[1];
rz(-pi) q[2];
rz(3.0487719) q[3];
sx q[3];
rz(-2.1448359) q[3];
sx q[3];
rz(1.012158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0121997) q[2];
sx q[2];
rz(-2.2012074) q[2];
sx q[2];
rz(-2.5811035) q[2];
rz(-3.1292606) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(-1.0906609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.43301582) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(2.3204455) q[0];
rz(0.87617809) q[1];
sx q[1];
rz(-2.2416302) q[1];
sx q[1];
rz(-1.7339773) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27027425) q[0];
sx q[0];
rz(-0.5373913) q[0];
sx q[0];
rz(-0.73211615) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1248963) q[2];
sx q[2];
rz(-2.5113528) q[2];
sx q[2];
rz(-1.5843887) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.447532) q[1];
sx q[1];
rz(-0.58528712) q[1];
sx q[1];
rz(-1.8156169) q[1];
rz(-1.2110662) q[3];
sx q[3];
rz(-1.779428) q[3];
sx q[3];
rz(1.9344575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6901107) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(3.0991128) q[2];
rz(0.63043198) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38934389) q[0];
sx q[0];
rz(-1.2279953) q[0];
sx q[0];
rz(-0.4831627) q[0];
rz(-1.0522316) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(-2.5767456) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72950596) q[0];
sx q[0];
rz(-2.5967263) q[0];
sx q[0];
rz(1.8338404) q[0];
rz(-1.9916612) q[2];
sx q[2];
rz(-1.5016342) q[2];
sx q[2];
rz(-1.3644621) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40183345) q[1];
sx q[1];
rz(-1.9976915) q[1];
sx q[1];
rz(3.0763807) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4156028) q[3];
sx q[3];
rz(-1.3421913) q[3];
sx q[3];
rz(0.58333635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.57006449) q[2];
sx q[2];
rz(-2.0740985) q[2];
sx q[2];
rz(-1.8035536) q[2];
rz(1.8367052) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(2.4664972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9682482) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(-2.5937953) q[0];
rz(0.785218) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(0.26842591) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7086605) q[0];
sx q[0];
rz(-1.4888568) q[0];
sx q[0];
rz(-2.8952778) q[0];
x q[1];
rz(-2.354291) q[2];
sx q[2];
rz(-2.1848218) q[2];
sx q[2];
rz(1.8067443) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5749579) q[1];
sx q[1];
rz(-2.1556427) q[1];
sx q[1];
rz(3.023874) q[1];
rz(-pi) q[2];
rz(-2.8397452) q[3];
sx q[3];
rz(-1.0511304) q[3];
sx q[3];
rz(2.2083851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61775529) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(-0.8141554) q[2];
rz(2.7653149) q[3];
sx q[3];
rz(-1.9777931) q[3];
sx q[3];
rz(-3.0686839) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5557264) q[0];
sx q[0];
rz(-1.7365475) q[0];
sx q[0];
rz(1.0193753) q[0];
rz(2.2881919) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(-2.6928435) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42496029) q[0];
sx q[0];
rz(-2.6052931) q[0];
sx q[0];
rz(2.0400356) q[0];
x q[1];
rz(2.8837187) q[2];
sx q[2];
rz(-1.2784625) q[2];
sx q[2];
rz(2.984798) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6905578) q[1];
sx q[1];
rz(-2.096855) q[1];
sx q[1];
rz(-1.2760389) q[1];
x q[2];
rz(-0.71172165) q[3];
sx q[3];
rz(-0.91152836) q[3];
sx q[3];
rz(3.1275415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.86924187) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(-0.31420079) q[2];
rz(2.3172486) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(-2.3468988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070351275) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(1.8810133) q[0];
rz(-0.64385995) q[1];
sx q[1];
rz(-1.9088129) q[1];
sx q[1];
rz(0.018928122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8185794) q[0];
sx q[0];
rz(-1.5886663) q[0];
sx q[0];
rz(-0.26364003) q[0];
rz(-pi) q[1];
rz(0.55982121) q[2];
sx q[2];
rz(-2.9580742) q[2];
sx q[2];
rz(0.57932094) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70871204) q[1];
sx q[1];
rz(-1.6212665) q[1];
sx q[1];
rz(-2.9570079) q[1];
x q[2];
rz(-0.67846672) q[3];
sx q[3];
rz(-2.7760091) q[3];
sx q[3];
rz(1.4640704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54946047) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(-2.4712759) q[2];
rz(-2.629225) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(1.4204773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5877514) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(-2.642139) q[0];
rz(-1.5746501) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(1.0826899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4944045) q[0];
sx q[0];
rz(-1.5218966) q[0];
sx q[0];
rz(-0.58736579) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.322299) q[2];
sx q[2];
rz(-1.5466006) q[2];
sx q[2];
rz(2.1404612) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.078552695) q[1];
sx q[1];
rz(-1.551997) q[1];
sx q[1];
rz(-2.5979554) q[1];
x q[2];
rz(-2.4389078) q[3];
sx q[3];
rz(-2.6583238) q[3];
sx q[3];
rz(-1.0815222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3988951) q[2];
sx q[2];
rz(-1.164914) q[2];
sx q[2];
rz(-2.6297074) q[2];
rz(0.39294696) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(2.0846562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95505161) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(-2.4004249) q[1];
sx q[1];
rz(-0.82294958) q[1];
sx q[1];
rz(-0.23946147) q[1];
rz(-0.95460931) q[2];
sx q[2];
rz(-1.90345) q[2];
sx q[2];
rz(-0.30751139) q[2];
rz(-0.74647222) q[3];
sx q[3];
rz(-1.6975879) q[3];
sx q[3];
rz(-3.0293037) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
