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
rz(-2.7741127) q[0];
sx q[0];
rz(-1.8125266) q[0];
sx q[0];
rz(-1.4867866) q[0];
rz(-1.6410671) q[1];
sx q[1];
rz(-0.60368109) q[1];
sx q[1];
rz(0.85890213) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6971815) q[0];
sx q[0];
rz(-1.0134122) q[0];
sx q[0];
rz(-1.0943221) q[0];
x q[1];
rz(-2.3043465) q[2];
sx q[2];
rz(-2.4105802) q[2];
sx q[2];
rz(1.8164509) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.907469) q[1];
sx q[1];
rz(-0.69141839) q[1];
sx q[1];
rz(-1.3983634) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1726491) q[3];
sx q[3];
rz(-2.7730983) q[3];
sx q[3];
rz(-0.70337765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.19771244) q[2];
sx q[2];
rz(-2.3494425) q[2];
sx q[2];
rz(2.0749178) q[2];
rz(-3.0471622) q[3];
sx q[3];
rz(-2.5240199) q[3];
sx q[3];
rz(0.42816952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18148947) q[0];
sx q[0];
rz(-0.26679376) q[0];
sx q[0];
rz(-1.9285404) q[0];
rz(0.45252291) q[1];
sx q[1];
rz(-2.8892543) q[1];
sx q[1];
rz(2.4620893) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0328355) q[0];
sx q[0];
rz(-1.8201314) q[0];
sx q[0];
rz(1.6052206) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29085906) q[2];
sx q[2];
rz(-1.4122987) q[2];
sx q[2];
rz(3.045451) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7391239) q[1];
sx q[1];
rz(-1.1796412) q[1];
sx q[1];
rz(0.7372589) q[1];
rz(-pi) q[2];
rz(2.4538058) q[3];
sx q[3];
rz(-1.4151598) q[3];
sx q[3];
rz(-0.17222541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.54000336) q[2];
sx q[2];
rz(-0.82401472) q[2];
sx q[2];
rz(2.4448815) q[2];
rz(-1.2434897) q[3];
sx q[3];
rz(-1.1949298) q[3];
sx q[3];
rz(0.5881601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531279) q[0];
sx q[0];
rz(-1.0364113) q[0];
sx q[0];
rz(2.3179407) q[0];
rz(2.9846687) q[1];
sx q[1];
rz(-1.7053968) q[1];
sx q[1];
rz(-2.7395111) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7347468) q[0];
sx q[0];
rz(-1.1752793) q[0];
sx q[0];
rz(2.9989373) q[0];
rz(-pi) q[1];
rz(-1.4060611) q[2];
sx q[2];
rz(-2.852612) q[2];
sx q[2];
rz(-1.8584205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3797219) q[1];
sx q[1];
rz(-2.5609697) q[1];
sx q[1];
rz(-0.94731829) q[1];
x q[2];
rz(1.9605432) q[3];
sx q[3];
rz(-1.7754835) q[3];
sx q[3];
rz(-1.2677416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1261525) q[2];
sx q[2];
rz(-0.93706477) q[2];
sx q[2];
rz(0.13119571) q[2];
rz(-2.0845856) q[3];
sx q[3];
rz(-1.1934692) q[3];
sx q[3];
rz(-1.2098275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(2.7159202) q[0];
sx q[0];
rz(-1.965006) q[0];
sx q[0];
rz(-1.7536989) q[0];
rz(-1.362494) q[1];
sx q[1];
rz(-1.8905996) q[1];
sx q[1];
rz(-1.9349792) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0232892) q[0];
sx q[0];
rz(-2.1739066) q[0];
sx q[0];
rz(1.009737) q[0];
rz(0.64995857) q[2];
sx q[2];
rz(-2.0043249) q[2];
sx q[2];
rz(-0.94160801) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6304508) q[1];
sx q[1];
rz(-1.8725553) q[1];
sx q[1];
rz(-1.6514227) q[1];
rz(-pi) q[2];
rz(-2.4033) q[3];
sx q[3];
rz(-1.2311934) q[3];
sx q[3];
rz(2.5552354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6917307) q[2];
sx q[2];
rz(-1.3351771) q[2];
sx q[2];
rz(2.195669) q[2];
rz(2.6876884) q[3];
sx q[3];
rz(-2.4920521) q[3];
sx q[3];
rz(0.18979931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40189633) q[0];
sx q[0];
rz(-1.3119768) q[0];
sx q[0];
rz(0.78347462) q[0];
rz(-0.90617323) q[1];
sx q[1];
rz(-2.2515191) q[1];
sx q[1];
rz(2.5772212) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2191129) q[0];
sx q[0];
rz(-1.4267491) q[0];
sx q[0];
rz(-3.0213468) q[0];
rz(-pi) q[1];
rz(1.7092114) q[2];
sx q[2];
rz(-0.40988806) q[2];
sx q[2];
rz(-2.5382535) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42309883) q[1];
sx q[1];
rz(-2.1717487) q[1];
sx q[1];
rz(-2.3312631) q[1];
rz(2.4268812) q[3];
sx q[3];
rz(-1.6922975) q[3];
sx q[3];
rz(-0.35706676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.924661) q[2];
sx q[2];
rz(-1.4035808) q[2];
sx q[2];
rz(-3.0972287) q[2];
rz(1.4102604) q[3];
sx q[3];
rz(-2.9328465) q[3];
sx q[3];
rz(1.1948208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.57589543) q[0];
sx q[0];
rz(-0.79224753) q[0];
sx q[0];
rz(0.76882291) q[0];
rz(2.5104751) q[1];
sx q[1];
rz(-1.2080071) q[1];
sx q[1];
rz(2.9879976) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0868743) q[0];
sx q[0];
rz(-1.2445893) q[0];
sx q[0];
rz(2.1534462) q[0];
rz(1.8472438) q[2];
sx q[2];
rz(-2.203456) q[2];
sx q[2];
rz(-1.1142434) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99912737) q[1];
sx q[1];
rz(-2.8395238) q[1];
sx q[1];
rz(2.7727454) q[1];
x q[2];
rz(0.77135857) q[3];
sx q[3];
rz(-1.5017209) q[3];
sx q[3];
rz(1.3265287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.067387335) q[2];
sx q[2];
rz(-1.5978483) q[2];
sx q[2];
rz(0.81573168) q[2];
rz(0.054281209) q[3];
sx q[3];
rz(-1.7547601) q[3];
sx q[3];
rz(-1.3666216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62189198) q[0];
sx q[0];
rz(-1.0655572) q[0];
sx q[0];
rz(-0.52072293) q[0];
rz(-2.0479274) q[1];
sx q[1];
rz(-0.92898527) q[1];
sx q[1];
rz(1.7589794) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0313697) q[0];
sx q[0];
rz(-0.58833226) q[0];
sx q[0];
rz(-1.3084685) q[0];
x q[1];
rz(-0.52806921) q[2];
sx q[2];
rz(-1.1748399) q[2];
sx q[2];
rz(-1.2483526) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7944889) q[1];
sx q[1];
rz(-2.1784049) q[1];
sx q[1];
rz(-2.6332864) q[1];
x q[2];
rz(1.3741854) q[3];
sx q[3];
rz(-1.1093321) q[3];
sx q[3];
rz(-2.9795071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19892056) q[2];
sx q[2];
rz(-0.3519381) q[2];
sx q[2];
rz(1.675763) q[2];
rz(0.28907019) q[3];
sx q[3];
rz(-1.7710268) q[3];
sx q[3];
rz(-1.2808965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11956231) q[0];
sx q[0];
rz(-0.95198315) q[0];
sx q[0];
rz(-0.49222487) q[0];
rz(0.64552632) q[1];
sx q[1];
rz(-1.5954834) q[1];
sx q[1];
rz(0.92799497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6715901) q[0];
sx q[0];
rz(-0.3151463) q[0];
sx q[0];
rz(-1.1130087) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83023425) q[2];
sx q[2];
rz(-0.50721675) q[2];
sx q[2];
rz(-0.00032256034) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.38233373) q[1];
sx q[1];
rz(-1.302726) q[1];
sx q[1];
rz(-2.2925582) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2755166) q[3];
sx q[3];
rz(-0.47086916) q[3];
sx q[3];
rz(2.4195723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8202028) q[2];
sx q[2];
rz(-1.5265744) q[2];
sx q[2];
rz(2.3825633) q[2];
rz(-1.8697033) q[3];
sx q[3];
rz(-1.8153518) q[3];
sx q[3];
rz(2.429764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7976545) q[0];
sx q[0];
rz(-0.54867083) q[0];
sx q[0];
rz(-0.77242533) q[0];
rz(1.7732357) q[1];
sx q[1];
rz(-1.1027579) q[1];
sx q[1];
rz(0.30430421) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40201449) q[0];
sx q[0];
rz(-2.9106044) q[0];
sx q[0];
rz(-0.84257479) q[0];
rz(3.1379406) q[2];
sx q[2];
rz(-1.7965791) q[2];
sx q[2];
rz(-0.89930764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9134679) q[1];
sx q[1];
rz(-2.2992837) q[1];
sx q[1];
rz(-1.8124652) q[1];
x q[2];
rz(-2.8465495) q[3];
sx q[3];
rz(-2.7310762) q[3];
sx q[3];
rz(-0.66995588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7336537) q[2];
sx q[2];
rz(-2.759178) q[2];
sx q[2];
rz(-1.9668874) q[2];
rz(0.63742739) q[3];
sx q[3];
rz(-1.8791684) q[3];
sx q[3];
rz(1.2229961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2735073) q[0];
sx q[0];
rz(-0.33184505) q[0];
sx q[0];
rz(-2.3688431) q[0];
rz(-2.6840774) q[1];
sx q[1];
rz(-1.7465778) q[1];
sx q[1];
rz(-1.4332019) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3210216) q[0];
sx q[0];
rz(-0.95527705) q[0];
sx q[0];
rz(-2.5825898) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3914621) q[2];
sx q[2];
rz(-2.7731189) q[2];
sx q[2];
rz(2.3725906) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0876574) q[1];
sx q[1];
rz(-1.1217351) q[1];
sx q[1];
rz(-0.29718558) q[1];
rz(-pi) q[2];
rz(2.1674311) q[3];
sx q[3];
rz(-2.0630368) q[3];
sx q[3];
rz(-2.8288768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92051238) q[2];
sx q[2];
rz(-0.94197333) q[2];
sx q[2];
rz(2.4648049) q[2];
rz(1.5306728) q[3];
sx q[3];
rz(-2.8407606) q[3];
sx q[3];
rz(-1.0845186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21657011) q[0];
sx q[0];
rz(-1.3022447) q[0];
sx q[0];
rz(-1.362823) q[0];
rz(2.2577747) q[1];
sx q[1];
rz(-0.73278905) q[1];
sx q[1];
rz(-0.97272452) q[1];
rz(-0.92166273) q[2];
sx q[2];
rz(-2.4646152) q[2];
sx q[2];
rz(1.2645417) q[2];
rz(-1.3317713) q[3];
sx q[3];
rz(-1.1186794) q[3];
sx q[3];
rz(0.787822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
