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
rz(-2.9838188) q[0];
sx q[0];
rz(4.3133419) q[0];
sx q[0];
rz(11.316909) q[0];
rz(-3.2818031) q[1];
sx q[1];
rz(1.6970716) q[1];
sx q[1];
rz(12.628218) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0573579) q[0];
sx q[0];
rz(-2.0664729) q[0];
sx q[0];
rz(2.2593014) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0038764406) q[2];
sx q[2];
rz(-3.0511694) q[2];
sx q[2];
rz(0.10744444) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7951491) q[1];
sx q[1];
rz(-0.7278053) q[1];
sx q[1];
rz(1.4239356) q[1];
rz(-0.40932406) q[3];
sx q[3];
rz(-2.217389) q[3];
sx q[3];
rz(-2.4149946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42741117) q[2];
sx q[2];
rz(-0.84250557) q[2];
sx q[2];
rz(0.25644914) q[2];
rz(-2.9668258) q[3];
sx q[3];
rz(-1.1656961) q[3];
sx q[3];
rz(-0.83785653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8137708) q[0];
sx q[0];
rz(-2.796266) q[0];
sx q[0];
rz(0.12369618) q[0];
rz(0.23873121) q[1];
sx q[1];
rz(-2.3614466) q[1];
sx q[1];
rz(3.0366268) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83451) q[0];
sx q[0];
rz(-1.7210467) q[0];
sx q[0];
rz(-0.86850496) q[0];
rz(-pi) q[1];
rz(-1.3388073) q[2];
sx q[2];
rz(-1.9769443) q[2];
sx q[2];
rz(-2.7628203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.87890307) q[1];
sx q[1];
rz(-0.57003747) q[1];
sx q[1];
rz(-2.0712584) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41173287) q[3];
sx q[3];
rz(-0.74854088) q[3];
sx q[3];
rz(-0.68621503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77439848) q[2];
sx q[2];
rz(-1.6509193) q[2];
sx q[2];
rz(1.6801838) q[2];
rz(-1.2857619) q[3];
sx q[3];
rz(-1.6971842) q[3];
sx q[3];
rz(2.8694966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1478145) q[0];
sx q[0];
rz(-0.21037924) q[0];
sx q[0];
rz(-1.0149957) q[0];
rz(2.2442832) q[1];
sx q[1];
rz(-1.6474479) q[1];
sx q[1];
rz(2.4651333) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2059442) q[0];
sx q[0];
rz(-2.1876011) q[0];
sx q[0];
rz(-0.031868462) q[0];
rz(-pi) q[1];
rz(-2.1968165) q[2];
sx q[2];
rz(-1.9817366) q[2];
sx q[2];
rz(2.0650244) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8415422) q[1];
sx q[1];
rz(-1.1917229) q[1];
sx q[1];
rz(2.9247523) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5514938) q[3];
sx q[3];
rz(-0.75452828) q[3];
sx q[3];
rz(-1.7280662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.93495381) q[2];
sx q[2];
rz(-2.166344) q[2];
sx q[2];
rz(-2.3812531) q[2];
rz(1.9350516) q[3];
sx q[3];
rz(-2.3307266) q[3];
sx q[3];
rz(0.84347239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.7059785) q[0];
sx q[0];
rz(-2.5137081) q[0];
sx q[0];
rz(1.5465558) q[0];
rz(0.71890038) q[1];
sx q[1];
rz(-1.8214106) q[1];
sx q[1];
rz(0.23153201) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084615413) q[0];
sx q[0];
rz(-1.9997134) q[0];
sx q[0];
rz(0.57022988) q[0];
rz(-pi) q[1];
rz(1.5668554) q[2];
sx q[2];
rz(-0.77814279) q[2];
sx q[2];
rz(-0.63150327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12688237) q[1];
sx q[1];
rz(-2.0538804) q[1];
sx q[1];
rz(1.4654069) q[1];
x q[2];
rz(-2.6014884) q[3];
sx q[3];
rz(-2.4632235) q[3];
sx q[3];
rz(-2.2578007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.71451688) q[2];
sx q[2];
rz(-0.38893739) q[2];
sx q[2];
rz(2.5330949) q[2];
rz(0.19615873) q[3];
sx q[3];
rz(-1.4213296) q[3];
sx q[3];
rz(2.1430446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.8683559) q[0];
sx q[0];
rz(-2.4899857) q[0];
sx q[0];
rz(2.3625145) q[0];
rz(0.92998663) q[1];
sx q[1];
rz(-1.9735186) q[1];
sx q[1];
rz(1.9680061) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1176276) q[0];
sx q[0];
rz(-1.574734) q[0];
sx q[0];
rz(-1.9088452) q[0];
rz(-3.0069238) q[2];
sx q[2];
rz(-0.95974082) q[2];
sx q[2];
rz(-0.5765644) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0231578) q[1];
sx q[1];
rz(-1.506532) q[1];
sx q[1];
rz(-1.8885018) q[1];
rz(-pi) q[2];
rz(1.7608123) q[3];
sx q[3];
rz(-1.992986) q[3];
sx q[3];
rz(-1.9334396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.32213) q[2];
sx q[2];
rz(-0.16699114) q[2];
sx q[2];
rz(2.4665311) q[2];
rz(2.2760462) q[3];
sx q[3];
rz(-1.6277438) q[3];
sx q[3];
rz(1.2077829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.753767) q[0];
sx q[0];
rz(-3.1238811) q[0];
sx q[0];
rz(-0.87131635) q[0];
rz(-0.21698347) q[1];
sx q[1];
rz(-1.5996409) q[1];
sx q[1];
rz(1.3380231) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7965391) q[0];
sx q[0];
rz(-0.65700475) q[0];
sx q[0];
rz(3.0057602) q[0];
x q[1];
rz(-2.2964444) q[2];
sx q[2];
rz(-1.6763699) q[2];
sx q[2];
rz(-1.1136063) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.034778193) q[1];
sx q[1];
rz(-0.15131525) q[1];
sx q[1];
rz(-1.8968015) q[1];
x q[2];
rz(-1.3806255) q[3];
sx q[3];
rz(-2.0479879) q[3];
sx q[3];
rz(-0.3482477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1285105) q[2];
sx q[2];
rz(-1.6435813) q[2];
sx q[2];
rz(2.3415372) q[2];
rz(-0.99177805) q[3];
sx q[3];
rz(-0.9477152) q[3];
sx q[3];
rz(0.94314027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25471383) q[0];
sx q[0];
rz(-2.7006221) q[0];
sx q[0];
rz(-2.9898341) q[0];
rz(1.4166547) q[1];
sx q[1];
rz(-2.2817426) q[1];
sx q[1];
rz(-0.83622611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1372414) q[0];
sx q[0];
rz(-0.81295952) q[0];
sx q[0];
rz(-1.0197958) q[0];
rz(-pi) q[1];
rz(-2.5649037) q[2];
sx q[2];
rz(-1.5778013) q[2];
sx q[2];
rz(3.0617065) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5150093) q[1];
sx q[1];
rz(-1.8191511) q[1];
sx q[1];
rz(3.0529313) q[1];
x q[2];
rz(-0.25612513) q[3];
sx q[3];
rz(-0.47893347) q[3];
sx q[3];
rz(2.0253945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4862711) q[2];
sx q[2];
rz(-0.063491193) q[2];
sx q[2];
rz(3.0625694) q[2];
rz(-0.71845734) q[3];
sx q[3];
rz(-1.6301165) q[3];
sx q[3];
rz(-0.92459905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35543168) q[0];
sx q[0];
rz(-2.5434255) q[0];
sx q[0];
rz(-0.86791903) q[0];
rz(-0.68880853) q[1];
sx q[1];
rz(-2.8354366) q[1];
sx q[1];
rz(0.93592962) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3132738) q[0];
sx q[0];
rz(-0.5567282) q[0];
sx q[0];
rz(0.3484881) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8008119) q[2];
sx q[2];
rz(-1.9592957) q[2];
sx q[2];
rz(2.4430371) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3080146) q[1];
sx q[1];
rz(-1.2460099) q[1];
sx q[1];
rz(-0.056450162) q[1];
x q[2];
rz(0.0015663107) q[3];
sx q[3];
rz(-0.74347444) q[3];
sx q[3];
rz(1.2921799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1615289) q[2];
sx q[2];
rz(-0.65902013) q[2];
sx q[2];
rz(2.8847983) q[2];
rz(-1.0181381) q[3];
sx q[3];
rz(-2.0223821) q[3];
sx q[3];
rz(-2.170678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7629906) q[0];
sx q[0];
rz(-2.584223) q[0];
sx q[0];
rz(0.85365224) q[0];
rz(-1.8866106) q[1];
sx q[1];
rz(-1.9957142) q[1];
sx q[1];
rz(-0.34057158) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69215323) q[0];
sx q[0];
rz(-1.3763577) q[0];
sx q[0];
rz(-2.9783335) q[0];
rz(-pi) q[1];
rz(2.8355424) q[2];
sx q[2];
rz(-0.69717625) q[2];
sx q[2];
rz(-2.9070284) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3601968) q[1];
sx q[1];
rz(-2.2193647) q[1];
sx q[1];
rz(2.1482723) q[1];
rz(-pi) q[2];
rz(0.27517516) q[3];
sx q[3];
rz(-2.0180297) q[3];
sx q[3];
rz(-2.3768611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7178932) q[2];
sx q[2];
rz(-1.3807978) q[2];
sx q[2];
rz(0.37810668) q[2];
rz(0.2374436) q[3];
sx q[3];
rz(-2.0707371) q[3];
sx q[3];
rz(-0.28089359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045192748) q[0];
sx q[0];
rz(-1.0843596) q[0];
sx q[0];
rz(2.4943446) q[0];
rz(-1.9735533) q[1];
sx q[1];
rz(-2.2868575) q[1];
sx q[1];
rz(0.38356575) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3733394) q[0];
sx q[0];
rz(-0.037469286) q[0];
sx q[0];
rz(-0.61385198) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9012544) q[2];
sx q[2];
rz(-2.0784194) q[2];
sx q[2];
rz(2.3445545) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4216363) q[1];
sx q[1];
rz(-1.7767748) q[1];
sx q[1];
rz(0.27011807) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0553352) q[3];
sx q[3];
rz(-1.4773507) q[3];
sx q[3];
rz(-2.6938113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26455227) q[2];
sx q[2];
rz(-2.2525747) q[2];
sx q[2];
rz(1.5569347) q[2];
rz(-1.6282188) q[3];
sx q[3];
rz(-1.7729365) q[3];
sx q[3];
rz(-2.7148066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3601111) q[0];
sx q[0];
rz(-1.7088912) q[0];
sx q[0];
rz(0.29722469) q[0];
rz(1.1813286) q[1];
sx q[1];
rz(-1.344463) q[1];
sx q[1];
rz(-0.32483473) q[1];
rz(2.8988373) q[2];
sx q[2];
rz(-2.0363672) q[2];
sx q[2];
rz(-2.80574) q[2];
rz(-1.6024797) q[3];
sx q[3];
rz(-1.0864232) q[3];
sx q[3];
rz(-0.3530799) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
