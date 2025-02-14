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
rz(-1.9698434) q[0];
sx q[0];
rz(1.8921312) q[0];
rz(-3.2818031) q[1];
sx q[1];
rz(1.6970716) q[1];
sx q[1];
rz(12.628218) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0573579) q[0];
sx q[0];
rz(-1.0751197) q[0];
sx q[0];
rz(-2.2593014) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5704448) q[2];
sx q[2];
rz(-1.6612189) q[2];
sx q[2];
rz(3.0302559) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0272408) q[1];
sx q[1];
rz(-1.6682965) q[1];
sx q[1];
rz(-2.2932294) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0855882) q[3];
sx q[3];
rz(-0.74923979) q[3];
sx q[3];
rz(3.0390714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42741117) q[2];
sx q[2];
rz(-0.84250557) q[2];
sx q[2];
rz(-0.25644914) q[2];
rz(-0.17476684) q[3];
sx q[3];
rz(-1.1656961) q[3];
sx q[3];
rz(-2.3037361) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8137708) q[0];
sx q[0];
rz(-2.796266) q[0];
sx q[0];
rz(0.12369618) q[0];
rz(-0.23873121) q[1];
sx q[1];
rz(-2.3614466) q[1];
sx q[1];
rz(-3.0366268) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86228131) q[0];
sx q[0];
rz(-0.87800069) q[0];
sx q[0];
rz(-0.19578085) q[0];
rz(-pi) q[1];
rz(-1.8027854) q[2];
sx q[2];
rz(-1.1646484) q[2];
sx q[2];
rz(-2.7628203) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1233842) q[1];
sx q[1];
rz(-1.3088641) q[1];
sx q[1];
rz(2.0831152) q[1];
rz(-pi) q[2];
rz(-0.70521783) q[3];
sx q[3];
rz(-1.2949484) q[3];
sx q[3];
rz(-0.57491702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77439848) q[2];
sx q[2];
rz(-1.6509193) q[2];
sx q[2];
rz(-1.4614089) q[2];
rz(-1.8558308) q[3];
sx q[3];
rz(-1.6971842) q[3];
sx q[3];
rz(0.27209601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1478145) q[0];
sx q[0];
rz(-2.9312134) q[0];
sx q[0];
rz(-2.126597) q[0];
rz(0.89730942) q[1];
sx q[1];
rz(-1.6474479) q[1];
sx q[1];
rz(-2.4651333) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9356485) q[0];
sx q[0];
rz(-2.1876011) q[0];
sx q[0];
rz(-0.031868462) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2102506) q[2];
sx q[2];
rz(-2.4081875) q[2];
sx q[2];
rz(-2.1426107) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8415422) q[1];
sx q[1];
rz(-1.1917229) q[1];
sx q[1];
rz(-0.21684034) q[1];
rz(-2.3252316) q[3];
sx q[3];
rz(-1.5575756) q[3];
sx q[3];
rz(0.17133443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93495381) q[2];
sx q[2];
rz(-0.97524869) q[2];
sx q[2];
rz(2.3812531) q[2];
rz(-1.9350516) q[3];
sx q[3];
rz(-2.3307266) q[3];
sx q[3];
rz(2.2981203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43561414) q[0];
sx q[0];
rz(-0.62788457) q[0];
sx q[0];
rz(1.5950369) q[0];
rz(0.71890038) q[1];
sx q[1];
rz(-1.3201821) q[1];
sx q[1];
rz(-0.23153201) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2255428) q[0];
sx q[0];
rz(-2.0839491) q[0];
sx q[0];
rz(1.0731359) q[0];
rz(-pi) q[1];
rz(-0.79265742) q[2];
sx q[2];
rz(-1.5680299) q[2];
sx q[2];
rz(-0.94209988) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.097024767) q[1];
sx q[1];
rz(-2.6480354) q[1];
sx q[1];
rz(-2.943671) q[1];
rz(-2.6014884) q[3];
sx q[3];
rz(-0.67836919) q[3];
sx q[3];
rz(2.2578007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.71451688) q[2];
sx q[2];
rz(-2.7526553) q[2];
sx q[2];
rz(2.5330949) q[2];
rz(2.9454339) q[3];
sx q[3];
rz(-1.4213296) q[3];
sx q[3];
rz(0.99854809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.8683559) q[0];
sx q[0];
rz(-2.4899857) q[0];
sx q[0];
rz(0.77907816) q[0];
rz(-0.92998663) q[1];
sx q[1];
rz(-1.168074) q[1];
sx q[1];
rz(-1.1735865) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54544696) q[0];
sx q[0];
rz(-1.2327502) q[0];
sx q[0];
rz(0.0041739504) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9554602) q[2];
sx q[2];
rz(-1.6809856) q[2];
sx q[2];
rz(-0.9166536) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25949892) q[1];
sx q[1];
rz(-2.8176687) q[1];
sx q[1];
rz(-1.367635) q[1];
x q[2];
rz(0.39802528) q[3];
sx q[3];
rz(-0.46062352) q[3];
sx q[3];
rz(-1.4945791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.32213) q[2];
sx q[2];
rz(-2.9746015) q[2];
sx q[2];
rz(-2.4665311) q[2];
rz(2.2760462) q[3];
sx q[3];
rz(-1.5138488) q[3];
sx q[3];
rz(-1.2077829) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38782564) q[0];
sx q[0];
rz(-3.1238811) q[0];
sx q[0];
rz(0.87131635) q[0];
rz(-0.21698347) q[1];
sx q[1];
rz(-1.5419518) q[1];
sx q[1];
rz(1.8035696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8080475) q[0];
sx q[0];
rz(-1.6535954) q[0];
sx q[0];
rz(0.6525349) q[0];
x q[1];
rz(0.84514825) q[2];
sx q[2];
rz(-1.6763699) q[2];
sx q[2];
rz(2.0279864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1068145) q[1];
sx q[1];
rz(-0.15131525) q[1];
sx q[1];
rz(-1.2447912) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3806255) q[3];
sx q[3];
rz(-1.0936048) q[3];
sx q[3];
rz(0.3482477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1285105) q[2];
sx q[2];
rz(-1.4980114) q[2];
sx q[2];
rz(2.3415372) q[2];
rz(-2.1498146) q[3];
sx q[3];
rz(-0.9477152) q[3];
sx q[3];
rz(-0.94314027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25471383) q[0];
sx q[0];
rz(-2.7006221) q[0];
sx q[0];
rz(-0.15175858) q[0];
rz(-1.4166547) q[1];
sx q[1];
rz(-2.2817426) q[1];
sx q[1];
rz(0.83622611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4164616) q[0];
sx q[0];
rz(-2.2380479) q[0];
sx q[0];
rz(-2.6362772) q[0];
rz(-3.128746) q[2];
sx q[2];
rz(-2.564866) q[2];
sx q[2];
rz(-1.6614514) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.077634634) q[1];
sx q[1];
rz(-1.4848621) q[1];
sx q[1];
rz(-1.3215022) q[1];
rz(1.4400021) q[3];
sx q[3];
rz(-2.0328641) q[3];
sx q[3];
rz(-0.82926428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6553216) q[2];
sx q[2];
rz(-3.0781015) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.786161) q[0];
sx q[0];
rz(-0.59816718) q[0];
sx q[0];
rz(0.86791903) q[0];
rz(2.4527841) q[1];
sx q[1];
rz(-2.8354366) q[1];
sx q[1];
rz(-2.205663) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42396069) q[0];
sx q[0];
rz(-2.0905295) q[0];
sx q[0];
rz(-1.7802159) q[0];
rz(-pi) q[1];
rz(2.8008119) q[2];
sx q[2];
rz(-1.182297) q[2];
sx q[2];
rz(-0.69855554) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8607801) q[1];
sx q[1];
rz(-1.6242923) q[1];
sx q[1];
rz(-1.2455275) q[1];
rz(-1.5722365) q[3];
sx q[3];
rz(-0.82732302) q[3];
sx q[3];
rz(1.2900521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1615289) q[2];
sx q[2];
rz(-0.65902013) q[2];
sx q[2];
rz(0.25679437) q[2];
rz(-1.0181381) q[3];
sx q[3];
rz(-2.0223821) q[3];
sx q[3];
rz(-2.170678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7629906) q[0];
sx q[0];
rz(-2.584223) q[0];
sx q[0];
rz(0.85365224) q[0];
rz(1.8866106) q[1];
sx q[1];
rz(-1.9957142) q[1];
sx q[1];
rz(0.34057158) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4494394) q[0];
sx q[0];
rz(-1.3763577) q[0];
sx q[0];
rz(2.9783335) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8179632) q[2];
sx q[2];
rz(-0.91193141) q[2];
sx q[2];
rz(2.5160922) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6046227) q[1];
sx q[1];
rz(-2.3019362) q[1];
sx q[1];
rz(0.62418681) q[1];
rz(2.8664175) q[3];
sx q[3];
rz(-1.1235629) q[3];
sx q[3];
rz(0.76473151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42369947) q[2];
sx q[2];
rz(-1.3807978) q[2];
sx q[2];
rz(2.763486) q[2];
rz(-2.9041491) q[3];
sx q[3];
rz(-1.0708555) q[3];
sx q[3];
rz(0.28089359) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045192748) q[0];
sx q[0];
rz(-2.0572331) q[0];
sx q[0];
rz(-2.4943446) q[0];
rz(1.1680394) q[1];
sx q[1];
rz(-0.85473514) q[1];
sx q[1];
rz(-0.38356575) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7591568) q[0];
sx q[0];
rz(-1.6014227) q[0];
sx q[0];
rz(1.5492065) q[0];
rz(-0.53159376) q[2];
sx q[2];
rz(-1.2832912) q[2];
sx q[2];
rz(2.2026041) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.79427045) q[1];
sx q[1];
rz(-1.8350661) q[1];
sx q[1];
rz(1.7842954) q[1];
rz(-pi) q[2];
rz(-1.7586772) q[3];
sx q[3];
rz(-2.6184824) q[3];
sx q[3];
rz(2.1818102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26455227) q[2];
sx q[2];
rz(-0.88901797) q[2];
sx q[2];
rz(1.5846579) q[2];
rz(-1.6282188) q[3];
sx q[3];
rz(-1.3686562) q[3];
sx q[3];
rz(2.7148066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7814816) q[0];
sx q[0];
rz(-1.7088912) q[0];
sx q[0];
rz(0.29722469) q[0];
rz(1.1813286) q[1];
sx q[1];
rz(-1.344463) q[1];
sx q[1];
rz(-0.32483473) q[1];
rz(1.0931821) q[2];
sx q[2];
rz(-1.3543159) q[2];
sx q[2];
rz(-1.3456624) q[2];
rz(3.0814617) q[3];
sx q[3];
rz(-2.6562666) q[3];
sx q[3];
rz(-0.2851214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
