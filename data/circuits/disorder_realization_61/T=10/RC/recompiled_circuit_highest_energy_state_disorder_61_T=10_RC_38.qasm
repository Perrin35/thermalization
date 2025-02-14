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
rz(-1.6725809) q[0];
sx q[0];
rz(-2.2218158) q[0];
sx q[0];
rz(2.4866009) q[0];
rz(-1.4870149) q[1];
sx q[1];
rz(-1.4227285) q[1];
sx q[1];
rz(1.3230327) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9441863) q[0];
sx q[0];
rz(-1.8776769) q[0];
sx q[0];
rz(2.5635864) q[0];
rz(1.8048067) q[2];
sx q[2];
rz(-1.26097) q[2];
sx q[2];
rz(1.3739623) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6055687) q[1];
sx q[1];
rz(-1.6704511) q[1];
sx q[1];
rz(-1.8479363) q[1];
rz(-0.58889525) q[3];
sx q[3];
rz(-2.3213904) q[3];
sx q[3];
rz(-0.39862788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.46372867) q[2];
sx q[2];
rz(-2.4552796) q[2];
sx q[2];
rz(2.4832671) q[2];
rz(2.5074734) q[3];
sx q[3];
rz(-1.442797) q[3];
sx q[3];
rz(1.770795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.051801) q[0];
sx q[0];
rz(-0.3312411) q[0];
sx q[0];
rz(-3.064503) q[0];
rz(0.58473051) q[1];
sx q[1];
rz(-1.3405776) q[1];
sx q[1];
rz(-1.1999406) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8611885) q[0];
sx q[0];
rz(-1.1139835) q[0];
sx q[0];
rz(2.9218319) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4064155) q[2];
sx q[2];
rz(-1.5417678) q[2];
sx q[2];
rz(2.8433702) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.40545826) q[1];
sx q[1];
rz(-1.9455487) q[1];
sx q[1];
rz(2.2103851) q[1];
x q[2];
rz(-2.582483) q[3];
sx q[3];
rz(-1.5190795) q[3];
sx q[3];
rz(-2.8962108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8889019) q[2];
sx q[2];
rz(-2.1609047) q[2];
sx q[2];
rz(-2.4661031) q[2];
rz(3.1318956) q[3];
sx q[3];
rz(-2.8986425) q[3];
sx q[3];
rz(-3.1233136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49552396) q[0];
sx q[0];
rz(-1.4718453) q[0];
sx q[0];
rz(-0.75463265) q[0];
rz(3.1309639) q[1];
sx q[1];
rz(-1.3898712) q[1];
sx q[1];
rz(-2.1307814) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36463144) q[0];
sx q[0];
rz(-1.8753795) q[0];
sx q[0];
rz(1.0582032) q[0];
rz(-pi) q[1];
rz(-3.1132142) q[2];
sx q[2];
rz(-2.252223) q[2];
sx q[2];
rz(2.4383557) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7422223) q[1];
sx q[1];
rz(-0.74811799) q[1];
sx q[1];
rz(3.0518603) q[1];
rz(2.378024) q[3];
sx q[3];
rz(-1.1022304) q[3];
sx q[3];
rz(-0.46043049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2406771) q[2];
sx q[2];
rz(-0.42833504) q[2];
sx q[2];
rz(-2.6348616) q[2];
rz(-1.2725376) q[3];
sx q[3];
rz(-1.7839909) q[3];
sx q[3];
rz(2.8633964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5469359) q[0];
sx q[0];
rz(-1.2242067) q[0];
sx q[0];
rz(2.0261672) q[0];
rz(1.8755272) q[1];
sx q[1];
rz(-1.5770117) q[1];
sx q[1];
rz(0.87475264) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6653671) q[0];
sx q[0];
rz(-1.3194808) q[0];
sx q[0];
rz(0.091500207) q[0];
rz(-pi) q[1];
rz(1.5953996) q[2];
sx q[2];
rz(-0.65381351) q[2];
sx q[2];
rz(-3.1179414) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0759039) q[1];
sx q[1];
rz(-1.3049647) q[1];
sx q[1];
rz(-2.7194383) q[1];
x q[2];
rz(2.9625499) q[3];
sx q[3];
rz(-2.2692625) q[3];
sx q[3];
rz(-1.2579789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0227585) q[2];
sx q[2];
rz(-2.6648882) q[2];
sx q[2];
rz(-1.9258707) q[2];
rz(-1.3974961) q[3];
sx q[3];
rz(-1.5236866) q[3];
sx q[3];
rz(-1.4309179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8656411) q[0];
sx q[0];
rz(-3.0396099) q[0];
sx q[0];
rz(1.6708466) q[0];
rz(1.7631081) q[1];
sx q[1];
rz(-2.3536317) q[1];
sx q[1];
rz(-1.4102304) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84274693) q[0];
sx q[0];
rz(-2.0618467) q[0];
sx q[0];
rz(2.9058855) q[0];
rz(1.8115787) q[2];
sx q[2];
rz(-2.2301939) q[2];
sx q[2];
rz(1.3032152) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8961337) q[1];
sx q[1];
rz(-1.2583548) q[1];
sx q[1];
rz(1.5691084) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0038807) q[3];
sx q[3];
rz(-1.2473462) q[3];
sx q[3];
rz(-2.5057305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0146279) q[2];
sx q[2];
rz(-0.85863272) q[2];
sx q[2];
rz(1.8178168) q[2];
rz(1.7823559) q[3];
sx q[3];
rz(-1.2991354) q[3];
sx q[3];
rz(2.7544379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0303665) q[0];
sx q[0];
rz(-1.68196) q[0];
sx q[0];
rz(0.02221814) q[0];
rz(-2.4413595) q[1];
sx q[1];
rz(-1.2513688) q[1];
sx q[1];
rz(-2.1846695) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9807463) q[0];
sx q[0];
rz(-2.2637746) q[0];
sx q[0];
rz(-1.3816383) q[0];
rz(-1.7526723) q[2];
sx q[2];
rz(-1.7851014) q[2];
sx q[2];
rz(2.8150812) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.8808525) q[1];
sx q[1];
rz(-1.2965974) q[1];
sx q[1];
rz(1.9874279) q[1];
x q[2];
rz(0.33958667) q[3];
sx q[3];
rz(-1.9882147) q[3];
sx q[3];
rz(-1.0610895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2172829) q[2];
sx q[2];
rz(-1.1223015) q[2];
sx q[2];
rz(-1.699532) q[2];
rz(-2.0349272) q[3];
sx q[3];
rz(-2.536074) q[3];
sx q[3];
rz(-1.1520011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7149413) q[0];
sx q[0];
rz(-1.8639257) q[0];
sx q[0];
rz(-0.46698025) q[0];
rz(1.8309719) q[1];
sx q[1];
rz(-0.95342031) q[1];
sx q[1];
rz(-0.39042815) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434131) q[0];
sx q[0];
rz(-1.5382086) q[0];
sx q[0];
rz(3.0010953) q[0];
x q[1];
rz(-1.2762345) q[2];
sx q[2];
rz(-1.4678363) q[2];
sx q[2];
rz(-1.1338794) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2749719) q[1];
sx q[1];
rz(-0.83217794) q[1];
sx q[1];
rz(2.2029331) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14777811) q[3];
sx q[3];
rz(-0.96333233) q[3];
sx q[3];
rz(2.1097418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19162576) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(1.0170271) q[2];
rz(0.93938604) q[3];
sx q[3];
rz(-0.87299577) q[3];
sx q[3];
rz(0.92923195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7153213) q[0];
sx q[0];
rz(-2.4763698) q[0];
sx q[0];
rz(1.9287047) q[0];
rz(2.5384278) q[1];
sx q[1];
rz(-2.0096571) q[1];
sx q[1];
rz(2.718198) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42453897) q[0];
sx q[0];
rz(-1.6485212) q[0];
sx q[0];
rz(-0.86941289) q[0];
rz(-pi) q[1];
rz(0.076831623) q[2];
sx q[2];
rz(-1.4203826) q[2];
sx q[2];
rz(-2.7422649) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.054476995) q[1];
sx q[1];
rz(-1.3605788) q[1];
sx q[1];
rz(-1.3940548) q[1];
rz(2.9717507) q[3];
sx q[3];
rz(-1.5780996) q[3];
sx q[3];
rz(-0.23582349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49070552) q[2];
sx q[2];
rz(-2.8287973) q[2];
sx q[2];
rz(0.97839626) q[2];
rz(0.69495106) q[3];
sx q[3];
rz(-1.9415104) q[3];
sx q[3];
rz(-1.8470701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0438743) q[0];
sx q[0];
rz(-0.23625034) q[0];
sx q[0];
rz(3.0978715) q[0];
rz(-1.9678736) q[1];
sx q[1];
rz(-0.81467384) q[1];
sx q[1];
rz(-2.3497605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4872015) q[0];
sx q[0];
rz(-0.7672317) q[0];
sx q[0];
rz(2.5228398) q[0];
rz(-pi) q[1];
rz(1.2709684) q[2];
sx q[2];
rz(-1.764503) q[2];
sx q[2];
rz(1.1739588) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3085295) q[1];
sx q[1];
rz(-1.2385784) q[1];
sx q[1];
rz(-2.019472) q[1];
rz(0.47747647) q[3];
sx q[3];
rz(-1.9100185) q[3];
sx q[3];
rz(-0.71996688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9110079) q[2];
sx q[2];
rz(-2.8318996) q[2];
sx q[2];
rz(1.9045551) q[2];
rz(-2.6593995) q[3];
sx q[3];
rz(-2.2174406) q[3];
sx q[3];
rz(-2.6528416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1173387) q[0];
sx q[0];
rz(-2.2145705) q[0];
sx q[0];
rz(-1.3379958) q[0];
rz(2.2894739) q[1];
sx q[1];
rz(-1.2282649) q[1];
sx q[1];
rz(-1.9911912) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43088461) q[0];
sx q[0];
rz(-1.0477433) q[0];
sx q[0];
rz(2.5780748) q[0];
rz(-pi) q[1];
rz(-1.2537635) q[2];
sx q[2];
rz(-1.46508) q[2];
sx q[2];
rz(-0.52717613) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89429606) q[1];
sx q[1];
rz(-1.1793696) q[1];
sx q[1];
rz(2.9904234) q[1];
rz(-pi) q[2];
rz(0.38693736) q[3];
sx q[3];
rz(-1.0067356) q[3];
sx q[3];
rz(1.2880985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.368025) q[2];
sx q[2];
rz(-2.1358392) q[2];
sx q[2];
rz(-2.8311938) q[2];
rz(-2.5144905) q[3];
sx q[3];
rz(-0.98395625) q[3];
sx q[3];
rz(-1.1712801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2720168) q[0];
sx q[0];
rz(-1.5675114) q[0];
sx q[0];
rz(-1.5935224) q[0];
rz(-1.2177474) q[1];
sx q[1];
rz(-1.4532614) q[1];
sx q[1];
rz(-0.90167602) q[1];
rz(1.1989087) q[2];
sx q[2];
rz(-2.6116284) q[2];
sx q[2];
rz(-2.2378599) q[2];
rz(-0.18219215) q[3];
sx q[3];
rz(-2.8145418) q[3];
sx q[3];
rz(-2.9180632) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
