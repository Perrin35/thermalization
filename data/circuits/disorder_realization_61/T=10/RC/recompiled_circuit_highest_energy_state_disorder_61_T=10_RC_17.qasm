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
rz(-0.65499175) q[0];
rz(-1.4870149) q[1];
sx q[1];
rz(-1.4227285) q[1];
sx q[1];
rz(-1.8185599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0812936) q[0];
sx q[0];
rz(-2.4954688) q[0];
sx q[0];
rz(2.6160014) q[0];
rz(-0.31794117) q[2];
sx q[2];
rz(-1.7934718) q[2];
sx q[2];
rz(0.26938619) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2046584) q[1];
sx q[1];
rz(-1.2950674) q[1];
sx q[1];
rz(0.10358056) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72812702) q[3];
sx q[3];
rz(-1.9890729) q[3];
sx q[3];
rz(0.74467105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46372867) q[2];
sx q[2];
rz(-0.68631309) q[2];
sx q[2];
rz(0.65832552) q[2];
rz(2.5074734) q[3];
sx q[3];
rz(-1.442797) q[3];
sx q[3];
rz(1.770795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089791678) q[0];
sx q[0];
rz(-2.8103516) q[0];
sx q[0];
rz(-0.077089699) q[0];
rz(-0.58473051) q[1];
sx q[1];
rz(-1.3405776) q[1];
sx q[1];
rz(1.1999406) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18835078) q[0];
sx q[0];
rz(-2.6380499) q[0];
sx q[0];
rz(1.153323) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7464094) q[2];
sx q[2];
rz(-2.9746911) q[2];
sx q[2];
rz(1.6957972) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2421094) q[1];
sx q[1];
rz(-2.1596908) q[1];
sx q[1];
rz(-0.45581006) q[1];
rz(0.097278519) q[3];
sx q[3];
rz(-2.5803498) q[3];
sx q[3];
rz(1.2429855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8889019) q[2];
sx q[2];
rz(-0.98068792) q[2];
sx q[2];
rz(0.67548951) q[2];
rz(0.0096970079) q[3];
sx q[3];
rz(-2.8986425) q[3];
sx q[3];
rz(3.1233136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6460687) q[0];
sx q[0];
rz(-1.6697474) q[0];
sx q[0];
rz(-0.75463265) q[0];
rz(-0.010628788) q[1];
sx q[1];
rz(-1.7517215) q[1];
sx q[1];
rz(2.1307814) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0389688) q[0];
sx q[0];
rz(-1.083923) q[0];
sx q[0];
rz(0.34619934) q[0];
rz(-pi) q[1];
x q[1];
rz(0.028378475) q[2];
sx q[2];
rz(-0.88936964) q[2];
sx q[2];
rz(-2.4383557) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.52151187) q[1];
sx q[1];
rz(-0.82640582) q[1];
sx q[1];
rz(-1.6537731) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63186462) q[3];
sx q[3];
rz(-0.87040983) q[3];
sx q[3];
rz(-2.4720342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2406771) q[2];
sx q[2];
rz(-0.42833504) q[2];
sx q[2];
rz(0.50673103) q[2];
rz(-1.2725376) q[3];
sx q[3];
rz(-1.7839909) q[3];
sx q[3];
rz(-0.27819628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5469359) q[0];
sx q[0];
rz(-1.2242067) q[0];
sx q[0];
rz(1.1154255) q[0];
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
rz(-2.311888) q[0];
sx q[0];
rz(-0.26712298) q[0];
sx q[0];
rz(-1.2288837) q[0];
x q[1];
rz(-1.5953996) q[2];
sx q[2];
rz(-0.65381351) q[2];
sx q[2];
rz(-0.023651274) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.03434788) q[1];
sx q[1];
rz(-0.49458359) q[1];
sx q[1];
rz(0.58652189) q[1];
rz(2.277209) q[3];
sx q[3];
rz(-1.7076075) q[3];
sx q[3];
rz(-2.9446326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11883417) q[2];
sx q[2];
rz(-2.6648882) q[2];
sx q[2];
rz(1.215722) q[2];
rz(1.7440965) q[3];
sx q[3];
rz(-1.5236866) q[3];
sx q[3];
rz(1.7106748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8656411) q[0];
sx q[0];
rz(-0.10198274) q[0];
sx q[0];
rz(-1.6708466) q[0];
rz(1.3784846) q[1];
sx q[1];
rz(-2.3536317) q[1];
sx q[1];
rz(1.4102304) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8277692) q[0];
sx q[0];
rz(-2.6010989) q[0];
sx q[0];
rz(1.9825516) q[0];
rz(-1.8115787) q[2];
sx q[2];
rz(-2.2301939) q[2];
sx q[2];
rz(-1.3032152) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2454589) q[1];
sx q[1];
rz(-1.2583548) q[1];
sx q[1];
rz(-1.5724843) q[1];
rz(-pi) q[2];
rz(2.7878109) q[3];
sx q[3];
rz(-1.980034) q[3];
sx q[3];
rz(-2.0607467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0146279) q[2];
sx q[2];
rz(-2.2829599) q[2];
sx q[2];
rz(1.8178168) q[2];
rz(-1.7823559) q[3];
sx q[3];
rz(-1.8424572) q[3];
sx q[3];
rz(-0.38715473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(1.0303665) q[0];
sx q[0];
rz(-1.68196) q[0];
sx q[0];
rz(3.1193745) q[0];
rz(-0.70023316) q[1];
sx q[1];
rz(-1.8902238) q[1];
sx q[1];
rz(0.95692316) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2882521) q[0];
sx q[0];
rz(-1.4256251) q[0];
sx q[0];
rz(2.4397544) q[0];
rz(-2.9238052) q[2];
sx q[2];
rz(-1.7484669) q[2];
sx q[2];
rz(-1.9363994) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0007504) q[1];
sx q[1];
rz(-0.4943119) q[1];
sx q[1];
rz(-0.96338455) q[1];
rz(-pi) q[2];
rz(-2.2150008) q[3];
sx q[3];
rz(-0.53172382) q[3];
sx q[3];
rz(0.34429541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9243098) q[2];
sx q[2];
rz(-2.0192912) q[2];
sx q[2];
rz(-1.4420606) q[2];
rz(-1.1066655) q[3];
sx q[3];
rz(-0.60551867) q[3];
sx q[3];
rz(-1.1520011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098179558) q[0];
sx q[0];
rz(-1.603384) q[0];
sx q[0];
rz(0.14049732) q[0];
rz(-pi) q[1];
rz(3.0340334) q[2];
sx q[2];
rz(-1.8637519) q[2];
sx q[2];
rz(2.7358472) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7539333) q[1];
sx q[1];
rz(-2.0229335) q[1];
sx q[1];
rz(-2.2958295) q[1];
rz(-pi) q[2];
rz(-1.3620699) q[3];
sx q[3];
rz(-0.62297076) q[3];
sx q[3];
rz(-1.8546212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19162576) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(2.1245655) q[2];
rz(2.2022066) q[3];
sx q[3];
rz(-2.2685969) q[3];
sx q[3];
rz(-2.2123607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4262714) q[0];
sx q[0];
rz(-0.66522288) q[0];
sx q[0];
rz(1.2128879) q[0];
rz(2.5384278) q[1];
sx q[1];
rz(-2.0096571) q[1];
sx q[1];
rz(2.718198) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0807665) q[0];
sx q[0];
rz(-0.87196022) q[0];
sx q[0];
rz(3.0399975) q[0];
rz(-pi) q[1];
x q[1];
rz(0.076831623) q[2];
sx q[2];
rz(-1.4203826) q[2];
sx q[2];
rz(0.39932775) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.054476995) q[1];
sx q[1];
rz(-1.3605788) q[1];
sx q[1];
rz(1.3940548) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5782062) q[3];
sx q[3];
rz(-1.7406337) q[3];
sx q[3];
rz(1.8078723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49070552) q[2];
sx q[2];
rz(-0.31279534) q[2];
sx q[2];
rz(-0.97839626) q[2];
rz(-0.69495106) q[3];
sx q[3];
rz(-1.2000822) q[3];
sx q[3];
rz(1.2945226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0438743) q[0];
sx q[0];
rz(-2.9053423) q[0];
sx q[0];
rz(-0.043721113) q[0];
rz(1.9678736) q[1];
sx q[1];
rz(-2.3269188) q[1];
sx q[1];
rz(0.79183212) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2671473) q[0];
sx q[0];
rz(-2.1717779) q[0];
sx q[0];
rz(-2.0807666) q[0];
rz(-0.20251198) q[2];
sx q[2];
rz(-1.8648476) q[2];
sx q[2];
rz(-2.685315) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83306317) q[1];
sx q[1];
rz(-1.9030142) q[1];
sx q[1];
rz(-2.019472) q[1];
rz(-pi) q[2];
x q[2];
rz(1.94897) q[3];
sx q[3];
rz(-2.0190051) q[3];
sx q[3];
rz(-1.0213272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9110079) q[2];
sx q[2];
rz(-0.30969301) q[2];
sx q[2];
rz(1.9045551) q[2];
rz(2.6593995) q[3];
sx q[3];
rz(-2.2174406) q[3];
sx q[3];
rz(2.6528416) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.024254) q[0];
sx q[0];
rz(-2.2145705) q[0];
sx q[0];
rz(-1.3379958) q[0];
rz(-0.85211873) q[1];
sx q[1];
rz(-1.2282649) q[1];
sx q[1];
rz(1.1504014) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.710708) q[0];
sx q[0];
rz(-2.0938494) q[0];
sx q[0];
rz(-2.5780748) q[0];
x q[1];
rz(1.2537635) q[2];
sx q[2];
rz(-1.6765127) q[2];
sx q[2];
rz(-0.52717613) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.51440367) q[1];
sx q[1];
rz(-0.41819388) q[1];
sx q[1];
rz(1.2209284) q[1];
x q[2];
rz(-0.97148599) q[3];
sx q[3];
rz(-1.8953634) q[3];
sx q[3];
rz(2.6443986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77356768) q[2];
sx q[2];
rz(-2.1358392) q[2];
sx q[2];
rz(2.8311938) q[2];
rz(0.6271022) q[3];
sx q[3];
rz(-2.1576364) q[3];
sx q[3];
rz(-1.9703126) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2720168) q[0];
sx q[0];
rz(-1.5675114) q[0];
sx q[0];
rz(-1.5935224) q[0];
rz(1.9238453) q[1];
sx q[1];
rz(-1.4532614) q[1];
sx q[1];
rz(-0.90167602) q[1];
rz(1.1989087) q[2];
sx q[2];
rz(-2.6116284) q[2];
sx q[2];
rz(-2.2378599) q[2];
rz(2.8195856) q[3];
sx q[3];
rz(-1.6290355) q[3];
sx q[3];
rz(1.6215948) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
