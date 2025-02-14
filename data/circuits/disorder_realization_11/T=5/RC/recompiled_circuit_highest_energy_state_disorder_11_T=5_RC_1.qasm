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
rz(-1.7186681) q[0];
sx q[0];
rz(-1.0942425) q[0];
sx q[0];
rz(-2.8835468) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(-1.300783) q[1];
sx q[1];
rz(0.25564495) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6545171) q[0];
sx q[0];
rz(-2.1571113) q[0];
sx q[0];
rz(0.70588995) q[0];
rz(-pi) q[1];
rz(-2.7529703) q[2];
sx q[2];
rz(-1.6291233) q[2];
sx q[2];
rz(1.7494534) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.62754831) q[1];
sx q[1];
rz(-1.5020292) q[1];
sx q[1];
rz(1.7538944) q[1];
rz(-2.643804) q[3];
sx q[3];
rz(-1.7588132) q[3];
sx q[3];
rz(1.1527485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51097441) q[2];
sx q[2];
rz(-1.7642085) q[2];
sx q[2];
rz(1.1154491) q[2];
rz(-0.014178064) q[3];
sx q[3];
rz(-1.3445798) q[3];
sx q[3];
rz(-0.43629638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(1.9344591) q[0];
sx q[0];
rz(-2.5196228) q[0];
sx q[0];
rz(-1.2472664) q[0];
rz(0.44218749) q[1];
sx q[1];
rz(-1.7555883) q[1];
sx q[1];
rz(-2.3242548) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2674448) q[0];
sx q[0];
rz(-2.5782881) q[0];
sx q[0];
rz(-1.2059709) q[0];
rz(2.8645682) q[2];
sx q[2];
rz(-1.9961341) q[2];
sx q[2];
rz(-0.1859196) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0603588) q[1];
sx q[1];
rz(-0.66141719) q[1];
sx q[1];
rz(2.8702186) q[1];
rz(0.91897398) q[3];
sx q[3];
rz(-0.92286829) q[3];
sx q[3];
rz(2.1369262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3680129) q[2];
sx q[2];
rz(-2.7648338) q[2];
sx q[2];
rz(-0.22029857) q[2];
rz(1.9593272) q[3];
sx q[3];
rz(-1.1743098) q[3];
sx q[3];
rz(-3.0784472) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0910864) q[0];
sx q[0];
rz(-1.8244705) q[0];
sx q[0];
rz(1.1385981) q[0];
rz(2.7298722) q[1];
sx q[1];
rz(-2.3717334) q[1];
sx q[1];
rz(2.128111) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47501144) q[0];
sx q[0];
rz(-1.555976) q[0];
sx q[0];
rz(2.8885452) q[0];
rz(-pi) q[1];
rz(3.0092952) q[2];
sx q[2];
rz(-1.1331285) q[2];
sx q[2];
rz(-0.79298151) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.24559969) q[1];
sx q[1];
rz(-2.2618544) q[1];
sx q[1];
rz(2.5697903) q[1];
rz(0.8441505) q[3];
sx q[3];
rz(-1.6682079) q[3];
sx q[3];
rz(1.0507492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98075214) q[2];
sx q[2];
rz(-1.9858805) q[2];
sx q[2];
rz(-1.2551003) q[2];
rz(0.27215019) q[3];
sx q[3];
rz(-0.77747074) q[3];
sx q[3];
rz(-3.0009771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.960152) q[0];
sx q[0];
rz(-0.93638268) q[0];
sx q[0];
rz(-0.61035672) q[0];
rz(2.4329674) q[1];
sx q[1];
rz(-2.1253864) q[1];
sx q[1];
rz(-1.5708539) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40442586) q[0];
sx q[0];
rz(-0.20658399) q[0];
sx q[0];
rz(-0.73295672) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4518634) q[2];
sx q[2];
rz(-1.2858675) q[2];
sx q[2];
rz(-0.97078427) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.849087) q[1];
sx q[1];
rz(-1.6189112) q[1];
sx q[1];
rz(2.8699257) q[1];
rz(2.7486984) q[3];
sx q[3];
rz(-0.84019444) q[3];
sx q[3];
rz(0.073176633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9307956) q[2];
sx q[2];
rz(-2.1288629) q[2];
sx q[2];
rz(-0.55541682) q[2];
rz(0.72426116) q[3];
sx q[3];
rz(-1.1490425) q[3];
sx q[3];
rz(2.485062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50437462) q[0];
sx q[0];
rz(-1.9509622) q[0];
sx q[0];
rz(0.36886886) q[0];
rz(-2.8952307) q[1];
sx q[1];
rz(-1.8135704) q[1];
sx q[1];
rz(1.7074283) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61489366) q[0];
sx q[0];
rz(-1.5656359) q[0];
sx q[0];
rz(-2.3847488) q[0];
rz(-1.3952012) q[2];
sx q[2];
rz(-2.0101974) q[2];
sx q[2];
rz(1.1281769) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.37232698) q[1];
sx q[1];
rz(-0.40180909) q[1];
sx q[1];
rz(-0.010784464) q[1];
x q[2];
rz(-3.0415972) q[3];
sx q[3];
rz(-1.3444364) q[3];
sx q[3];
rz(2.5842427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.65115923) q[2];
sx q[2];
rz(-1.3110524) q[2];
sx q[2];
rz(-1.0820214) q[2];
rz(-2.3467482) q[3];
sx q[3];
rz(-1.4719897) q[3];
sx q[3];
rz(-1.2580416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.865888) q[0];
sx q[0];
rz(-0.10144932) q[0];
sx q[0];
rz(-0.24359447) q[0];
rz(-2.1547735) q[1];
sx q[1];
rz(-0.74816626) q[1];
sx q[1];
rz(-2.2056244) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2631131) q[0];
sx q[0];
rz(-1.8917483) q[0];
sx q[0];
rz(1.9892938) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34910874) q[2];
sx q[2];
rz(-0.85452467) q[2];
sx q[2];
rz(-2.8385065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82560357) q[1];
sx q[1];
rz(-2.2554419) q[1];
sx q[1];
rz(-2.0988093) q[1];
rz(-pi) q[2];
rz(3.1201911) q[3];
sx q[3];
rz(-0.97107065) q[3];
sx q[3];
rz(2.3779526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4794856) q[2];
sx q[2];
rz(-1.138849) q[2];
sx q[2];
rz(0.34995079) q[2];
rz(-0.55772603) q[3];
sx q[3];
rz(-1.9316659) q[3];
sx q[3];
rz(2.2902655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6996985) q[0];
sx q[0];
rz(-2.5040369) q[0];
sx q[0];
rz(0.30174524) q[0];
rz(0.31983495) q[1];
sx q[1];
rz(-1.539307) q[1];
sx q[1];
rz(1.1357657) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58590305) q[0];
sx q[0];
rz(-1.5169889) q[0];
sx q[0];
rz(2.4967628) q[0];
rz(0.93923969) q[2];
sx q[2];
rz(-1.9983091) q[2];
sx q[2];
rz(1.3434501) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.144995) q[1];
sx q[1];
rz(-2.1690344) q[1];
sx q[1];
rz(3.0451751) q[1];
x q[2];
rz(-3.0301827) q[3];
sx q[3];
rz(-2.3476331) q[3];
sx q[3];
rz(-0.68796989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7319506) q[2];
sx q[2];
rz(-2.4739517) q[2];
sx q[2];
rz(0.14275924) q[2];
rz(-1.7536633) q[3];
sx q[3];
rz(-1.9526491) q[3];
sx q[3];
rz(0.68283844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.1801572) q[0];
sx q[0];
rz(-0.016594369) q[0];
sx q[0];
rz(0.55602443) q[0];
rz(2.9122638) q[1];
sx q[1];
rz(-1.8792968) q[1];
sx q[1];
rz(1.3831327) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62263238) q[0];
sx q[0];
rz(-0.83370249) q[0];
sx q[0];
rz(1.5234768) q[0];
x q[1];
rz(-3.0223614) q[2];
sx q[2];
rz(-2.6802353) q[2];
sx q[2];
rz(-0.90864321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.38249967) q[1];
sx q[1];
rz(-1.9803932) q[1];
sx q[1];
rz(-0.31633693) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4501958) q[3];
sx q[3];
rz(-2.5775839) q[3];
sx q[3];
rz(-0.15182748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6250299) q[2];
sx q[2];
rz(-2.8690858) q[2];
sx q[2];
rz(1.2410835) q[2];
rz(0.062601335) q[3];
sx q[3];
rz(-1.7187985) q[3];
sx q[3];
rz(-0.82714287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1271707) q[0];
sx q[0];
rz(-1.7272471) q[0];
sx q[0];
rz(0.14778368) q[0];
rz(-0.5087018) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(-0.37857372) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4157279) q[0];
sx q[0];
rz(-1.8047389) q[0];
sx q[0];
rz(-2.7153003) q[0];
x q[1];
rz(-0.97855391) q[2];
sx q[2];
rz(-1.2662953) q[2];
sx q[2];
rz(1.026012) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.55146688) q[1];
sx q[1];
rz(-0.46537897) q[1];
sx q[1];
rz(-2.3193633) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6949953) q[3];
sx q[3];
rz(-0.81362766) q[3];
sx q[3];
rz(-0.99909335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1754237) q[2];
sx q[2];
rz(-2.1899026) q[2];
sx q[2];
rz(1.8625205) q[2];
rz(1.7133948) q[3];
sx q[3];
rz(-2.1396075) q[3];
sx q[3];
rz(2.9960347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3928423) q[0];
sx q[0];
rz(-2.5635283) q[0];
sx q[0];
rz(3.0774935) q[0];
rz(0.52866689) q[1];
sx q[1];
rz(-1.268498) q[1];
sx q[1];
rz(-2.166523) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9720358) q[0];
sx q[0];
rz(-1.0825048) q[0];
sx q[0];
rz(1.2704865) q[0];
rz(-0.56985241) q[2];
sx q[2];
rz(-2.0790711) q[2];
sx q[2];
rz(-3.0076671) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8243557) q[1];
sx q[1];
rz(-2.0065161) q[1];
sx q[1];
rz(1.3292666) q[1];
x q[2];
rz(0.75210877) q[3];
sx q[3];
rz(-2.8127713) q[3];
sx q[3];
rz(-3.0260835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.14572445) q[2];
sx q[2];
rz(-0.66103649) q[2];
sx q[2];
rz(-2.5346942) q[2];
rz(0.25035826) q[3];
sx q[3];
rz(-1.7253877) q[3];
sx q[3];
rz(2.6355766) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.09457) q[0];
sx q[0];
rz(-1.9245514) q[0];
sx q[0];
rz(-1.2522329) q[0];
rz(2.6429214) q[1];
sx q[1];
rz(-1.5892727) q[1];
sx q[1];
rz(-1.7472063) q[1];
rz(-2.5468536) q[2];
sx q[2];
rz(-0.95437106) q[2];
sx q[2];
rz(0.30827733) q[2];
rz(-1.7767033) q[3];
sx q[3];
rz(-1.4510703) q[3];
sx q[3];
rz(2.5616796) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
