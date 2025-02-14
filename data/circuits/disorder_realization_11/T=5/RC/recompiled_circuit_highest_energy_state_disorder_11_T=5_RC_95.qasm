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
rz(1.4229245) q[0];
sx q[0];
rz(1.0942425) q[0];
sx q[0];
rz(9.6828238) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(-1.300783) q[1];
sx q[1];
rz(0.25564495) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6171744) q[0];
sx q[0];
rz(-0.99992311) q[0];
sx q[0];
rz(-0.85321315) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7529703) q[2];
sx q[2];
rz(-1.5124694) q[2];
sx q[2];
rz(-1.3921392) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8430184) q[1];
sx q[1];
rz(-2.9461423) q[1];
sx q[1];
rz(1.9324383) q[1];
rz(-1.3575451) q[3];
sx q[3];
rz(-2.0590326) q[3];
sx q[3];
rz(0.31682107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6306182) q[2];
sx q[2];
rz(-1.3773842) q[2];
sx q[2];
rz(1.1154491) q[2];
rz(-3.1274146) q[3];
sx q[3];
rz(-1.7970128) q[3];
sx q[3];
rz(2.7052963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9344591) q[0];
sx q[0];
rz(-0.62196982) q[0];
sx q[0];
rz(-1.8943262) q[0];
rz(0.44218749) q[1];
sx q[1];
rz(-1.3860044) q[1];
sx q[1];
rz(-0.8173379) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2674448) q[0];
sx q[0];
rz(-0.56330452) q[0];
sx q[0];
rz(1.2059709) q[0];
rz(-pi) q[1];
rz(-0.27702443) q[2];
sx q[2];
rz(-1.9961341) q[2];
sx q[2];
rz(2.9556731) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8838447) q[1];
sx q[1];
rz(-2.2040329) q[1];
sx q[1];
rz(-1.3650989) q[1];
x q[2];
rz(-2.2226187) q[3];
sx q[3];
rz(-2.2187244) q[3];
sx q[3];
rz(-2.1369262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3680129) q[2];
sx q[2];
rz(-2.7648338) q[2];
sx q[2];
rz(2.9212941) q[2];
rz(-1.9593272) q[3];
sx q[3];
rz(-1.9672829) q[3];
sx q[3];
rz(0.063145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0910864) q[0];
sx q[0];
rz(-1.8244705) q[0];
sx q[0];
rz(-2.0029946) q[0];
rz(-0.41172045) q[1];
sx q[1];
rz(-0.76985923) q[1];
sx q[1];
rz(1.0134816) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47501144) q[0];
sx q[0];
rz(-1.555976) q[0];
sx q[0];
rz(0.25304742) q[0];
rz(1.8455681) q[2];
sx q[2];
rz(-2.6856075) q[2];
sx q[2];
rz(0.48874654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24559969) q[1];
sx q[1];
rz(-0.87973824) q[1];
sx q[1];
rz(0.57180239) q[1];
rz(-0.8441505) q[3];
sx q[3];
rz(-1.4733847) q[3];
sx q[3];
rz(1.0507492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98075214) q[2];
sx q[2];
rz(-1.1557121) q[2];
sx q[2];
rz(1.8864924) q[2];
rz(-2.8694425) q[3];
sx q[3];
rz(-0.77747074) q[3];
sx q[3];
rz(-3.0009771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18144064) q[0];
sx q[0];
rz(-2.20521) q[0];
sx q[0];
rz(-2.5312359) q[0];
rz(-2.4329674) q[1];
sx q[1];
rz(-1.0162063) q[1];
sx q[1];
rz(-1.5708539) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7371668) q[0];
sx q[0];
rz(-2.9350087) q[0];
sx q[0];
rz(-0.73295672) q[0];
rz(-pi) q[1];
rz(-0.38489401) q[2];
sx q[2];
rz(-2.8334624) q[2];
sx q[2];
rz(2.5727814) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8766986) q[1];
sx q[1];
rz(-1.8421409) q[1];
sx q[1];
rz(-1.62074) q[1];
rz(-pi) q[2];
rz(-2.3409178) q[3];
sx q[3];
rz(-1.281637) q[3];
sx q[3];
rz(-1.9137933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9307956) q[2];
sx q[2];
rz(-1.0127298) q[2];
sx q[2];
rz(2.5861758) q[2];
rz(-0.72426116) q[3];
sx q[3];
rz(-1.1490425) q[3];
sx q[3];
rz(-2.485062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50437462) q[0];
sx q[0];
rz(-1.1906304) q[0];
sx q[0];
rz(-0.36886886) q[0];
rz(-2.8952307) q[1];
sx q[1];
rz(-1.3280222) q[1];
sx q[1];
rz(-1.7074283) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61489366) q[0];
sx q[0];
rz(-1.5656359) q[0];
sx q[0];
rz(-2.3847488) q[0];
x q[1];
rz(2.7857615) q[2];
sx q[2];
rz(-0.47105481) q[2];
sx q[2];
rz(-1.5233153) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1885438) q[1];
sx q[1];
rz(-1.5665788) q[1];
sx q[1];
rz(0.40178816) q[1];
x q[2];
rz(-3.0415972) q[3];
sx q[3];
rz(-1.7971562) q[3];
sx q[3];
rz(-2.5842427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.65115923) q[2];
sx q[2];
rz(-1.8305402) q[2];
sx q[2];
rz(1.0820214) q[2];
rz(2.3467482) q[3];
sx q[3];
rz(-1.4719897) q[3];
sx q[3];
rz(-1.883551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27570462) q[0];
sx q[0];
rz(-0.10144932) q[0];
sx q[0];
rz(-2.8979982) q[0];
rz(0.98681915) q[1];
sx q[1];
rz(-2.3934264) q[1];
sx q[1];
rz(2.2056244) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075386062) q[0];
sx q[0];
rz(-2.6200326) q[0];
sx q[0];
rz(-0.88514502) q[0];
x q[1];
rz(-2.7924839) q[2];
sx q[2];
rz(-2.287068) q[2];
sx q[2];
rz(2.8385065) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3159891) q[1];
sx q[1];
rz(-2.2554419) q[1];
sx q[1];
rz(1.0427834) q[1];
rz(-pi) q[2];
rz(3.1201911) q[3];
sx q[3];
rz(-2.170522) q[3];
sx q[3];
rz(0.76364005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.662107) q[2];
sx q[2];
rz(-1.138849) q[2];
sx q[2];
rz(0.34995079) q[2];
rz(0.55772603) q[3];
sx q[3];
rz(-1.2099268) q[3];
sx q[3];
rz(2.2902655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6996985) q[0];
sx q[0];
rz(-2.5040369) q[0];
sx q[0];
rz(0.30174524) q[0];
rz(-2.8217577) q[1];
sx q[1];
rz(-1.539307) q[1];
sx q[1];
rz(1.1357657) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1162735) q[0];
sx q[0];
rz(-2.2145382) q[0];
sx q[0];
rz(-1.6380861) q[0];
rz(-pi) q[1];
rz(-2.2280424) q[2];
sx q[2];
rz(-0.74591178) q[2];
sx q[2];
rz(2.8532956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9965976) q[1];
sx q[1];
rz(-0.97255822) q[1];
sx q[1];
rz(-0.096417565) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3507422) q[3];
sx q[3];
rz(-1.491427) q[3];
sx q[3];
rz(-2.1805003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7319506) q[2];
sx q[2];
rz(-0.66764098) q[2];
sx q[2];
rz(-2.9988334) q[2];
rz(-1.7536633) q[3];
sx q[3];
rz(-1.9526491) q[3];
sx q[3];
rz(0.68283844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9614354) q[0];
sx q[0];
rz(-0.016594369) q[0];
sx q[0];
rz(2.5855682) q[0];
rz(2.9122638) q[1];
sx q[1];
rz(-1.2622958) q[1];
sx q[1];
rz(1.75846) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1616104) q[0];
sx q[0];
rz(-1.605827) q[0];
sx q[0];
rz(2.4039414) q[0];
rz(-pi) q[1];
rz(2.6830693) q[2];
sx q[2];
rz(-1.5178198) q[2];
sx q[2];
rz(-2.3725739) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3179211) q[1];
sx q[1];
rz(-1.8601728) q[1];
sx q[1];
rz(-1.1422864) q[1];
rz(-pi) q[2];
rz(-1.1874299) q[3];
sx q[3];
rz(-1.1463506) q[3];
sx q[3];
rz(-2.2147873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6250299) q[2];
sx q[2];
rz(-2.8690858) q[2];
sx q[2];
rz(-1.9005091) q[2];
rz(-0.062601335) q[3];
sx q[3];
rz(-1.4227941) q[3];
sx q[3];
rz(-0.82714287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0144219) q[0];
sx q[0];
rz(-1.7272471) q[0];
sx q[0];
rz(-0.14778368) q[0];
rz(-0.5087018) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(-0.37857372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62703778) q[0];
sx q[0];
rz(-2.6588024) q[0];
sx q[0];
rz(-0.52282368) q[0];
rz(-2.1630387) q[2];
sx q[2];
rz(-1.8752974) q[2];
sx q[2];
rz(-2.1155807) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8151103) q[1];
sx q[1];
rz(-1.2604144) q[1];
sx q[1];
rz(1.923345) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76102961) q[3];
sx q[3];
rz(-1.4806403) q[3];
sx q[3];
rz(-0.48616274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96616894) q[2];
sx q[2];
rz(-0.95169008) q[2];
sx q[2];
rz(-1.2790722) q[2];
rz(-1.7133948) q[3];
sx q[3];
rz(-2.1396075) q[3];
sx q[3];
rz(0.14555791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3928423) q[0];
sx q[0];
rz(-0.57806438) q[0];
sx q[0];
rz(-3.0774935) q[0];
rz(-0.52866689) q[1];
sx q[1];
rz(-1.268498) q[1];
sx q[1];
rz(-0.97506964) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.753016) q[0];
sx q[0];
rz(-0.56682366) q[0];
sx q[0];
rz(-2.6334769) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98627536) q[2];
sx q[2];
rz(-2.0615675) q[2];
sx q[2];
rz(-1.402439) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9916315) q[1];
sx q[1];
rz(-1.35222) q[1];
sx q[1];
rz(-2.6944955) q[1];
x q[2];
rz(-0.24419489) q[3];
sx q[3];
rz(-1.3483485) q[3];
sx q[3];
rz(-0.73058587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9958682) q[2];
sx q[2];
rz(-0.66103649) q[2];
sx q[2];
rz(-0.60689849) q[2];
rz(0.25035826) q[3];
sx q[3];
rz(-1.4162049) q[3];
sx q[3];
rz(0.50601602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.09457) q[0];
sx q[0];
rz(-1.2170412) q[0];
sx q[0];
rz(1.8893597) q[0];
rz(0.49867123) q[1];
sx q[1];
rz(-1.55232) q[1];
sx q[1];
rz(1.3943863) q[1];
rz(2.278419) q[2];
sx q[2];
rz(-2.045608) q[2];
sx q[2];
rz(1.506293) q[2];
rz(1.0389497) q[3];
sx q[3];
rz(-0.23775756) q[3];
sx q[3];
rz(0.47142117) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
