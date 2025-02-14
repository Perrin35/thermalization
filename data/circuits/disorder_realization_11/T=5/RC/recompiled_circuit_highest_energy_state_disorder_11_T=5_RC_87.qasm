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
rz(-2.0473502) q[0];
sx q[0];
rz(2.8835468) q[0];
rz(-0.99335042) q[1];
sx q[1];
rz(-1.8408096) q[1];
sx q[1];
rz(2.8859477) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50798049) q[0];
sx q[0];
rz(-2.2574212) q[0];
sx q[0];
rz(0.79721862) q[0];
rz(-pi) q[1];
rz(-1.5077816) q[2];
sx q[2];
rz(-1.9587226) q[2];
sx q[2];
rz(-0.20252075) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2985743) q[1];
sx q[1];
rz(-2.9461423) q[1];
sx q[1];
rz(-1.9324383) q[1];
rz(-pi) q[2];
rz(-2.643804) q[3];
sx q[3];
rz(-1.3827795) q[3];
sx q[3];
rz(1.9888442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51097441) q[2];
sx q[2];
rz(-1.7642085) q[2];
sx q[2];
rz(1.1154491) q[2];
rz(-3.1274146) q[3];
sx q[3];
rz(-1.7970128) q[3];
sx q[3];
rz(-0.43629638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9344591) q[0];
sx q[0];
rz(-0.62196982) q[0];
sx q[0];
rz(-1.8943262) q[0];
rz(-0.44218749) q[1];
sx q[1];
rz(-1.7555883) q[1];
sx q[1];
rz(2.3242548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2674448) q[0];
sx q[0];
rz(-0.56330452) q[0];
sx q[0];
rz(1.2059709) q[0];
rz(-pi) q[1];
rz(-1.1306612) q[2];
sx q[2];
rz(-1.3190184) q[2];
sx q[2];
rz(-1.2680858) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8838447) q[1];
sx q[1];
rz(-0.93755975) q[1];
sx q[1];
rz(-1.3650989) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2226187) q[3];
sx q[3];
rz(-2.2187244) q[3];
sx q[3];
rz(1.0046665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77357972) q[2];
sx q[2];
rz(-0.37675884) q[2];
sx q[2];
rz(-0.22029857) q[2];
rz(-1.9593272) q[3];
sx q[3];
rz(-1.9672829) q[3];
sx q[3];
rz(0.063145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0910864) q[0];
sx q[0];
rz(-1.3171221) q[0];
sx q[0];
rz(2.0029946) q[0];
rz(0.41172045) q[1];
sx q[1];
rz(-2.3717334) q[1];
sx q[1];
rz(-2.128111) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0419755) q[0];
sx q[0];
rz(-1.3177773) q[0];
sx q[0];
rz(1.5861041) q[0];
rz(-pi) q[1];
rz(-1.2960245) q[2];
sx q[2];
rz(-2.6856075) q[2];
sx q[2];
rz(-2.6528461) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0358082) q[1];
sx q[1];
rz(-0.86584751) q[1];
sx q[1];
rz(-2.1501599) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2974422) q[3];
sx q[3];
rz(-1.6682079) q[3];
sx q[3];
rz(1.0507492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98075214) q[2];
sx q[2];
rz(-1.1557121) q[2];
sx q[2];
rz(-1.8864924) q[2];
rz(-2.8694425) q[3];
sx q[3];
rz(-0.77747074) q[3];
sx q[3];
rz(0.14061558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18144064) q[0];
sx q[0];
rz(-0.93638268) q[0];
sx q[0];
rz(2.5312359) q[0];
rz(-0.70862526) q[1];
sx q[1];
rz(-1.0162063) q[1];
sx q[1];
rz(-1.5707387) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7371668) q[0];
sx q[0];
rz(-2.9350087) q[0];
sx q[0];
rz(2.4086359) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38489401) q[2];
sx q[2];
rz(-2.8334624) q[2];
sx q[2];
rz(0.56881126) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29250568) q[1];
sx q[1];
rz(-1.5226814) q[1];
sx q[1];
rz(0.27166697) q[1];
x q[2];
rz(0.39289423) q[3];
sx q[3];
rz(-2.3013982) q[3];
sx q[3];
rz(-3.068416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2107971) q[2];
sx q[2];
rz(-2.1288629) q[2];
sx q[2];
rz(-0.55541682) q[2];
rz(-2.4173315) q[3];
sx q[3];
rz(-1.1490425) q[3];
sx q[3];
rz(-0.65653062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.637218) q[0];
sx q[0];
rz(-1.9509622) q[0];
sx q[0];
rz(-2.7727238) q[0];
rz(-2.8952307) q[1];
sx q[1];
rz(-1.3280222) q[1];
sx q[1];
rz(1.4341644) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95043889) q[0];
sx q[0];
rz(-2.3847347) q[0];
sx q[0];
rz(-3.1340772) q[0];
x q[1];
rz(-0.35583115) q[2];
sx q[2];
rz(-0.47105481) q[2];
sx q[2];
rz(-1.5233153) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.37232698) q[1];
sx q[1];
rz(-2.7397836) q[1];
sx q[1];
rz(3.1308082) q[1];
x q[2];
rz(-1.1617817) q[3];
sx q[3];
rz(-2.8944765) q[3];
sx q[3];
rz(2.1638526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4904334) q[2];
sx q[2];
rz(-1.8305402) q[2];
sx q[2];
rz(2.0595713) q[2];
rz(0.79484445) q[3];
sx q[3];
rz(-1.4719897) q[3];
sx q[3];
rz(1.883551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27570462) q[0];
sx q[0];
rz(-0.10144932) q[0];
sx q[0];
rz(2.8979982) q[0];
rz(-2.1547735) q[1];
sx q[1];
rz(-2.3934264) q[1];
sx q[1];
rz(2.2056244) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0662066) q[0];
sx q[0];
rz(-0.52156007) q[0];
sx q[0];
rz(-0.88514502) q[0];
rz(2.317993) q[2];
sx q[2];
rz(-1.3098426) q[2];
sx q[2];
rz(1.5023155) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0609801) q[1];
sx q[1];
rz(-0.83773936) q[1];
sx q[1];
rz(0.55292801) q[1];
x q[2];
rz(-0.021401568) q[3];
sx q[3];
rz(-2.170522) q[3];
sx q[3];
rz(0.76364005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.662107) q[2];
sx q[2];
rz(-2.0027436) q[2];
sx q[2];
rz(-2.7916419) q[2];
rz(-2.5838666) q[3];
sx q[3];
rz(-1.9316659) q[3];
sx q[3];
rz(-2.2902655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4418942) q[0];
sx q[0];
rz(-2.5040369) q[0];
sx q[0];
rz(0.30174524) q[0];
rz(-0.31983495) q[1];
sx q[1];
rz(-1.6022857) q[1];
sx q[1];
rz(1.1357657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9135043) q[0];
sx q[0];
rz(-2.4948409) q[0];
sx q[0];
rz(0.089368377) q[0];
x q[1];
rz(-2.202353) q[2];
sx q[2];
rz(-1.1432836) q[2];
sx q[2];
rz(1.7981426) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7702076) q[1];
sx q[1];
rz(-1.4911629) q[1];
sx q[1];
rz(-2.1712028) q[1];
x q[2];
rz(0.79085042) q[3];
sx q[3];
rz(-1.6501657) q[3];
sx q[3];
rz(-2.1805003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7319506) q[2];
sx q[2];
rz(-0.66764098) q[2];
sx q[2];
rz(0.14275924) q[2];
rz(1.3879294) q[3];
sx q[3];
rz(-1.9526491) q[3];
sx q[3];
rz(0.68283844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1801572) q[0];
sx q[0];
rz(-0.016594369) q[0];
sx q[0];
rz(-2.5855682) q[0];
rz(0.22932886) q[1];
sx q[1];
rz(-1.8792968) q[1];
sx q[1];
rz(1.75846) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.589298) q[0];
sx q[0];
rz(-0.73832608) q[0];
sx q[0];
rz(-0.05206042) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45852335) q[2];
sx q[2];
rz(-1.6237729) q[2];
sx q[2];
rz(-0.76901877) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.759093) q[1];
sx q[1];
rz(-1.1611995) q[1];
sx q[1];
rz(0.31633693) q[1];
rz(0.45342584) q[3];
sx q[3];
rz(-1.2229706) q[3];
sx q[3];
rz(-2.3330101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51656276) q[2];
sx q[2];
rz(-0.27250686) q[2];
sx q[2];
rz(1.9005091) q[2];
rz(-3.0789913) q[3];
sx q[3];
rz(-1.7187985) q[3];
sx q[3];
rz(2.3144498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
rz(2.6328909) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(-0.37857372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5145549) q[0];
sx q[0];
rz(-0.48279027) q[0];
sx q[0];
rz(-0.52282368) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0835593) q[2];
sx q[2];
rz(-2.4840925) q[2];
sx q[2];
rz(-3.0160144) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3562153) q[1];
sx q[1];
rz(-1.9058203) q[1];
sx q[1];
rz(0.32932333) q[1];
rz(-pi) q[2];
rz(3.0112565) q[3];
sx q[3];
rz(-0.76528434) q[3];
sx q[3];
rz(1.9627067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.96616894) q[2];
sx q[2];
rz(-2.1899026) q[2];
sx q[2];
rz(-1.2790722) q[2];
rz(-1.4281979) q[3];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7487504) q[0];
sx q[0];
rz(-0.57806438) q[0];
sx q[0];
rz(0.064099126) q[0];
rz(0.52866689) q[1];
sx q[1];
rz(-1.268498) q[1];
sx q[1];
rz(0.97506964) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1695568) q[0];
sx q[0];
rz(-2.0590879) q[0];
sx q[0];
rz(-1.8711062) q[0];
x q[1];
rz(-0.80143546) q[2];
sx q[2];
rz(-2.3972627) q[2];
sx q[2];
rz(-0.78730415) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1499612) q[1];
sx q[1];
rz(-1.35222) q[1];
sx q[1];
rz(-2.6944955) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75210877) q[3];
sx q[3];
rz(-0.32882133) q[3];
sx q[3];
rz(-0.11550918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.14572445) q[2];
sx q[2];
rz(-0.66103649) q[2];
sx q[2];
rz(2.5346942) q[2];
rz(-0.25035826) q[3];
sx q[3];
rz(-1.4162049) q[3];
sx q[3];
rz(-0.50601602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0470227) q[0];
sx q[0];
rz(-1.9245514) q[0];
sx q[0];
rz(-1.2522329) q[0];
rz(-2.6429214) q[1];
sx q[1];
rz(-1.55232) q[1];
sx q[1];
rz(1.3943863) q[1];
rz(-2.2398938) q[2];
sx q[2];
rz(-0.82868262) q[2];
sx q[2];
rz(2.5862624) q[2];
rz(-1.3648894) q[3];
sx q[3];
rz(-1.6905224) q[3];
sx q[3];
rz(-0.57991309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
