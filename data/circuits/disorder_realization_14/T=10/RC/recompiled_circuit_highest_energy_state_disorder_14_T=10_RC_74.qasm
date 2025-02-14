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
rz(0.90647107) q[0];
sx q[0];
rz(-1.4426761) q[0];
sx q[0];
rz(0.28191167) q[0];
rz(-2.6126722) q[1];
sx q[1];
rz(-1.6493874) q[1];
sx q[1];
rz(1.563501) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.603005) q[0];
sx q[0];
rz(-1.9360124) q[0];
sx q[0];
rz(2.9486548) q[0];
rz(3.1225948) q[2];
sx q[2];
rz(-0.27488959) q[2];
sx q[2];
rz(1.887448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2927478) q[1];
sx q[1];
rz(-0.94460058) q[1];
sx q[1];
rz(1.4153773) q[1];
rz(0.36564499) q[3];
sx q[3];
rz(-0.72178255) q[3];
sx q[3];
rz(-2.4810098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6866744) q[2];
sx q[2];
rz(-2.8439971) q[2];
sx q[2];
rz(2.5285524) q[2];
rz(2.6679299) q[3];
sx q[3];
rz(-1.1981755) q[3];
sx q[3];
rz(-1.7090428) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7365725) q[0];
sx q[0];
rz(-2.9614145) q[0];
sx q[0];
rz(-0.75102425) q[0];
rz(0.48149064) q[1];
sx q[1];
rz(-1.0844237) q[1];
sx q[1];
rz(0.96985936) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78208047) q[0];
sx q[0];
rz(-2.1954698) q[0];
sx q[0];
rz(2.0332574) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18377797) q[2];
sx q[2];
rz(-2.922195) q[2];
sx q[2];
rz(-2.8156026) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9958933) q[1];
sx q[1];
rz(-0.70317344) q[1];
sx q[1];
rz(-1.9759167) q[1];
rz(-0.59680802) q[3];
sx q[3];
rz(-1.1973901) q[3];
sx q[3];
rz(1.4781836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2227309) q[2];
sx q[2];
rz(-1.6657882) q[2];
sx q[2];
rz(0.53885031) q[2];
rz(-3.1239964) q[3];
sx q[3];
rz(-0.16418695) q[3];
sx q[3];
rz(0.9084475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4410412) q[0];
sx q[0];
rz(-2.7563162) q[0];
sx q[0];
rz(0.061263099) q[0];
rz(-1.6368658) q[1];
sx q[1];
rz(-2.442339) q[1];
sx q[1];
rz(-1.1963074) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.969097) q[0];
sx q[0];
rz(-0.75706702) q[0];
sx q[0];
rz(-1.9226546) q[0];
x q[1];
rz(-0.26232403) q[2];
sx q[2];
rz(-1.6379426) q[2];
sx q[2];
rz(0.1131499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.34431008) q[1];
sx q[1];
rz(-2.4786665) q[1];
sx q[1];
rz(1.0271038) q[1];
rz(-0.34516224) q[3];
sx q[3];
rz(-2.2335839) q[3];
sx q[3];
rz(-1.8927285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4521744) q[2];
sx q[2];
rz(-1.3070561) q[2];
sx q[2];
rz(2.7090731) q[2];
rz(-2.050926) q[3];
sx q[3];
rz(-2.5011823) q[3];
sx q[3];
rz(1.0961016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5525621) q[0];
sx q[0];
rz(-2.2046389) q[0];
sx q[0];
rz(-0.20137782) q[0];
rz(-2.6248113) q[1];
sx q[1];
rz(-0.38025451) q[1];
sx q[1];
rz(0.94863272) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0364591) q[0];
sx q[0];
rz(-2.3984342) q[0];
sx q[0];
rz(1.0418329) q[0];
rz(-2.2791512) q[2];
sx q[2];
rz(-1.7956327) q[2];
sx q[2];
rz(0.98029691) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3740179) q[1];
sx q[1];
rz(-1.5574871) q[1];
sx q[1];
rz(-0.50434169) q[1];
x q[2];
rz(-2.9734083) q[3];
sx q[3];
rz(-0.52709377) q[3];
sx q[3];
rz(0.2179365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2739233) q[2];
sx q[2];
rz(-0.40931585) q[2];
sx q[2];
rz(-2.5841827) q[2];
rz(3.0607767) q[3];
sx q[3];
rz(-1.9010474) q[3];
sx q[3];
rz(-1.935299) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4163365) q[0];
sx q[0];
rz(-1.4755604) q[0];
sx q[0];
rz(-2.1112554) q[0];
rz(0.15580767) q[1];
sx q[1];
rz(-1.067433) q[1];
sx q[1];
rz(-2.9373998) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4725114) q[0];
sx q[0];
rz(-1.8127499) q[0];
sx q[0];
rz(-1.4999963) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7014808) q[2];
sx q[2];
rz(-1.3595264) q[2];
sx q[2];
rz(-1.4071161) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0818644) q[1];
sx q[1];
rz(-1.3646185) q[1];
sx q[1];
rz(2.6463406) q[1];
x q[2];
rz(-1.3315807) q[3];
sx q[3];
rz(-0.81221928) q[3];
sx q[3];
rz(-1.9977457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2121409) q[2];
sx q[2];
rz(-1.6484304) q[2];
sx q[2];
rz(-2.7280651) q[2];
rz(2.2996969) q[3];
sx q[3];
rz(-1.7588408) q[3];
sx q[3];
rz(-2.1690185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095116422) q[0];
sx q[0];
rz(-1.0914509) q[0];
sx q[0];
rz(-3.1360151) q[0];
rz(1.7976409) q[1];
sx q[1];
rz(-0.19897142) q[1];
sx q[1];
rz(-0.70095789) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85421766) q[0];
sx q[0];
rz(-1.1297884) q[0];
sx q[0];
rz(0.99325755) q[0];
rz(-0.92971071) q[2];
sx q[2];
rz(-1.2107009) q[2];
sx q[2];
rz(1.4907229) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.13853874) q[1];
sx q[1];
rz(-1.8459311) q[1];
sx q[1];
rz(1.5281926) q[1];
x q[2];
rz(2.5684729) q[3];
sx q[3];
rz(-1.1240715) q[3];
sx q[3];
rz(-2.5271811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.31412101) q[2];
sx q[2];
rz(-2.8077329) q[2];
sx q[2];
rz(1.1673048) q[2];
rz(2.2589034) q[3];
sx q[3];
rz(-0.76789951) q[3];
sx q[3];
rz(-2.816443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90671396) q[0];
sx q[0];
rz(-2.0218847) q[0];
sx q[0];
rz(0.085302189) q[0];
rz(-2.3872497) q[1];
sx q[1];
rz(-2.6383548) q[1];
sx q[1];
rz(2.1139961) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32815702) q[0];
sx q[0];
rz(-1.3582843) q[0];
sx q[0];
rz(-1.5603609) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.565404) q[2];
sx q[2];
rz(-2.1674334) q[2];
sx q[2];
rz(2.3062458) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4484549) q[1];
sx q[1];
rz(-1.4488646) q[1];
sx q[1];
rz(1.8961402) q[1];
x q[2];
rz(2.1291162) q[3];
sx q[3];
rz(-1.1577679) q[3];
sx q[3];
rz(-2.922204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9895642) q[2];
sx q[2];
rz(-0.98735183) q[2];
sx q[2];
rz(-2.7222471) q[2];
rz(-0.63465214) q[3];
sx q[3];
rz(-2.3366163) q[3];
sx q[3];
rz(2.0161207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3318876) q[0];
sx q[0];
rz(-2.2655847) q[0];
sx q[0];
rz(-2.4494655) q[0];
rz(-0.57811111) q[1];
sx q[1];
rz(-2.7222241) q[1];
sx q[1];
rz(1.3828166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0771187) q[0];
sx q[0];
rz(-0.070264272) q[0];
sx q[0];
rz(1.8514567) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4695669) q[2];
sx q[2];
rz(-1.0254854) q[2];
sx q[2];
rz(1.3915075) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.18305138) q[1];
sx q[1];
rz(-1.019283) q[1];
sx q[1];
rz(0.27490487) q[1];
rz(-pi) q[2];
rz(2.6666497) q[3];
sx q[3];
rz(-0.90427665) q[3];
sx q[3];
rz(2.8996244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77070037) q[2];
sx q[2];
rz(-2.6110677) q[2];
sx q[2];
rz(1.4250866) q[2];
rz(0.51618451) q[3];
sx q[3];
rz(-0.97847146) q[3];
sx q[3];
rz(1.2396575) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6462964) q[0];
sx q[0];
rz(-1.7072059) q[0];
sx q[0];
rz(0.37149757) q[0];
rz(2.3067572) q[1];
sx q[1];
rz(-2.3284262) q[1];
sx q[1];
rz(3.1239948) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69372256) q[0];
sx q[0];
rz(-1.2945064) q[0];
sx q[0];
rz(-1.0280861) q[0];
rz(-pi) q[1];
rz(-0.27522343) q[2];
sx q[2];
rz(-2.38678) q[2];
sx q[2];
rz(1.7856904) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4092651) q[1];
sx q[1];
rz(-2.1346993) q[1];
sx q[1];
rz(-2.0134408) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6952088) q[3];
sx q[3];
rz(-0.65731114) q[3];
sx q[3];
rz(-0.039358649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1060433) q[2];
sx q[2];
rz(-2.273166) q[2];
sx q[2];
rz(3.0412728) q[2];
rz(-2.912168) q[3];
sx q[3];
rz(-1.0474297) q[3];
sx q[3];
rz(0.44982287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4556274) q[0];
sx q[0];
rz(-0.28814155) q[0];
sx q[0];
rz(2.567754) q[0];
rz(-2.0876743) q[1];
sx q[1];
rz(-2.0951447) q[1];
sx q[1];
rz(-2.4651249) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5887816) q[0];
sx q[0];
rz(-1.2924606) q[0];
sx q[0];
rz(1.8235444) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25679882) q[2];
sx q[2];
rz(-1.9635824) q[2];
sx q[2];
rz(1.2744255) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0770196) q[1];
sx q[1];
rz(-1.4823682) q[1];
sx q[1];
rz(2.6131373) q[1];
rz(-2.9432137) q[3];
sx q[3];
rz(-1.0470445) q[3];
sx q[3];
rz(-0.74568053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4544868) q[2];
sx q[2];
rz(-0.7578308) q[2];
sx q[2];
rz(1.6472316) q[2];
rz(-2.2729661) q[3];
sx q[3];
rz(-0.70610154) q[3];
sx q[3];
rz(2.6715265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070521991) q[0];
sx q[0];
rz(-1.1916397) q[0];
sx q[0];
rz(1.2217039) q[0];
rz(-1.8078049) q[1];
sx q[1];
rz(-1.318327) q[1];
sx q[1];
rz(-0.64073906) q[1];
rz(-2.376198) q[2];
sx q[2];
rz(-1.5049878) q[2];
sx q[2];
rz(2.4743248) q[2];
rz(1.7464433) q[3];
sx q[3];
rz(-1.5271389) q[3];
sx q[3];
rz(2.3300119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
