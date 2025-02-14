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
rz(0.60740745) q[0];
sx q[0];
rz(2.6976801) q[0];
sx q[0];
rz(9.1023268) q[0];
rz(-1.9880265) q[1];
sx q[1];
rz(-1.1918951) q[1];
sx q[1];
rz(-2.9236887) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76989749) q[0];
sx q[0];
rz(-1.8718693) q[0];
sx q[0];
rz(0.45905827) q[0];
rz(-pi) q[1];
x q[1];
rz(1.409265) q[2];
sx q[2];
rz(-1.7551975) q[2];
sx q[2];
rz(-0.34689515) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3609546) q[1];
sx q[1];
rz(-1.8285969) q[1];
sx q[1];
rz(1.4514489) q[1];
rz(-pi) q[2];
rz(2.3090906) q[3];
sx q[3];
rz(-0.41838405) q[3];
sx q[3];
rz(2.7184021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.13088642) q[2];
sx q[2];
rz(-2.71038) q[2];
sx q[2];
rz(-0.19403379) q[2];
rz(0.38328299) q[3];
sx q[3];
rz(-2.2563939) q[3];
sx q[3];
rz(1.5017728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5633157) q[0];
sx q[0];
rz(-1.0551772) q[0];
sx q[0];
rz(-1.1042327) q[0];
rz(2.4597994) q[1];
sx q[1];
rz(-1.7742523) q[1];
sx q[1];
rz(-1.7368447) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2629614) q[0];
sx q[0];
rz(-1.6978953) q[0];
sx q[0];
rz(-0.2574347) q[0];
rz(-0.67595903) q[2];
sx q[2];
rz(-1.1164719) q[2];
sx q[2];
rz(-0.76269645) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8190143) q[1];
sx q[1];
rz(-2.0657263) q[1];
sx q[1];
rz(-1.9274345) q[1];
rz(-pi) q[2];
rz(-1.4561171) q[3];
sx q[3];
rz(-1.3175829) q[3];
sx q[3];
rz(2.6270783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3746987) q[2];
sx q[2];
rz(-2.1465492) q[2];
sx q[2];
rz(1.9083171) q[2];
rz(0.79902664) q[3];
sx q[3];
rz(-2.7965751) q[3];
sx q[3];
rz(-0.12921216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98053539) q[0];
sx q[0];
rz(-2.0687456) q[0];
sx q[0];
rz(2.5824353) q[0];
rz(2.0237538) q[1];
sx q[1];
rz(-1.585377) q[1];
sx q[1];
rz(0.1112172) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98830279) q[0];
sx q[0];
rz(-2.7493101) q[0];
sx q[0];
rz(2.8610703) q[0];
rz(-pi) q[1];
rz(2.5170277) q[2];
sx q[2];
rz(-1.7918158) q[2];
sx q[2];
rz(2.2159037) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8168252) q[1];
sx q[1];
rz(-2.2681464) q[1];
sx q[1];
rz(0.24705418) q[1];
rz(-pi) q[2];
rz(1.4239156) q[3];
sx q[3];
rz(-1.1814756) q[3];
sx q[3];
rz(1.9922272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3126276) q[2];
sx q[2];
rz(-1.6964922) q[2];
sx q[2];
rz(-1.9435389) q[2];
rz(1.4231921) q[3];
sx q[3];
rz(-1.6953902) q[3];
sx q[3];
rz(0.47553441) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6296185) q[0];
sx q[0];
rz(-0.29933512) q[0];
sx q[0];
rz(-0.83628118) q[0];
rz(-1.5760999) q[1];
sx q[1];
rz(-1.0944159) q[1];
sx q[1];
rz(2.2907168) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.068965) q[0];
sx q[0];
rz(-0.57372813) q[0];
sx q[0];
rz(-1.5092586) q[0];
rz(-2.3281872) q[2];
sx q[2];
rz(-2.4508173) q[2];
sx q[2];
rz(-1.9936313) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2712114) q[1];
sx q[1];
rz(-0.79194259) q[1];
sx q[1];
rz(-0.19510896) q[1];
rz(-pi) q[2];
rz(-0.20765813) q[3];
sx q[3];
rz(-1.1829783) q[3];
sx q[3];
rz(0.70630276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1026844) q[2];
sx q[2];
rz(-2.2816198) q[2];
sx q[2];
rz(0.50722185) q[2];
rz(-2.475259) q[3];
sx q[3];
rz(-2.3907876) q[3];
sx q[3];
rz(-0.91608086) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12793334) q[0];
sx q[0];
rz(-1.6935231) q[0];
sx q[0];
rz(1.5284982) q[0];
rz(-1.5616034) q[1];
sx q[1];
rz(-1.0630307) q[1];
sx q[1];
rz(-2.9772421) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2737758) q[0];
sx q[0];
rz(-1.4901596) q[0];
sx q[0];
rz(-0.095726526) q[0];
rz(-pi) q[1];
rz(-0.88824026) q[2];
sx q[2];
rz(-1.054371) q[2];
sx q[2];
rz(2.7009855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1406255) q[1];
sx q[1];
rz(-0.78636175) q[1];
sx q[1];
rz(-1.0746075) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2772395) q[3];
sx q[3];
rz(-1.3428215) q[3];
sx q[3];
rz(-0.10693947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5932172) q[2];
sx q[2];
rz(-2.2943353) q[2];
sx q[2];
rz(2.4753921) q[2];
rz(-0.14144746) q[3];
sx q[3];
rz(-1.5902767) q[3];
sx q[3];
rz(0.9945873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87404609) q[0];
sx q[0];
rz(-0.62804896) q[0];
sx q[0];
rz(-2.6603267) q[0];
rz(0.73356837) q[1];
sx q[1];
rz(-1.8736519) q[1];
sx q[1];
rz(3.0487294) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32321445) q[0];
sx q[0];
rz(-2.0131454) q[0];
sx q[0];
rz(1.044892) q[0];
x q[1];
rz(3.0888482) q[2];
sx q[2];
rz(-3.102902) q[2];
sx q[2];
rz(0.82088137) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9290849) q[1];
sx q[1];
rz(-1.9783535) q[1];
sx q[1];
rz(-2.4465438) q[1];
rz(-pi) q[2];
rz(-2.8792297) q[3];
sx q[3];
rz(-1.2195346) q[3];
sx q[3];
rz(-2.2134804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4197454) q[2];
sx q[2];
rz(-1.5200619) q[2];
sx q[2];
rz(-1.2255555) q[2];
rz(-0.078977481) q[3];
sx q[3];
rz(-2.0073399) q[3];
sx q[3];
rz(1.8592161) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52800286) q[0];
sx q[0];
rz(-0.71857518) q[0];
sx q[0];
rz(2.3336616) q[0];
rz(1.5241874) q[1];
sx q[1];
rz(-0.15847358) q[1];
sx q[1];
rz(2.4375516) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8600143) q[0];
sx q[0];
rz(-1.5509971) q[0];
sx q[0];
rz(0.0088282882) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1374409) q[2];
sx q[2];
rz(-1.5791681) q[2];
sx q[2];
rz(-0.75941759) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0800533) q[1];
sx q[1];
rz(-0.94990094) q[1];
sx q[1];
rz(3.1136334) q[1];
x q[2];
rz(-2.9765997) q[3];
sx q[3];
rz(-2.1793302) q[3];
sx q[3];
rz(-0.42770619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0713221) q[2];
sx q[2];
rz(-2.4532048) q[2];
sx q[2];
rz(2.2570611) q[2];
rz(0.5136579) q[3];
sx q[3];
rz(-1.1762041) q[3];
sx q[3];
rz(2.6562712) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270444) q[0];
sx q[0];
rz(-0.7651279) q[0];
sx q[0];
rz(2.2807518) q[0];
rz(-0.035482081) q[1];
sx q[1];
rz(-1.2018964) q[1];
sx q[1];
rz(-1.4883581) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5527715) q[0];
sx q[0];
rz(-2.5591922) q[0];
sx q[0];
rz(0.78314535) q[0];
rz(-pi) q[1];
rz(-2.3242438) q[2];
sx q[2];
rz(-0.32830445) q[2];
sx q[2];
rz(0.71336929) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70362857) q[1];
sx q[1];
rz(-1.1878106) q[1];
sx q[1];
rz(0.18937892) q[1];
rz(-0.32522301) q[3];
sx q[3];
rz(-2.8517493) q[3];
sx q[3];
rz(-1.0018283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.615768) q[2];
sx q[2];
rz(-1.1349698) q[2];
sx q[2];
rz(1.5110678) q[2];
rz(-0.75119558) q[3];
sx q[3];
rz(-2.5978751) q[3];
sx q[3];
rz(0.80120075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57665956) q[0];
sx q[0];
rz(-1.7486005) q[0];
sx q[0];
rz(0.15889731) q[0];
rz(-1.502602) q[1];
sx q[1];
rz(-1.8105806) q[1];
sx q[1];
rz(1.9749036) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1142857) q[0];
sx q[0];
rz(-2.5059359) q[0];
sx q[0];
rz(-0.42814769) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8620637) q[2];
sx q[2];
rz(-1.88481) q[2];
sx q[2];
rz(-1.2316126) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9433971) q[1];
sx q[1];
rz(-1.5540197) q[1];
sx q[1];
rz(1.6086964) q[1];
rz(-pi) q[2];
rz(-0.95440475) q[3];
sx q[3];
rz(-0.036374854) q[3];
sx q[3];
rz(-2.006435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73929536) q[2];
sx q[2];
rz(-0.6604971) q[2];
sx q[2];
rz(2.8759549) q[2];
rz(1.9499251) q[3];
sx q[3];
rz(-1.8040413) q[3];
sx q[3];
rz(0.55408293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83248508) q[0];
sx q[0];
rz(-2.5256248) q[0];
sx q[0];
rz(-0.68514222) q[0];
rz(2.5849672) q[1];
sx q[1];
rz(-1.4645422) q[1];
sx q[1];
rz(2.305078) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13930411) q[0];
sx q[0];
rz(-1.909585) q[0];
sx q[0];
rz(2.0071507) q[0];
rz(-pi) q[1];
rz(-0.30843251) q[2];
sx q[2];
rz(-1.0717725) q[2];
sx q[2];
rz(-2.4633046) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.54933465) q[1];
sx q[1];
rz(-0.95227949) q[1];
sx q[1];
rz(1.0299512) q[1];
x q[2];
rz(0.1565069) q[3];
sx q[3];
rz(-1.8409074) q[3];
sx q[3];
rz(0.47299415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4461925) q[2];
sx q[2];
rz(-1.4245028) q[2];
sx q[2];
rz(0.85644537) q[2];
rz(-0.19927464) q[3];
sx q[3];
rz(-2.1592185) q[3];
sx q[3];
rz(2.2299531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0675426) q[0];
sx q[0];
rz(-1.0150801) q[0];
sx q[0];
rz(-1.6652921) q[0];
rz(0.53786565) q[1];
sx q[1];
rz(-1.270351) q[1];
sx q[1];
rz(-1.7958633) q[1];
rz(-0.79277586) q[2];
sx q[2];
rz(-1.0083099) q[2];
sx q[2];
rz(-0.56854962) q[2];
rz(1.8798697) q[3];
sx q[3];
rz(-2.3936987) q[3];
sx q[3];
rz(1.0868418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
