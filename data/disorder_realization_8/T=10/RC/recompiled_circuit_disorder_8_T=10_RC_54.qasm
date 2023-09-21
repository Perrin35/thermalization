OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8060057) q[0];
sx q[0];
rz(-0.94526362) q[0];
sx q[0];
rz(2.6160016) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(-0.90484172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5072767) q[0];
sx q[0];
rz(-2.0024558) q[0];
sx q[0];
rz(-1.0213724) q[0];
x q[1];
rz(2.9154539) q[2];
sx q[2];
rz(-1.7625426) q[2];
sx q[2];
rz(2.6297671) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9368254) q[1];
sx q[1];
rz(-1.966404) q[1];
sx q[1];
rz(-0.25524615) q[1];
x q[2];
rz(-2.9526887) q[3];
sx q[3];
rz(-1.8763262) q[3];
sx q[3];
rz(-2.4019935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9378172) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(-0.096244372) q[2];
rz(-2.105666) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(0.15371418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44089833) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(-2.3764215) q[0];
rz(-1.2922497) q[1];
sx q[1];
rz(-2.6563829) q[1];
sx q[1];
rz(0.66295019) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0239319) q[0];
sx q[0];
rz(-1.5045325) q[0];
sx q[0];
rz(-0.06225417) q[0];
rz(2.2234369) q[2];
sx q[2];
rz(-1.895972) q[2];
sx q[2];
rz(0.52415372) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6389097) q[1];
sx q[1];
rz(-2.0109634) q[1];
sx q[1];
rz(1.6619976) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8421721) q[3];
sx q[3];
rz(-0.50588183) q[3];
sx q[3];
rz(-0.29380709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.092676) q[2];
sx q[2];
rz(-2.0930223) q[2];
sx q[2];
rz(0.32901397) q[2];
rz(-2.4760903) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(-1.8238508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5927785) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(-2.9192525) q[0];
rz(1.0173343) q[1];
sx q[1];
rz(-2.4203114) q[1];
sx q[1];
rz(0.51868784) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.737239) q[0];
sx q[0];
rz(-0.81439942) q[0];
sx q[0];
rz(1.0712207) q[0];
rz(-1.6665886) q[2];
sx q[2];
rz(-1.3057858) q[2];
sx q[2];
rz(-1.4694627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.31836244) q[1];
sx q[1];
rz(-1.743411) q[1];
sx q[1];
rz(-0.77962064) q[1];
x q[2];
rz(1.5716001) q[3];
sx q[3];
rz(-1.5161848) q[3];
sx q[3];
rz(0.87517232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2805933) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(0.60244256) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(-3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4908726) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(2.9673476) q[0];
rz(-0.53025591) q[1];
sx q[1];
rz(-1.6590051) q[1];
sx q[1];
rz(2.6285016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4671191) q[0];
sx q[0];
rz(-1.4917443) q[0];
sx q[0];
rz(0.62069686) q[0];
rz(-pi) q[1];
rz(-0.80715837) q[2];
sx q[2];
rz(-2.0370738) q[2];
sx q[2];
rz(1.5392787) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0780371) q[1];
sx q[1];
rz(-2.6647898) q[1];
sx q[1];
rz(2.8366691) q[1];
rz(1.3644049) q[3];
sx q[3];
rz(-1.6009814) q[3];
sx q[3];
rz(0.15011945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6667368) q[2];
sx q[2];
rz(-0.36617827) q[2];
sx q[2];
rz(-0.22988698) q[2];
rz(2.722548) q[3];
sx q[3];
rz(-1.7925526) q[3];
sx q[3];
rz(-2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4693562) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(-2.4647734) q[0];
rz(0.49304402) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(-0.61606032) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0115833) q[0];
sx q[0];
rz(-1.6066215) q[0];
sx q[0];
rz(-1.623276) q[0];
rz(1.6904171) q[2];
sx q[2];
rz(-1.9336485) q[2];
sx q[2];
rz(2.5809443) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38356009) q[1];
sx q[1];
rz(-0.57764232) q[1];
sx q[1];
rz(-3.1314965) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2474399) q[3];
sx q[3];
rz(-1.3291306) q[3];
sx q[3];
rz(-1.5622996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4074576) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(-0.25536728) q[2];
rz(-1.6051965) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(2.3674964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42246321) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(-0.47250026) q[0];
rz(0.52608144) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(0.88476673) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038267604) q[0];
sx q[0];
rz(-1.190435) q[0];
sx q[0];
rz(0.079770712) q[0];
rz(1.556649) q[2];
sx q[2];
rz(-1.3607585) q[2];
sx q[2];
rz(-1.998979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6217475) q[1];
sx q[1];
rz(-1.8583082) q[1];
sx q[1];
rz(1.1430172) q[1];
rz(-pi) q[2];
rz(-1.1690508) q[3];
sx q[3];
rz(-2.5707977) q[3];
sx q[3];
rz(-1.4626383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3970967) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(-0.3113783) q[2];
rz(1.3686251) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(-0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235274) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(-3.080522) q[0];
rz(-0.04018499) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(0.73289245) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7948579) q[0];
sx q[0];
rz(-1.1746527) q[0];
sx q[0];
rz(-2.5651155) q[0];
rz(-pi) q[1];
x q[1];
rz(2.272846) q[2];
sx q[2];
rz(-1.1940496) q[2];
sx q[2];
rz(1.6422611) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.43436189) q[1];
sx q[1];
rz(-1.5642628) q[1];
sx q[1];
rz(-2.1965501) q[1];
x q[2];
rz(1.93768) q[3];
sx q[3];
rz(-0.11493472) q[3];
sx q[3];
rz(-0.68655187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0682893) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(-2.627009) q[2];
rz(-1.2375281) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056203689) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(-2.432166) q[0];
rz(1.6363232) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(-0.27871305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.404971) q[0];
sx q[0];
rz(-1.4355037) q[0];
sx q[0];
rz(2.9297329) q[0];
rz(0.3785554) q[2];
sx q[2];
rz(-2.5494529) q[2];
sx q[2];
rz(0.90422599) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7396001) q[1];
sx q[1];
rz(-0.38312373) q[1];
sx q[1];
rz(0.95462228) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9023499) q[3];
sx q[3];
rz(-2.9058876) q[3];
sx q[3];
rz(1.2329668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8401106) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(-2.2864443) q[3];
sx q[3];
rz(-2.6931098) q[3];
sx q[3];
rz(0.40518951) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99496019) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(-0.96889281) q[0];
rz(2.6682207) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(1.999058) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5794967) q[0];
sx q[0];
rz(-1.8651433) q[0];
sx q[0];
rz(-0.094898183) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4382486) q[2];
sx q[2];
rz(-1.5726657) q[2];
sx q[2];
rz(1.6714931) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7361703) q[1];
sx q[1];
rz(-2.1076964) q[1];
sx q[1];
rz(1.7169397) q[1];
x q[2];
rz(-2.9386018) q[3];
sx q[3];
rz(-0.91068017) q[3];
sx q[3];
rz(-2.4329894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.896495) q[2];
sx q[2];
rz(-1.921804) q[2];
sx q[2];
rz(1.0207821) q[2];
rz(-2.8178689) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.1616515) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(-0.7014057) q[0];
rz(0.91570634) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(-1.9030301) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68574821) q[0];
sx q[0];
rz(-1.0554753) q[0];
sx q[0];
rz(-2.8932163) q[0];
rz(0.054758666) q[2];
sx q[2];
rz(-2.3878532) q[2];
sx q[2];
rz(-2.4278305) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.726767) q[1];
sx q[1];
rz(-1.6884202) q[1];
sx q[1];
rz(-0.91314258) q[1];
x q[2];
rz(-0.11390399) q[3];
sx q[3];
rz(-0.98131991) q[3];
sx q[3];
rz(-1.9437499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68676585) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(-1.94058) q[3];
sx q[3];
rz(-2.4062556) q[3];
sx q[3];
rz(-2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6702406) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(-0.025370601) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(1.17169) q[2];
sx q[2];
rz(-1.8532955) q[2];
sx q[2];
rz(2.4886139) q[2];
rz(1.3958037) q[3];
sx q[3];
rz(-2.3621515) q[3];
sx q[3];
rz(3.0099517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
