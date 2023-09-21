OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(-1.3347081) q[0];
rz(2.788738) q[1];
sx q[1];
rz(-2.9810413) q[1];
sx q[1];
rz(-0.97595739) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87603509) q[0];
sx q[0];
rz(-1.3296488) q[0];
sx q[0];
rz(-0.59224706) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8085254) q[2];
sx q[2];
rz(-1.0330079) q[2];
sx q[2];
rz(-1.2531467) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.576697) q[1];
sx q[1];
rz(-1.7045867) q[1];
sx q[1];
rz(-2.1810075) q[1];
rz(-pi) q[2];
rz(-0.8043886) q[3];
sx q[3];
rz(-1.8612923) q[3];
sx q[3];
rz(2.1936072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41539899) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(2.3577918) q[2];
rz(2.8090254) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(1.3403085) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587104) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(0.17856199) q[0];
rz(1.3372955) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(0.006342412) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77942383) q[0];
sx q[0];
rz(-2.9075025) q[0];
sx q[0];
rz(-2.1650044) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6132934) q[2];
sx q[2];
rz(-2.106278) q[2];
sx q[2];
rz(-2.2250125) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.983922) q[1];
sx q[1];
rz(-1.5092641) q[1];
sx q[1];
rz(-2.8308949) q[1];
x q[2];
rz(2.5793521) q[3];
sx q[3];
rz(-1.7842818) q[3];
sx q[3];
rz(1.436304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.94770849) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(0.87810278) q[2];
rz(0.86205035) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5851615) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(-0.01097824) q[0];
rz(0.36704656) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(3.045851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84155267) q[0];
sx q[0];
rz(-1.3944355) q[0];
sx q[0];
rz(0.15533133) q[0];
rz(-pi) q[1];
rz(2.0289621) q[2];
sx q[2];
rz(-1.1095699) q[2];
sx q[2];
rz(-2.0756276) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.84903753) q[1];
sx q[1];
rz(-1.109086) q[1];
sx q[1];
rz(-2.1828116) q[1];
rz(0.5511958) q[3];
sx q[3];
rz(-1.5407729) q[3];
sx q[3];
rz(-2.4195645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0720955) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(-2.3068008) q[2];
rz(-0.21162027) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(-0.37477469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9469706) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(-0.34657493) q[0];
rz(-0.52571458) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(-2.1077164) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945275) q[0];
sx q[0];
rz(-0.80104242) q[0];
sx q[0];
rz(0.69181504) q[0];
rz(-0.526555) q[2];
sx q[2];
rz(-2.0806081) q[2];
sx q[2];
rz(-1.1916135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4775866) q[1];
sx q[1];
rz(-2.9395182) q[1];
sx q[1];
rz(1.9007773) q[1];
x q[2];
rz(-0.28169607) q[3];
sx q[3];
rz(-1.2558189) q[3];
sx q[3];
rz(-1.7236934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6870849) q[2];
sx q[2];
rz(-1.4322174) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(2.7205617) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(-0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43276697) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(-2.72686) q[0];
rz(1.746009) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(2.5674852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60526472) q[0];
sx q[0];
rz(-2.5308373) q[0];
sx q[0];
rz(1.3461793) q[0];
rz(-pi) q[1];
rz(-1.5933883) q[2];
sx q[2];
rz(-2.2859757) q[2];
sx q[2];
rz(1.1720282) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.638006) q[1];
sx q[1];
rz(-1.2212911) q[1];
sx q[1];
rz(-1.867678) q[1];
rz(0.19458171) q[3];
sx q[3];
rz(-2.2268725) q[3];
sx q[3];
rz(1.7723099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85270143) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(-1.627702) q[2];
rz(0.04143516) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(-3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067327499) q[0];
sx q[0];
rz(-1.3621618) q[0];
sx q[0];
rz(2.4434027) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(0.73227698) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4740144) q[0];
sx q[0];
rz(-1.3587553) q[0];
sx q[0];
rz(1.7330806) q[0];
rz(-pi) q[1];
rz(0.68533021) q[2];
sx q[2];
rz(-2.2517859) q[2];
sx q[2];
rz(0.53390098) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1440891) q[1];
sx q[1];
rz(-2.2170076) q[1];
sx q[1];
rz(-1.8540107) q[1];
rz(-pi) q[2];
x q[2];
rz(0.332571) q[3];
sx q[3];
rz(-2.5453574) q[3];
sx q[3];
rz(2.3308844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3884864) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(-1.9801271) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(-1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.5724065) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(2.9851595) q[0];
rz(-0.45174831) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(0.77004534) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2419287) q[0];
sx q[0];
rz(-2.1673492) q[0];
sx q[0];
rz(-1.3616189) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1427878) q[2];
sx q[2];
rz(-0.6711798) q[2];
sx q[2];
rz(0.83516781) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3502096) q[1];
sx q[1];
rz(-2.0767127) q[1];
sx q[1];
rz(-2.4398068) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6386912) q[3];
sx q[3];
rz(-1.8637878) q[3];
sx q[3];
rz(-2.6754232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6440755) q[2];
sx q[2];
rz(-0.13921177) q[2];
sx q[2];
rz(-1.0151803) q[2];
rz(-1.4366359) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(0.59593433) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5148233) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(0.24169895) q[0];
rz(-2.4027951) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(2.7752005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05224932) q[0];
sx q[0];
rz(-1.0567259) q[0];
sx q[0];
rz(-0.71787562) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56717746) q[2];
sx q[2];
rz(-0.71225538) q[2];
sx q[2];
rz(-0.30356193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.061325039) q[1];
sx q[1];
rz(-0.72239164) q[1];
sx q[1];
rz(-2.8894043) q[1];
rz(2.1226235) q[3];
sx q[3];
rz(-1.6120211) q[3];
sx q[3];
rz(0.036389694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1239132) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(0.29385847) q[2];
rz(-3.1270694) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(-0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095187) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(-0.88395399) q[0];
rz(2.4813095) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(-2.3502137) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8732849) q[0];
sx q[0];
rz(-2.377016) q[0];
sx q[0];
rz(-1.1029878) q[0];
x q[1];
rz(-0.4653761) q[2];
sx q[2];
rz(-2.010979) q[2];
sx q[2];
rz(2.344775) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2453354) q[1];
sx q[1];
rz(-0.63071139) q[1];
sx q[1];
rz(-0.74951042) q[1];
x q[2];
rz(2.1095554) q[3];
sx q[3];
rz(-1.7146535) q[3];
sx q[3];
rz(2.0335576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.066594921) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(-0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(-0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026697712) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(-2.5337906) q[0];
rz(0.23884493) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(-1.3806608) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5138807) q[0];
sx q[0];
rz(-2.4988334) q[0];
sx q[0];
rz(1.7895262) q[0];
rz(-3.1157687) q[2];
sx q[2];
rz(-2.0061473) q[2];
sx q[2];
rz(-0.045407427) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7271991) q[1];
sx q[1];
rz(-2.6135923) q[1];
sx q[1];
rz(2.2382733) q[1];
rz(-pi) q[2];
rz(1.4019743) q[3];
sx q[3];
rz(-1.7572548) q[3];
sx q[3];
rz(0.25690119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3537139) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-2.805368) q[2];
rz(-2.0119038) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(-2.0000134) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0040141) q[0];
sx q[0];
rz(-1.4083569) q[0];
sx q[0];
rz(1.2271723) q[0];
rz(-2.5972988) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(0.33967321) q[2];
sx q[2];
rz(-0.85396955) q[2];
sx q[2];
rz(0.11243482) q[2];
rz(-2.9059698) q[3];
sx q[3];
rz(-1.0141254) q[3];
sx q[3];
rz(0.49401382) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
