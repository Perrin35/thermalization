OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0394734) q[0];
sx q[0];
rz(-1.4730299) q[0];
sx q[0];
rz(0.13134512) q[0];
rz(-3.1047473) q[1];
sx q[1];
rz(3.6689833) q[1];
sx q[1];
rz(13.875516) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86517559) q[0];
sx q[0];
rz(-2.2593479) q[0];
sx q[0];
rz(0.49746969) q[0];
rz(-pi) q[1];
rz(-2.1151343) q[2];
sx q[2];
rz(-1.9690459) q[2];
sx q[2];
rz(-1.2521225) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.943191) q[1];
sx q[1];
rz(-1.3326338) q[1];
sx q[1];
rz(-2.8203301) q[1];
rz(-pi) q[2];
rz(1.3996436) q[3];
sx q[3];
rz(-1.389909) q[3];
sx q[3];
rz(0.76598841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6068136) q[2];
sx q[2];
rz(-2.7853192) q[2];
sx q[2];
rz(-0.66205364) q[2];
rz(-2.3946297) q[3];
sx q[3];
rz(-2.2253939) q[3];
sx q[3];
rz(1.0158739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58716431) q[0];
sx q[0];
rz(-2.9968408) q[0];
sx q[0];
rz(1.1821049) q[0];
rz(2.3960522) q[1];
sx q[1];
rz(-0.59919557) q[1];
sx q[1];
rz(0.10339698) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1348159) q[0];
sx q[0];
rz(-0.49044436) q[0];
sx q[0];
rz(2.6314526) q[0];
x q[1];
rz(2.75704) q[2];
sx q[2];
rz(-2.1332624) q[2];
sx q[2];
rz(0.99316521) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.34161257) q[1];
sx q[1];
rz(-2.3368521) q[1];
sx q[1];
rz(2.0384203) q[1];
rz(-pi) q[2];
rz(2.313226) q[3];
sx q[3];
rz(-1.4423372) q[3];
sx q[3];
rz(1.0297694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5891002) q[2];
sx q[2];
rz(-1.8742259) q[2];
sx q[2];
rz(0.44431552) q[2];
rz(-2.0878504) q[3];
sx q[3];
rz(-0.070662347) q[3];
sx q[3];
rz(1.1963371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.062773) q[0];
sx q[0];
rz(-1.2508996) q[0];
sx q[0];
rz(-0.66251063) q[0];
rz(-1.484681) q[1];
sx q[1];
rz(-2.1187481) q[1];
sx q[1];
rz(-2.9248765) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0652567) q[0];
sx q[0];
rz(-0.88533339) q[0];
sx q[0];
rz(-0.61409593) q[0];
rz(2.4091085) q[2];
sx q[2];
rz(-1.9062796) q[2];
sx q[2];
rz(-1.6215289) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2487434) q[1];
sx q[1];
rz(-2.4676552) q[1];
sx q[1];
rz(-1.6350063) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6469798) q[3];
sx q[3];
rz(-1.7124767) q[3];
sx q[3];
rz(2.6050267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4776939) q[2];
sx q[2];
rz(-0.52637664) q[2];
sx q[2];
rz(-2.9065175) q[2];
rz(-2.7663686) q[3];
sx q[3];
rz(-2.2347361) q[3];
sx q[3];
rz(-2.4231329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5821871) q[0];
sx q[0];
rz(-0.087787293) q[0];
sx q[0];
rz(-1.0002332) q[0];
rz(1.2031215) q[1];
sx q[1];
rz(-1.5088046) q[1];
sx q[1];
rz(-2.5968754) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9207912) q[0];
sx q[0];
rz(-1.5259597) q[0];
sx q[0];
rz(-0.71323552) q[0];
rz(-2.4248872) q[2];
sx q[2];
rz(-0.96895987) q[2];
sx q[2];
rz(0.81317893) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.57883731) q[1];
sx q[1];
rz(-1.75995) q[1];
sx q[1];
rz(3.021043) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8417743) q[3];
sx q[3];
rz(-2.162419) q[3];
sx q[3];
rz(-0.064379582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41877052) q[2];
sx q[2];
rz(-2.3493769) q[2];
sx q[2];
rz(-2.344632) q[2];
rz(-2.8523417) q[3];
sx q[3];
rz(-0.97024337) q[3];
sx q[3];
rz(0.42118916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5274984) q[0];
sx q[0];
rz(-0.088936381) q[0];
sx q[0];
rz(-1.8726789) q[0];
rz(2.1753963) q[1];
sx q[1];
rz(-1.7485917) q[1];
sx q[1];
rz(0.46636811) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.120935) q[0];
sx q[0];
rz(-1.0758721) q[0];
sx q[0];
rz(1.0878956) q[0];
x q[1];
rz(-0.63299243) q[2];
sx q[2];
rz(-2.5883753) q[2];
sx q[2];
rz(-0.67942441) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1524017) q[1];
sx q[1];
rz(-0.70211239) q[1];
sx q[1];
rz(1.0662119) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2339646) q[3];
sx q[3];
rz(-1.8844127) q[3];
sx q[3];
rz(-0.37330353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40631488) q[2];
sx q[2];
rz(-2.3390528) q[2];
sx q[2];
rz(-0.80424133) q[2];
rz(2.6546226) q[3];
sx q[3];
rz(-2.099497) q[3];
sx q[3];
rz(3.13412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93953472) q[0];
sx q[0];
rz(-2.845293) q[0];
sx q[0];
rz(2.1837088) q[0];
rz(1.5127888) q[1];
sx q[1];
rz(-2.084338) q[1];
sx q[1];
rz(1.588795) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5361547) q[0];
sx q[0];
rz(-3.0585318) q[0];
sx q[0];
rz(-1.2916748) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27811173) q[2];
sx q[2];
rz(-1.8459003) q[2];
sx q[2];
rz(-2.846707) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54974906) q[1];
sx q[1];
rz(-1.7477682) q[1];
sx q[1];
rz(-0.66284499) q[1];
x q[2];
rz(2.5143322) q[3];
sx q[3];
rz(-1.2377487) q[3];
sx q[3];
rz(1.8740538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2019041) q[2];
sx q[2];
rz(-1.0558015) q[2];
sx q[2];
rz(2.4064257) q[2];
rz(2.1926664) q[3];
sx q[3];
rz(-2.2831423) q[3];
sx q[3];
rz(2.2775876) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3801124) q[0];
sx q[0];
rz(-1.004847) q[0];
sx q[0];
rz(2.5349706) q[0];
rz(-2.3573549) q[1];
sx q[1];
rz(-1.8083068) q[1];
sx q[1];
rz(-1.3444791) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79216751) q[0];
sx q[0];
rz(-1.0972404) q[0];
sx q[0];
rz(2.4132632) q[0];
rz(2.6980618) q[2];
sx q[2];
rz(-0.54033989) q[2];
sx q[2];
rz(0.72393723) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.575702) q[1];
sx q[1];
rz(-1.9676349) q[1];
sx q[1];
rz(1.8116519) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7502039) q[3];
sx q[3];
rz(-0.47166892) q[3];
sx q[3];
rz(-2.6149132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6091696) q[2];
sx q[2];
rz(-0.50749856) q[2];
sx q[2];
rz(1.7519105) q[2];
rz(-0.17942795) q[3];
sx q[3];
rz(-0.98589412) q[3];
sx q[3];
rz(-2.7348886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0591902) q[0];
sx q[0];
rz(-2.817509) q[0];
sx q[0];
rz(-2.3865336) q[0];
rz(0.45626196) q[1];
sx q[1];
rz(-1.4249964) q[1];
sx q[1];
rz(1.0770575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2307325) q[0];
sx q[0];
rz(-2.5344779) q[0];
sx q[0];
rz(-0.40672983) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12142373) q[2];
sx q[2];
rz(-2.3540263) q[2];
sx q[2];
rz(1.9483669) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39246074) q[1];
sx q[1];
rz(-1.8549671) q[1];
sx q[1];
rz(-0.47537132) q[1];
x q[2];
rz(-1.3102688) q[3];
sx q[3];
rz(-0.55580492) q[3];
sx q[3];
rz(2.1422841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.72851744) q[2];
sx q[2];
rz(-2.04144) q[2];
sx q[2];
rz(-2.4033578) q[2];
rz(-1.5089367) q[3];
sx q[3];
rz(-1.3239219) q[3];
sx q[3];
rz(1.9807321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4974834) q[0];
sx q[0];
rz(-1.4583541) q[0];
sx q[0];
rz(0.64895502) q[0];
rz(-0.68967462) q[1];
sx q[1];
rz(-0.70976218) q[1];
sx q[1];
rz(-0.41742691) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2608503) q[0];
sx q[0];
rz(-2.1584903) q[0];
sx q[0];
rz(1.2686612) q[0];
rz(2.8418188) q[2];
sx q[2];
rz(-1.6437141) q[2];
sx q[2];
rz(-2.656183) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6747804) q[1];
sx q[1];
rz(-1.1622475) q[1];
sx q[1];
rz(1.4571587) q[1];
rz(-1.9661918) q[3];
sx q[3];
rz(-1.1987276) q[3];
sx q[3];
rz(1.7898066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6680341) q[2];
sx q[2];
rz(-1.7201951) q[2];
sx q[2];
rz(-0.045698015) q[2];
rz(-0.33347305) q[3];
sx q[3];
rz(-2.5662751) q[3];
sx q[3];
rz(-3.0157109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4561653) q[0];
sx q[0];
rz(-0.5394772) q[0];
sx q[0];
rz(1.9239377) q[0];
rz(2.894891) q[1];
sx q[1];
rz(-1.291899) q[1];
sx q[1];
rz(-2.6620679) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3772904) q[0];
sx q[0];
rz(-0.78702292) q[0];
sx q[0];
rz(-2.2020409) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2662451) q[2];
sx q[2];
rz(-1.9672868) q[2];
sx q[2];
rz(-3.0975395) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3563167) q[1];
sx q[1];
rz(-2.0802167) q[1];
sx q[1];
rz(1.9059577) q[1];
x q[2];
rz(2.6761618) q[3];
sx q[3];
rz(-1.7836535) q[3];
sx q[3];
rz(1.928291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55629998) q[2];
sx q[2];
rz(-1.2703398) q[2];
sx q[2];
rz(1.4053819) q[2];
rz(0.65226883) q[3];
sx q[3];
rz(-2.8758719) q[3];
sx q[3];
rz(-2.4238267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2300867) q[0];
sx q[0];
rz(-1.7510887) q[0];
sx q[0];
rz(-0.44816309) q[0];
rz(-2.3975092) q[1];
sx q[1];
rz(-2.4241445) q[1];
sx q[1];
rz(-1.4345899) q[1];
rz(-2.5581735) q[2];
sx q[2];
rz(-1.7023682) q[2];
sx q[2];
rz(-2.792991) q[2];
rz(-0.35318315) q[3];
sx q[3];
rz(-2.1174869) q[3];
sx q[3];
rz(-0.01275851) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
