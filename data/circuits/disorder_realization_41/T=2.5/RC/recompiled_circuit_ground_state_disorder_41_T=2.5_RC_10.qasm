OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95076686) q[0];
sx q[0];
rz(-1.8523676) q[0];
sx q[0];
rz(-3.0486795) q[0];
rz(-0.095280401) q[1];
sx q[1];
rz(2.4089101) q[1];
sx q[1];
rz(10.664193) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5372972) q[0];
sx q[0];
rz(-1.822246) q[0];
sx q[0];
rz(-0.23668134) q[0];
rz(2.8777166) q[2];
sx q[2];
rz(-1.3812764) q[2];
sx q[2];
rz(-0.60631982) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.69529205) q[1];
sx q[1];
rz(-1.9009034) q[1];
sx q[1];
rz(1.5395201) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69783437) q[3];
sx q[3];
rz(-2.4437332) q[3];
sx q[3];
rz(-2.7891911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4702845) q[2];
sx q[2];
rz(-1.6613864) q[2];
sx q[2];
rz(-1.3489464) q[2];
rz(0.74364439) q[3];
sx q[3];
rz(-2.8642004) q[3];
sx q[3];
rz(-2.9644137) q[3];
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
rz(-1.0002366) q[0];
sx q[0];
rz(-1.5445222) q[0];
sx q[0];
rz(-0.29362383) q[0];
rz(-1.0307505) q[1];
sx q[1];
rz(-1.3615969) q[1];
sx q[1];
rz(-0.34293276) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021504952) q[0];
sx q[0];
rz(-1.3653127) q[0];
sx q[0];
rz(2.1028825) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8937102) q[2];
sx q[2];
rz(-0.39034778) q[2];
sx q[2];
rz(0.54803145) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98951116) q[1];
sx q[1];
rz(-0.73484269) q[1];
sx q[1];
rz(-2.8192725) q[1];
rz(-pi) q[2];
rz(-1.9092122) q[3];
sx q[3];
rz(-2.4840601) q[3];
sx q[3];
rz(-0.38207182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.12884101) q[2];
sx q[2];
rz(-1.7495456) q[2];
sx q[2];
rz(0.7592321) q[2];
rz(0.020708474) q[3];
sx q[3];
rz(-1.2660675) q[3];
sx q[3];
rz(-0.43829632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52221209) q[0];
sx q[0];
rz(-2.5396357) q[0];
sx q[0];
rz(0.0044599175) q[0];
rz(-2.4941173) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(-0.78537816) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68092184) q[0];
sx q[0];
rz(-0.17853949) q[0];
sx q[0];
rz(0.4692678) q[0];
rz(-pi) q[1];
rz(1.0713646) q[2];
sx q[2];
rz(-0.3541666) q[2];
sx q[2];
rz(-3.0234697) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0177893) q[1];
sx q[1];
rz(-2.5964156) q[1];
sx q[1];
rz(1.7032318) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3397869) q[3];
sx q[3];
rz(-0.84322819) q[3];
sx q[3];
rz(-1.3096732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7711827) q[2];
sx q[2];
rz(-0.9202756) q[2];
sx q[2];
rz(-1.7837589) q[2];
rz(0.34058288) q[3];
sx q[3];
rz(-0.95652306) q[3];
sx q[3];
rz(0.79184872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1444645) q[0];
sx q[0];
rz(-2.0576532) q[0];
sx q[0];
rz(0.88515627) q[0];
rz(2.3947233) q[1];
sx q[1];
rz(-0.62364686) q[1];
sx q[1];
rz(-1.0964099) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089450739) q[0];
sx q[0];
rz(-2.6656287) q[0];
sx q[0];
rz(2.3030998) q[0];
x q[1];
rz(-1.5113513) q[2];
sx q[2];
rz(-1.5402093) q[2];
sx q[2];
rz(0.93418834) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8453522) q[1];
sx q[1];
rz(-1.3120323) q[1];
sx q[1];
rz(-0.47611632) q[1];
rz(0.55620749) q[3];
sx q[3];
rz(-2.3283151) q[3];
sx q[3];
rz(-0.63268703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5235644) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(-3.1169685) q[2];
rz(-2.5060182) q[3];
sx q[3];
rz(-0.98819757) q[3];
sx q[3];
rz(-2.2583101) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6898952) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(-2.6746993) q[0];
rz(-0.11058552) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(-1.0708403) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2617874) q[0];
sx q[0];
rz(-2.6222485) q[0];
sx q[0];
rz(-2.065361) q[0];
rz(-pi) q[1];
rz(0.51009615) q[2];
sx q[2];
rz(-0.3613216) q[2];
sx q[2];
rz(-2.9053743) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9164898) q[1];
sx q[1];
rz(-1.1943294) q[1];
sx q[1];
rz(2.589499) q[1];
x q[2];
rz(-2.4399906) q[3];
sx q[3];
rz(-2.1013341) q[3];
sx q[3];
rz(0.70487937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.11833) q[2];
sx q[2];
rz(-1.7795965) q[2];
sx q[2];
rz(-2.2892717) q[2];
rz(-3.1039589) q[3];
sx q[3];
rz(-1.2331542) q[3];
sx q[3];
rz(-0.5948624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0090050176) q[0];
sx q[0];
rz(-1.1786893) q[0];
sx q[0];
rz(1.8527385) q[0];
rz(-1.5124849) q[1];
sx q[1];
rz(-1.1237203) q[1];
sx q[1];
rz(1.1423133) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70002551) q[0];
sx q[0];
rz(-1.4730994) q[0];
sx q[0];
rz(-0.35432451) q[0];
rz(-pi) q[1];
rz(-3.0794607) q[2];
sx q[2];
rz(-1.7781864) q[2];
sx q[2];
rz(-2.1211989) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76872494) q[1];
sx q[1];
rz(-1.8187675) q[1];
sx q[1];
rz(2.7889113) q[1];
rz(2.8303538) q[3];
sx q[3];
rz(-1.4025619) q[3];
sx q[3];
rz(-2.5995035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8196572) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(-0.12040559) q[2];
rz(1.985792) q[3];
sx q[3];
rz(-1.5294231) q[3];
sx q[3];
rz(-0.22404484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8026546) q[0];
sx q[0];
rz(-1.8010362) q[0];
sx q[0];
rz(1.4539723) q[0];
rz(-1.5974207) q[1];
sx q[1];
rz(-1.495196) q[1];
sx q[1];
rz(2.6905751) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16674834) q[0];
sx q[0];
rz(-1.6303501) q[0];
sx q[0];
rz(-0.3015428) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1233929) q[2];
sx q[2];
rz(-0.42358735) q[2];
sx q[2];
rz(0.46679631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1731733) q[1];
sx q[1];
rz(-2.0035158) q[1];
sx q[1];
rz(1.4374025) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0093098442) q[3];
sx q[3];
rz(-0.81039372) q[3];
sx q[3];
rz(-1.1610462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94577998) q[2];
sx q[2];
rz(-1.4142298) q[2];
sx q[2];
rz(1.7808524) q[2];
rz(2.1360548) q[3];
sx q[3];
rz(-1.9660299) q[3];
sx q[3];
rz(-1.7520693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1239531) q[0];
sx q[0];
rz(-0.12549505) q[0];
sx q[0];
rz(2.4355167) q[0];
rz(1.5977244) q[1];
sx q[1];
rz(-1.7192625) q[1];
sx q[1];
rz(-2.2854038) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7995351) q[0];
sx q[0];
rz(-1.5781796) q[0];
sx q[0];
rz(1.5902341) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82243408) q[2];
sx q[2];
rz(-1.2456128) q[2];
sx q[2];
rz(-0.48516824) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.44196196) q[1];
sx q[1];
rz(-2.5285427) q[1];
sx q[1];
rz(-0.35787257) q[1];
x q[2];
rz(1.0821277) q[3];
sx q[3];
rz(-1.4525982) q[3];
sx q[3];
rz(1.4335872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84264821) q[2];
sx q[2];
rz(-1.7405258) q[2];
sx q[2];
rz(-2.974158) q[2];
rz(-1.5542479) q[3];
sx q[3];
rz(-2.3685679) q[3];
sx q[3];
rz(2.1412444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3779959) q[0];
sx q[0];
rz(-1.4507699) q[0];
sx q[0];
rz(2.7698621) q[0];
rz(1.4338214) q[1];
sx q[1];
rz(-1.5377518) q[1];
sx q[1];
rz(2.8415714) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8512631) q[0];
sx q[0];
rz(-1.4545049) q[0];
sx q[0];
rz(-0.28684464) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57711592) q[2];
sx q[2];
rz(-1.5367998) q[2];
sx q[2];
rz(-0.94552065) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2310861) q[1];
sx q[1];
rz(-0.66201895) q[1];
sx q[1];
rz(-0.20259095) q[1];
x q[2];
rz(-1.7172708) q[3];
sx q[3];
rz(-2.4333242) q[3];
sx q[3];
rz(-1.1290916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.47664777) q[2];
sx q[2];
rz(-2.6740394) q[2];
sx q[2];
rz(-2.7139968) q[2];
rz(1.556373) q[3];
sx q[3];
rz(-1.1170324) q[3];
sx q[3];
rz(-0.75470406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4204191) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(2.6721201) q[0];
rz(-0.92533127) q[1];
sx q[1];
rz(-2.053849) q[1];
sx q[1];
rz(0.023177711) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7849451) q[0];
sx q[0];
rz(-1.0482422) q[0];
sx q[0];
rz(1.3808668) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9018578) q[2];
sx q[2];
rz(-1.0791856) q[2];
sx q[2];
rz(0.41504809) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1851481) q[1];
sx q[1];
rz(-1.4063764) q[1];
sx q[1];
rz(-1.7471353) q[1];
rz(-pi) q[2];
rz(-0.18927197) q[3];
sx q[3];
rz(-2.2918275) q[3];
sx q[3];
rz(1.4355575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2716486) q[2];
sx q[2];
rz(-1.9359438) q[2];
sx q[2];
rz(0.49087697) q[2];
rz(2.0513746) q[3];
sx q[3];
rz(-2.1722983) q[3];
sx q[3];
rz(0.31392613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28521095) q[0];
sx q[0];
rz(-2.2638392) q[0];
sx q[0];
rz(1.0765156) q[0];
rz(-0.27676997) q[1];
sx q[1];
rz(-0.77817398) q[1];
sx q[1];
rz(-1.4336817) q[1];
rz(1.8770915) q[2];
sx q[2];
rz(-0.90521348) q[2];
sx q[2];
rz(0.95059849) q[2];
rz(3.1076486) q[3];
sx q[3];
rz(-0.53474075) q[3];
sx q[3];
rz(-0.026215601) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
