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
rz(-2.6540304) q[0];
sx q[0];
rz(-2.7609479) q[0];
sx q[0];
rz(-0.08925499) q[0];
rz(2.6788977) q[1];
sx q[1];
rz(-0.2806288) q[1];
sx q[1];
rz(-1.7791003) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4498964) q[0];
sx q[0];
rz(-1.5498156) q[0];
sx q[0];
rz(-2.7684661) q[0];
rz(-pi) q[1];
rz(-0.14855317) q[2];
sx q[2];
rz(-1.6069716) q[2];
sx q[2];
rz(0.52001563) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.465724) q[1];
sx q[1];
rz(-1.3480519) q[1];
sx q[1];
rz(-0.59091076) q[1];
rz(-pi) q[2];
rz(-2.9611601) q[3];
sx q[3];
rz(-0.46262925) q[3];
sx q[3];
rz(-0.53354665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6224299) q[2];
sx q[2];
rz(-2.3795655) q[2];
sx q[2];
rz(-2.2117174) q[2];
rz(-2.8376288) q[3];
sx q[3];
rz(-1.3320987) q[3];
sx q[3];
rz(-0.16873321) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8836483) q[0];
sx q[0];
rz(-1.5771663) q[0];
sx q[0];
rz(-1.2516578) q[0];
rz(-2.9257863) q[1];
sx q[1];
rz(-2.1025175) q[1];
sx q[1];
rz(3.0024517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3045791) q[0];
sx q[0];
rz(-2.4615335) q[0];
sx q[0];
rz(1.3152298) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6433442) q[2];
sx q[2];
rz(-1.8115269) q[2];
sx q[2];
rz(-1.6290851) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8976921) q[1];
sx q[1];
rz(-0.75112126) q[1];
sx q[1];
rz(0.35219197) q[1];
x q[2];
rz(2.4226133) q[3];
sx q[3];
rz(-2.1785979) q[3];
sx q[3];
rz(-1.833235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6824049) q[2];
sx q[2];
rz(-0.46899691) q[2];
sx q[2];
rz(2.0501308) q[2];
rz(-1.546321) q[3];
sx q[3];
rz(-1.073444) q[3];
sx q[3];
rz(2.6923164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0794373) q[0];
sx q[0];
rz(-1.2261483) q[0];
sx q[0];
rz(-3.1297744) q[0];
rz(-2.2305409) q[1];
sx q[1];
rz(-1.751519) q[1];
sx q[1];
rz(-1.3284838) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26052159) q[0];
sx q[0];
rz(-2.0876472) q[0];
sx q[0];
rz(-1.3840186) q[0];
x q[1];
rz(1.4558037) q[2];
sx q[2];
rz(-1.7604239) q[2];
sx q[2];
rz(-2.6791089) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40483002) q[1];
sx q[1];
rz(-1.1859845) q[1];
sx q[1];
rz(2.3865108) q[1];
rz(-1.4592917) q[3];
sx q[3];
rz(-1.6435677) q[3];
sx q[3];
rz(-2.7804216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38982424) q[2];
sx q[2];
rz(-1.6168892) q[2];
sx q[2];
rz(0.41979182) q[2];
rz(-1.7828364) q[3];
sx q[3];
rz(-2.8163781) q[3];
sx q[3];
rz(-1.9512008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6841986) q[0];
sx q[0];
rz(-1.1379108) q[0];
sx q[0];
rz(-0.5140636) q[0];
rz(-2.7534516) q[1];
sx q[1];
rz(-0.61182794) q[1];
sx q[1];
rz(1.2302037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83004763) q[0];
sx q[0];
rz(-1.8732949) q[0];
sx q[0];
rz(-1.8414458) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4172961) q[2];
sx q[2];
rz(-1.2822552) q[2];
sx q[2];
rz(-2.9674825) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9446805) q[1];
sx q[1];
rz(-0.80183235) q[1];
sx q[1];
rz(2.6242328) q[1];
x q[2];
rz(1.3421273) q[3];
sx q[3];
rz(-0.91311753) q[3];
sx q[3];
rz(2.2715037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.5019044) q[2];
sx q[2];
rz(-1.014726) q[2];
sx q[2];
rz(-0.91628966) q[2];
rz(0.4942975) q[3];
sx q[3];
rz(-1.1980134) q[3];
sx q[3];
rz(2.3673207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2275527) q[0];
sx q[0];
rz(-2.1844449) q[0];
sx q[0];
rz(-0.48602948) q[0];
rz(-2.5388429) q[1];
sx q[1];
rz(-1.175368) q[1];
sx q[1];
rz(-2.026162) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.225352) q[0];
sx q[0];
rz(-1.1778206) q[0];
sx q[0];
rz(-0.66182242) q[0];
x q[1];
rz(0.79430842) q[2];
sx q[2];
rz(-1.2758848) q[2];
sx q[2];
rz(-2.8821534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.045687606) q[1];
sx q[1];
rz(-2.3153911) q[1];
sx q[1];
rz(-0.5788486) q[1];
rz(-pi) q[2];
rz(2.0511469) q[3];
sx q[3];
rz(-1.3355888) q[3];
sx q[3];
rz(0.0071098162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1439765) q[2];
sx q[2];
rz(-2.82085) q[2];
sx q[2];
rz(2.0799267) q[2];
rz(1.6210506) q[3];
sx q[3];
rz(-1.942626) q[3];
sx q[3];
rz(-0.89402136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0735556) q[0];
sx q[0];
rz(-0.03980045) q[0];
sx q[0];
rz(-1.918248) q[0];
rz(-2.6262737) q[1];
sx q[1];
rz(-0.66910187) q[1];
sx q[1];
rz(-1.3656778) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66636234) q[0];
sx q[0];
rz(-1.7710553) q[0];
sx q[0];
rz(-3.0259575) q[0];
x q[1];
rz(-1.7530982) q[2];
sx q[2];
rz(-1.4788091) q[2];
sx q[2];
rz(-2.0192041) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6462145) q[1];
sx q[1];
rz(-1.7217024) q[1];
sx q[1];
rz(-0.47980984) q[1];
rz(-pi) q[2];
rz(2.8441098) q[3];
sx q[3];
rz(-1.5362036) q[3];
sx q[3];
rz(-0.57719798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0627275) q[2];
sx q[2];
rz(-1.3389503) q[2];
sx q[2];
rz(1.8143181) q[2];
rz(1.7257388) q[3];
sx q[3];
rz(-3.0684107) q[3];
sx q[3];
rz(2.2590526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13407229) q[0];
sx q[0];
rz(-0.53711689) q[0];
sx q[0];
rz(-2.9448331) q[0];
rz(2.9255731) q[1];
sx q[1];
rz(-1.0088423) q[1];
sx q[1];
rz(-0.78741995) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8926933) q[0];
sx q[0];
rz(-1.3266047) q[0];
sx q[0];
rz(0.1949744) q[0];
rz(-pi) q[1];
rz(1.2741267) q[2];
sx q[2];
rz(-2.0624522) q[2];
sx q[2];
rz(1.1664388) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43147165) q[1];
sx q[1];
rz(-2.1619051) q[1];
sx q[1];
rz(1.6256871) q[1];
rz(-pi) q[2];
rz(-1.0490956) q[3];
sx q[3];
rz(-0.36715436) q[3];
sx q[3];
rz(0.01361135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54062033) q[2];
sx q[2];
rz(-1.629402) q[2];
sx q[2];
rz(-0.38360325) q[2];
rz(2.2331734) q[3];
sx q[3];
rz(-2.3288265) q[3];
sx q[3];
rz(-0.15373357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8613794) q[0];
sx q[0];
rz(-0.52556831) q[0];
sx q[0];
rz(-0.63661611) q[0];
rz(2.5909297) q[1];
sx q[1];
rz(-1.6267585) q[1];
sx q[1];
rz(-0.47264636) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20339676) q[0];
sx q[0];
rz(-2.2455611) q[0];
sx q[0];
rz(1.4569253) q[0];
rz(0.093229712) q[2];
sx q[2];
rz(-1.5179253) q[2];
sx q[2];
rz(0.44651647) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1324942) q[1];
sx q[1];
rz(-2.9995697) q[1];
sx q[1];
rz(-0.36722398) q[1];
x q[2];
rz(-2.74284) q[3];
sx q[3];
rz(-0.67091432) q[3];
sx q[3];
rz(-0.89445597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7134646) q[2];
sx q[2];
rz(-1.2423923) q[2];
sx q[2];
rz(0.098527519) q[2];
rz(0.030700961) q[3];
sx q[3];
rz(-1.5714329) q[3];
sx q[3];
rz(-2.9876685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.6105662) q[0];
sx q[0];
rz(-2.1838146) q[0];
sx q[0];
rz(-0.73262334) q[0];
rz(-0.0064119617) q[1];
sx q[1];
rz(-2.1156204) q[1];
sx q[1];
rz(1.4220672) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5310282) q[0];
sx q[0];
rz(-2.8262666) q[0];
sx q[0];
rz(2.226693) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2071916) q[2];
sx q[2];
rz(-1.0997314) q[2];
sx q[2];
rz(-2.1944012) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8274424) q[1];
sx q[1];
rz(-1.536751) q[1];
sx q[1];
rz(-0.15383792) q[1];
rz(-pi) q[2];
rz(2.6220065) q[3];
sx q[3];
rz(-1.5056464) q[3];
sx q[3];
rz(-1.6217854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0473359) q[2];
sx q[2];
rz(-1.7862659) q[2];
sx q[2];
rz(2.6195841) q[2];
rz(3.1212854) q[3];
sx q[3];
rz(-0.69965196) q[3];
sx q[3];
rz(-2.844753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.2688399) q[0];
sx q[0];
rz(-1.6167384) q[0];
sx q[0];
rz(-2.010345) q[0];
rz(-2.3616683) q[1];
sx q[1];
rz(-1.9386539) q[1];
sx q[1];
rz(-2.6663229) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4085379) q[0];
sx q[0];
rz(-0.86472874) q[0];
sx q[0];
rz(2.6557794) q[0];
x q[1];
rz(1.8653706) q[2];
sx q[2];
rz(-2.4524052) q[2];
sx q[2];
rz(-2.1222494) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1058029) q[1];
sx q[1];
rz(-0.33638182) q[1];
sx q[1];
rz(-0.85806429) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9679139) q[3];
sx q[3];
rz(-2.4376737) q[3];
sx q[3];
rz(0.46797637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7705226) q[2];
sx q[2];
rz(-1.5757685) q[2];
sx q[2];
rz(1.2062262) q[2];
rz(1.3953588) q[3];
sx q[3];
rz(-0.54308707) q[3];
sx q[3];
rz(-0.079308184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64410011) q[0];
sx q[0];
rz(-2.946749) q[0];
sx q[0];
rz(-0.83644833) q[0];
rz(2.1438228) q[1];
sx q[1];
rz(-2.0496968) q[1];
sx q[1];
rz(-3.0539378) q[1];
rz(0.069167698) q[2];
sx q[2];
rz(-1.7887049) q[2];
sx q[2];
rz(-0.6520581) q[2];
rz(1.6818123) q[3];
sx q[3];
rz(-2.2465977) q[3];
sx q[3];
rz(-2.6096596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
