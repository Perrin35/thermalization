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
rz(1.9361629) q[0];
sx q[0];
rz(3.087145) q[0];
sx q[0];
rz(8.5589391) q[0];
rz(1.6021597) q[1];
sx q[1];
rz(-1.8973693) q[1];
sx q[1];
rz(3.08334) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6893495) q[0];
sx q[0];
rz(-1.6177243) q[0];
sx q[0];
rz(1.7353135) q[0];
rz(-pi) q[1];
rz(3.0676663) q[2];
sx q[2];
rz(-2.0652899) q[2];
sx q[2];
rz(1.6520776) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73535753) q[1];
sx q[1];
rz(-2.3182097) q[1];
sx q[1];
rz(0.75852345) q[1];
rz(-pi) q[2];
rz(-0.32367651) q[3];
sx q[3];
rz(-1.2403245) q[3];
sx q[3];
rz(-3.0728473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3129348) q[2];
sx q[2];
rz(-0.88064319) q[2];
sx q[2];
rz(-2.177296) q[2];
rz(-0.079484552) q[3];
sx q[3];
rz(-1.8389523) q[3];
sx q[3];
rz(0.77118072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-3.1281779) q[0];
sx q[0];
rz(-2.7981813) q[0];
sx q[0];
rz(-1.5403904) q[0];
rz(-0.84114289) q[1];
sx q[1];
rz(-1.7336188) q[1];
sx q[1];
rz(0.85743633) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090912772) q[0];
sx q[0];
rz(-1.5653531) q[0];
sx q[0];
rz(-1.132227) q[0];
x q[1];
rz(2.2067864) q[2];
sx q[2];
rz(-1.1804712) q[2];
sx q[2];
rz(-1.9644511) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.56164484) q[1];
sx q[1];
rz(-1.0542634) q[1];
sx q[1];
rz(-1.0622611) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4748092) q[3];
sx q[3];
rz(-1.848683) q[3];
sx q[3];
rz(0.052770719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.88605827) q[2];
sx q[2];
rz(-2.8958791) q[2];
sx q[2];
rz(-1.9319755) q[2];
rz(1.3813193) q[3];
sx q[3];
rz(-1.2608903) q[3];
sx q[3];
rz(2.5513726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7968314) q[0];
sx q[0];
rz(-1.3153356) q[0];
sx q[0];
rz(-1.3321846) q[0];
rz(1.6650797) q[1];
sx q[1];
rz(-1.8107199) q[1];
sx q[1];
rz(2.5217893) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7133742) q[0];
sx q[0];
rz(-3.0626903) q[0];
sx q[0];
rz(-2.9419241) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16850833) q[2];
sx q[2];
rz(-1.3187172) q[2];
sx q[2];
rz(-1.3244946) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.42592749) q[1];
sx q[1];
rz(-1.8879379) q[1];
sx q[1];
rz(1.6460544) q[1];
rz(0.26472802) q[3];
sx q[3];
rz(-2.1937498) q[3];
sx q[3];
rz(0.61443751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0009813112) q[2];
sx q[2];
rz(-2.8162214) q[2];
sx q[2];
rz(-0.78835431) q[2];
rz(1.2761448) q[3];
sx q[3];
rz(-1.2272464) q[3];
sx q[3];
rz(2.2787794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9901504) q[0];
sx q[0];
rz(-1.3862415) q[0];
sx q[0];
rz(-1.066712) q[0];
rz(-1.7881296) q[1];
sx q[1];
rz(-0.86694327) q[1];
sx q[1];
rz(-1.9852759) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10945548) q[0];
sx q[0];
rz(-0.14277139) q[0];
sx q[0];
rz(1.7032864) q[0];
x q[1];
rz(0.026211003) q[2];
sx q[2];
rz(-1.6115341) q[2];
sx q[2];
rz(2.7762976) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9383685) q[1];
sx q[1];
rz(-0.48048204) q[1];
sx q[1];
rz(-2.234747) q[1];
x q[2];
rz(-0.66439028) q[3];
sx q[3];
rz(-1.6229462) q[3];
sx q[3];
rz(-0.96168226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.452423) q[2];
sx q[2];
rz(-0.66340476) q[2];
sx q[2];
rz(3.0827674) q[2];
rz(-2.5022653) q[3];
sx q[3];
rz(-1.9202193) q[3];
sx q[3];
rz(-0.85734573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0462129) q[0];
sx q[0];
rz(-1.7261427) q[0];
sx q[0];
rz(0.010490622) q[0];
rz(1.563021) q[1];
sx q[1];
rz(-0.69067162) q[1];
sx q[1];
rz(2.6522327) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2374769) q[0];
sx q[0];
rz(-1.4129708) q[0];
sx q[0];
rz(0.62741168) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0919369) q[2];
sx q[2];
rz(-0.74981028) q[2];
sx q[2];
rz(2.9887226) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.64397821) q[1];
sx q[1];
rz(-2.3757739) q[1];
sx q[1];
rz(0.52180565) q[1];
rz(3.1108256) q[3];
sx q[3];
rz(-1.5602141) q[3];
sx q[3];
rz(-2.2969674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8420777) q[2];
sx q[2];
rz(-1.0966417) q[2];
sx q[2];
rz(3.1411324) q[2];
rz(-0.67356235) q[3];
sx q[3];
rz(-2.2431777) q[3];
sx q[3];
rz(0.7091929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29086581) q[0];
sx q[0];
rz(-1.8274266) q[0];
sx q[0];
rz(-1.8379743) q[0];
rz(-2.8438026) q[1];
sx q[1];
rz(-2.2460263) q[1];
sx q[1];
rz(-2.8477125) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93730629) q[0];
sx q[0];
rz(-1.0244475) q[0];
sx q[0];
rz(2.2541304) q[0];
rz(-pi) q[1];
rz(2.2511036) q[2];
sx q[2];
rz(-1.9142303) q[2];
sx q[2];
rz(2.3517017) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.77569796) q[1];
sx q[1];
rz(-2.0112717) q[1];
sx q[1];
rz(-1.0431402) q[1];
rz(1.3354662) q[3];
sx q[3];
rz(-1.2633861) q[3];
sx q[3];
rz(-3.0855872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46875724) q[2];
sx q[2];
rz(-1.7054649) q[2];
sx q[2];
rz(0.22077665) q[2];
rz(1.2729493) q[3];
sx q[3];
rz(-1.3947398) q[3];
sx q[3];
rz(-1.9934191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4918168) q[0];
sx q[0];
rz(-0.62115541) q[0];
sx q[0];
rz(2.5671) q[0];
rz(0.54221398) q[1];
sx q[1];
rz(-1.503399) q[1];
sx q[1];
rz(2.7508459) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6764561) q[0];
sx q[0];
rz(-2.1628597) q[0];
sx q[0];
rz(-1.3964064) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0894012) q[2];
sx q[2];
rz(-1.0269916) q[2];
sx q[2];
rz(-0.79629818) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4225715) q[1];
sx q[1];
rz(-0.41677058) q[1];
sx q[1];
rz(3.1279492) q[1];
rz(-pi) q[2];
x q[2];
rz(1.607702) q[3];
sx q[3];
rz(-1.2052943) q[3];
sx q[3];
rz(-0.17084641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6812402) q[2];
sx q[2];
rz(-1.1972903) q[2];
sx q[2];
rz(2.528842) q[2];
rz(-0.98226205) q[3];
sx q[3];
rz(-1.8270315) q[3];
sx q[3];
rz(2.6764892) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14923444) q[0];
sx q[0];
rz(-1.5667916) q[0];
sx q[0];
rz(-3.0862578) q[0];
rz(-1.5570359) q[1];
sx q[1];
rz(-2.0535856) q[1];
sx q[1];
rz(2.7016644) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81145168) q[0];
sx q[0];
rz(-2.6127671) q[0];
sx q[0];
rz(2.6374966) q[0];
x q[1];
rz(0.37005927) q[2];
sx q[2];
rz(-0.52071134) q[2];
sx q[2];
rz(-2.0441908) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7288558) q[1];
sx q[1];
rz(-1.2555033) q[1];
sx q[1];
rz(-0.58087491) q[1];
x q[2];
rz(-1.004435) q[3];
sx q[3];
rz(-1.3310588) q[3];
sx q[3];
rz(-2.9842003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.7679193) q[2];
sx q[2];
rz(-1.7934711) q[2];
sx q[2];
rz(-2.8295753) q[2];
rz(0.9196552) q[3];
sx q[3];
rz(-1.7731526) q[3];
sx q[3];
rz(2.9254204) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8298518) q[0];
sx q[0];
rz(-2.1599202) q[0];
sx q[0];
rz(3.1134636) q[0];
rz(1.4460538) q[1];
sx q[1];
rz(-1.1808993) q[1];
sx q[1];
rz(-0.45949724) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1618237) q[0];
sx q[0];
rz(-1.5422675) q[0];
sx q[0];
rz(1.606444) q[0];
rz(-pi) q[1];
rz(-1.2120281) q[2];
sx q[2];
rz(-0.78132403) q[2];
sx q[2];
rz(1.6143022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0952009) q[1];
sx q[1];
rz(-1.4262137) q[1];
sx q[1];
rz(-0.88674366) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0023381) q[3];
sx q[3];
rz(-1.1920658) q[3];
sx q[3];
rz(-1.1564573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10099899) q[2];
sx q[2];
rz(-1.7969635) q[2];
sx q[2];
rz(0.25516587) q[2];
rz(-0.98662871) q[3];
sx q[3];
rz(-0.41472236) q[3];
sx q[3];
rz(-1.7184947) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7607018) q[0];
sx q[0];
rz(-1.0699027) q[0];
sx q[0];
rz(1.799452) q[0];
rz(0.0019207151) q[1];
sx q[1];
rz(-2.0989959) q[1];
sx q[1];
rz(2.0916746) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0107897) q[0];
sx q[0];
rz(-1.783004) q[0];
sx q[0];
rz(-0.15168587) q[0];
rz(-1.5639864) q[2];
sx q[2];
rz(-0.17596888) q[2];
sx q[2];
rz(2.3962452) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.695076) q[1];
sx q[1];
rz(-2.7088266) q[1];
sx q[1];
rz(-2.5997735) q[1];
rz(2.3008651) q[3];
sx q[3];
rz(-1.4089214) q[3];
sx q[3];
rz(-2.9161798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92901015) q[2];
sx q[2];
rz(-1.2139823) q[2];
sx q[2];
rz(-2.6455961) q[2];
rz(1.4604733) q[3];
sx q[3];
rz(-1.3417599) q[3];
sx q[3];
rz(1.1869441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37954189) q[0];
sx q[0];
rz(-2.0370146) q[0];
sx q[0];
rz(1.4203352) q[0];
rz(-0.63915359) q[1];
sx q[1];
rz(-1.2624546) q[1];
sx q[1];
rz(0.75844567) q[1];
rz(0.47472246) q[2];
sx q[2];
rz(-1.793515) q[2];
sx q[2];
rz(-1.1303177) q[2];
rz(-0.81845508) q[3];
sx q[3];
rz(-2.0764805) q[3];
sx q[3];
rz(0.71972328) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
