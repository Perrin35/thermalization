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
rz(0.41843721) q[0];
sx q[0];
rz(5.3199407) q[0];
sx q[0];
rz(9.6286019) q[0];
rz(1.0526429) q[1];
sx q[1];
rz(3.9349603) q[1];
sx q[1];
rz(10.959672) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14102916) q[0];
sx q[0];
rz(-1.805843) q[0];
sx q[0];
rz(1.1147333) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34687545) q[2];
sx q[2];
rz(-0.87793186) q[2];
sx q[2];
rz(-1.3077867) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4272449) q[1];
sx q[1];
rz(-0.74810076) q[1];
sx q[1];
rz(-1.0393049) q[1];
rz(-pi) q[2];
rz(0.20333692) q[3];
sx q[3];
rz(-1.7032188) q[3];
sx q[3];
rz(1.9258969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7242929) q[2];
sx q[2];
rz(-2.5827926) q[2];
sx q[2];
rz(-2.144045) q[2];
rz(2.3857462) q[3];
sx q[3];
rz(-1.3563124) q[3];
sx q[3];
rz(2.1006179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7928829) q[0];
sx q[0];
rz(-3.0443158) q[0];
sx q[0];
rz(-1.8541699) q[0];
rz(0.48311326) q[1];
sx q[1];
rz(-2.099642) q[1];
sx q[1];
rz(-1.9226673) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6605777) q[0];
sx q[0];
rz(-1.7661375) q[0];
sx q[0];
rz(-0.16150772) q[0];
rz(-pi) q[1];
rz(-2.8701209) q[2];
sx q[2];
rz(-2.4515332) q[2];
sx q[2];
rz(-0.59218237) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9049553) q[1];
sx q[1];
rz(-1.7153011) q[1];
sx q[1];
rz(1.852024) q[1];
x q[2];
rz(-1.990149) q[3];
sx q[3];
rz(-2.4046728) q[3];
sx q[3];
rz(1.1518971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7634742) q[2];
sx q[2];
rz(-0.066630445) q[2];
sx q[2];
rz(-0.62180579) q[2];
rz(0.63831896) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(-0.84426713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.3062375) q[0];
sx q[0];
rz(-0.64202809) q[0];
sx q[0];
rz(0.21632347) q[0];
rz(-0.59858876) q[1];
sx q[1];
rz(-1.3052156) q[1];
sx q[1];
rz(-2.6187706) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7824088) q[0];
sx q[0];
rz(-2.6893927) q[0];
sx q[0];
rz(-2.3308999) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7758413) q[2];
sx q[2];
rz(-2.7257082) q[2];
sx q[2];
rz(-1.5378086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.93980689) q[1];
sx q[1];
rz(-0.56957376) q[1];
sx q[1];
rz(-0.42167432) q[1];
x q[2];
rz(-1.9013635) q[3];
sx q[3];
rz(-1.8617612) q[3];
sx q[3];
rz(1.7592848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.181695) q[2];
sx q[2];
rz(-1.8671702) q[2];
sx q[2];
rz(-1.6240906) q[2];
rz(2.6767139) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(-1.6200861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9972123) q[0];
sx q[0];
rz(-0.388044) q[0];
sx q[0];
rz(-1.5390747) q[0];
rz(-2.1082361) q[1];
sx q[1];
rz(-2.3732503) q[1];
sx q[1];
rz(1.8446406) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.813445) q[0];
sx q[0];
rz(-1.7792276) q[0];
sx q[0];
rz(2.3013902) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6455075) q[2];
sx q[2];
rz(-0.69715188) q[2];
sx q[2];
rz(0.32767228) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0721283) q[1];
sx q[1];
rz(-1.1934682) q[1];
sx q[1];
rz(2.0153785) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8861214) q[3];
sx q[3];
rz(-1.0389757) q[3];
sx q[3];
rz(1.1602311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3492655) q[2];
sx q[2];
rz(-2.5040864) q[2];
sx q[2];
rz(2.0959334) q[2];
rz(-2.9271434) q[3];
sx q[3];
rz(-0.72434536) q[3];
sx q[3];
rz(-1.9946056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8056718) q[0];
sx q[0];
rz(-1.2821953) q[0];
sx q[0];
rz(0.51934284) q[0];
rz(-2.1995811) q[1];
sx q[1];
rz(-2.8603034) q[1];
sx q[1];
rz(1.5325783) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62270861) q[0];
sx q[0];
rz(-1.2075338) q[0];
sx q[0];
rz(2.4469482) q[0];
rz(-pi) q[1];
rz(-2.4759859) q[2];
sx q[2];
rz(-0.22826787) q[2];
sx q[2];
rz(-0.98787243) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.457442) q[1];
sx q[1];
rz(-1.5760211) q[1];
sx q[1];
rz(2.0713776) q[1];
rz(-0.10974613) q[3];
sx q[3];
rz(-1.4767495) q[3];
sx q[3];
rz(1.4274244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4970826) q[2];
sx q[2];
rz(-2.3641059) q[2];
sx q[2];
rz(-2.8572148) q[2];
rz(0.55822462) q[3];
sx q[3];
rz(-1.3642718) q[3];
sx q[3];
rz(0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.072642) q[0];
sx q[0];
rz(-0.41256368) q[0];
sx q[0];
rz(-0.64055842) q[0];
rz(-1.6259646) q[1];
sx q[1];
rz(-2.0104355) q[1];
sx q[1];
rz(-0.21387771) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9616674) q[0];
sx q[0];
rz(-2.7583987) q[0];
sx q[0];
rz(-2.2201594) q[0];
x q[1];
rz(-2.7340545) q[2];
sx q[2];
rz(-0.78754163) q[2];
sx q[2];
rz(2.119273) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7678309) q[1];
sx q[1];
rz(-1.6400178) q[1];
sx q[1];
rz(0.94061416) q[1];
x q[2];
rz(-2.0908666) q[3];
sx q[3];
rz(-1.5087391) q[3];
sx q[3];
rz(2.7620535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3512257) q[2];
sx q[2];
rz(-0.48206097) q[2];
sx q[2];
rz(-2.9638929) q[2];
rz(1.7443582) q[3];
sx q[3];
rz(-1.1763562) q[3];
sx q[3];
rz(2.7241404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(3.0419256) q[0];
sx q[0];
rz(-0.61138994) q[0];
sx q[0];
rz(1.1055111) q[0];
rz(-2.8583177) q[1];
sx q[1];
rz(-2.6073644) q[1];
sx q[1];
rz(1.4615321) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7992226) q[0];
sx q[0];
rz(-2.2307751) q[0];
sx q[0];
rz(2.2837385) q[0];
rz(-pi) q[1];
rz(0.93223946) q[2];
sx q[2];
rz(-1.1503714) q[2];
sx q[2];
rz(-1.0300385) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6532306) q[1];
sx q[1];
rz(-1.4547336) q[1];
sx q[1];
rz(-0.56452063) q[1];
rz(-0.57638611) q[3];
sx q[3];
rz(-1.8462876) q[3];
sx q[3];
rz(2.798693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0854411) q[2];
sx q[2];
rz(-1.5957007) q[2];
sx q[2];
rz(-0.20480569) q[2];
rz(-2.0917995) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(2.7753196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7751223) q[0];
sx q[0];
rz(-2.6752052) q[0];
sx q[0];
rz(-3.108016) q[0];
rz(-1.6665943) q[1];
sx q[1];
rz(-2.3971403) q[1];
sx q[1];
rz(-1.9974476) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56267101) q[0];
sx q[0];
rz(-1.5454588) q[0];
sx q[0];
rz(-0.48669215) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2828243) q[2];
sx q[2];
rz(-2.8045553) q[2];
sx q[2];
rz(-2.2971643) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7209789) q[1];
sx q[1];
rz(-0.79809626) q[1];
sx q[1];
rz(-0.53687232) q[1];
rz(-0.28067067) q[3];
sx q[3];
rz(-1.2378527) q[3];
sx q[3];
rz(-1.6368293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7555776) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(-1.2516652) q[2];
rz(1.3517316) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(0.98021093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(0.0018472483) q[0];
sx q[0];
rz(-0.63848764) q[0];
sx q[0];
rz(-1.8852604) q[0];
rz(3.046335) q[1];
sx q[1];
rz(-0.88905159) q[1];
sx q[1];
rz(-2.6108066) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7266709) q[0];
sx q[0];
rz(-1.7309395) q[0];
sx q[0];
rz(-1.8216351) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31802788) q[2];
sx q[2];
rz(-1.5340759) q[2];
sx q[2];
rz(-1.5098234) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.6775178) q[1];
sx q[1];
rz(-1.4428268) q[1];
sx q[1];
rz(-0.24125464) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9050499) q[3];
sx q[3];
rz(-2.6841087) q[3];
sx q[3];
rz(-1.903423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8119767) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(0.54083332) q[2];
rz(0.81816188) q[3];
sx q[3];
rz(-1.6988138) q[3];
sx q[3];
rz(-2.5087859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.624991) q[0];
sx q[0];
rz(-1.7778968) q[0];
sx q[0];
rz(-2.7428108) q[0];
rz(1.7575691) q[1];
sx q[1];
rz(-0.75474352) q[1];
sx q[1];
rz(0.049364518) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.787034) q[0];
sx q[0];
rz(-1.6794717) q[0];
sx q[0];
rz(-3.0804481) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10714679) q[2];
sx q[2];
rz(-1.8320036) q[2];
sx q[2];
rz(0.10814737) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3393664) q[1];
sx q[1];
rz(-0.58103937) q[1];
sx q[1];
rz(-1.4895579) q[1];
rz(-pi) q[2];
rz(-0.63189854) q[3];
sx q[3];
rz(-1.7358586) q[3];
sx q[3];
rz(0.89323211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.12386879) q[2];
sx q[2];
rz(-1.8327291) q[2];
sx q[2];
rz(1.6309942) q[2];
rz(-2.2137568) q[3];
sx q[3];
rz(-1.3414914) q[3];
sx q[3];
rz(-1.9063037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8265726) q[0];
sx q[0];
rz(-0.44354225) q[0];
sx q[0];
rz(-1.0916239) q[0];
rz(-1.0646959) q[1];
sx q[1];
rz(-0.5263435) q[1];
sx q[1];
rz(1.975504) q[1];
rz(2.0666368) q[2];
sx q[2];
rz(-2.0315764) q[2];
sx q[2];
rz(-2.365664) q[2];
rz(2.8228685) q[3];
sx q[3];
rz(-2.095185) q[3];
sx q[3];
rz(-0.17879055) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
