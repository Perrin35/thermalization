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
rz(1.6744094) q[0];
sx q[0];
rz(-1.5404584) q[0];
sx q[0];
rz(1.2920446) q[0];
rz(-1.0979103) q[1];
sx q[1];
rz(-2.3854889) q[1];
sx q[1];
rz(0.71192157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5789509) q[0];
sx q[0];
rz(-1.6058086) q[0];
sx q[0];
rz(0.14346153) q[0];
x q[1];
rz(2.8991382) q[2];
sx q[2];
rz(-0.77390817) q[2];
sx q[2];
rz(0.36786181) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8375875) q[1];
sx q[1];
rz(-1.7978354) q[1];
sx q[1];
rz(0.63804896) q[1];
rz(-pi) q[2];
rz(-0.51070945) q[3];
sx q[3];
rz(-1.0282955) q[3];
sx q[3];
rz(1.6245019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0277675) q[2];
sx q[2];
rz(-2.1043468) q[2];
sx q[2];
rz(-0.85565957) q[2];
rz(-1.8850108) q[3];
sx q[3];
rz(-0.62658566) q[3];
sx q[3];
rz(-1.1945061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41334316) q[0];
sx q[0];
rz(-0.11390587) q[0];
sx q[0];
rz(0.86694992) q[0];
rz(-1.4504245) q[1];
sx q[1];
rz(-1.1879286) q[1];
sx q[1];
rz(-0.11122045) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1010872) q[0];
sx q[0];
rz(-1.2049647) q[0];
sx q[0];
rz(2.5554197) q[0];
x q[1];
rz(0.84663518) q[2];
sx q[2];
rz(-0.66613448) q[2];
sx q[2];
rz(-0.29633477) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31383816) q[1];
sx q[1];
rz(-2.1849473) q[1];
sx q[1];
rz(0.78734963) q[1];
x q[2];
rz(-2.2777249) q[3];
sx q[3];
rz(-1.4480531) q[3];
sx q[3];
rz(-1.8254759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1078681) q[2];
sx q[2];
rz(-0.76216951) q[2];
sx q[2];
rz(0.49125683) q[2];
rz(-1.2563502) q[3];
sx q[3];
rz(-1.7137073) q[3];
sx q[3];
rz(2.6826732) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1678109) q[0];
sx q[0];
rz(-1.351492) q[0];
sx q[0];
rz(0.0034045086) q[0];
rz(2.7963474) q[1];
sx q[1];
rz(-0.41989851) q[1];
sx q[1];
rz(-0.87510625) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9291139) q[0];
sx q[0];
rz(-1.0153595) q[0];
sx q[0];
rz(-1.655913) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3904949) q[2];
sx q[2];
rz(-1.688684) q[2];
sx q[2];
rz(2.8468612) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5179123) q[1];
sx q[1];
rz(-2.0225683) q[1];
sx q[1];
rz(2.6623515) q[1];
rz(-pi) q[2];
rz(-0.77121021) q[3];
sx q[3];
rz(-1.589984) q[3];
sx q[3];
rz(-2.8018708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.379999) q[2];
sx q[2];
rz(-0.21445175) q[2];
sx q[2];
rz(0.3271412) q[2];
rz(1.7548148) q[3];
sx q[3];
rz(-1.944647) q[3];
sx q[3];
rz(-0.073971204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4057994) q[0];
sx q[0];
rz(-1.0013591) q[0];
sx q[0];
rz(-1.5843947) q[0];
rz(2.7994075) q[1];
sx q[1];
rz(-2.1125968) q[1];
sx q[1];
rz(2.7142966) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7988069) q[0];
sx q[0];
rz(-1.6248934) q[0];
sx q[0];
rz(1.6038206) q[0];
rz(-0.37125606) q[2];
sx q[2];
rz(-1.0635738) q[2];
sx q[2];
rz(1.4812507) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4327185) q[1];
sx q[1];
rz(-2.1418835) q[1];
sx q[1];
rz(-2.4576748) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56663188) q[3];
sx q[3];
rz(-1.315009) q[3];
sx q[3];
rz(1.8961638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7067318) q[2];
sx q[2];
rz(-0.92919436) q[2];
sx q[2];
rz(-2.2959607) q[2];
rz(0.027806824) q[3];
sx q[3];
rz(-1.3549201) q[3];
sx q[3];
rz(1.4771247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61331767) q[0];
sx q[0];
rz(-1.0266101) q[0];
sx q[0];
rz(-0.051830526) q[0];
rz(1.9822281) q[1];
sx q[1];
rz(-2.4720188) q[1];
sx q[1];
rz(-2.5313098) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7399707) q[0];
sx q[0];
rz(-2.8729872) q[0];
sx q[0];
rz(0.69885077) q[0];
x q[1];
rz(2.8437623) q[2];
sx q[2];
rz(-1.1874842) q[2];
sx q[2];
rz(2.226086) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.38182872) q[1];
sx q[1];
rz(-2.0501158) q[1];
sx q[1];
rz(-1.7698395) q[1];
rz(0.60812833) q[3];
sx q[3];
rz(-1.593489) q[3];
sx q[3];
rz(2.3330027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74756885) q[2];
sx q[2];
rz(-1.3538066) q[2];
sx q[2];
rz(-2.5547011) q[2];
rz(1.759257) q[3];
sx q[3];
rz(-1.047225) q[3];
sx q[3];
rz(-0.97572485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1047644) q[0];
sx q[0];
rz(-2.6703175) q[0];
sx q[0];
rz(-0.18071827) q[0];
rz(1.3982999) q[1];
sx q[1];
rz(-1.2340052) q[1];
sx q[1];
rz(1.3999375) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4961131) q[0];
sx q[0];
rz(-2.3262657) q[0];
sx q[0];
rz(-2.8337571) q[0];
rz(-2.3467881) q[2];
sx q[2];
rz(-0.53404885) q[2];
sx q[2];
rz(0.42511031) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29076642) q[1];
sx q[1];
rz(-1.6572448) q[1];
sx q[1];
rz(-2.269901) q[1];
x q[2];
rz(-0.67541452) q[3];
sx q[3];
rz(-2.3590909) q[3];
sx q[3];
rz(-0.43849643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3900628) q[2];
sx q[2];
rz(-2.2921102) q[2];
sx q[2];
rz(-1.0706104) q[2];
rz(-1.6603598) q[3];
sx q[3];
rz(-2.345572) q[3];
sx q[3];
rz(-1.7695561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(2.2299131) q[0];
sx q[0];
rz(-0.80753082) q[0];
sx q[0];
rz(0.56021488) q[0];
rz(-2.920976) q[1];
sx q[1];
rz(-2.1897774) q[1];
sx q[1];
rz(2.2645948) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79898873) q[0];
sx q[0];
rz(-1.5681802) q[0];
sx q[0];
rz(-1.5883587) q[0];
rz(0.39773293) q[2];
sx q[2];
rz(-1.3215725) q[2];
sx q[2];
rz(-0.13532369) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8940398) q[1];
sx q[1];
rz(-2.8553989) q[1];
sx q[1];
rz(-1.2003837) q[1];
rz(-pi) q[2];
rz(-1.4424004) q[3];
sx q[3];
rz(-0.69046016) q[3];
sx q[3];
rz(2.4336876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.666854) q[2];
sx q[2];
rz(-0.89484221) q[2];
sx q[2];
rz(-1.8043) q[2];
rz(-0.21833359) q[3];
sx q[3];
rz(-1.0206914) q[3];
sx q[3];
rz(-1.4214424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6364994) q[0];
sx q[0];
rz(-1.2525083) q[0];
sx q[0];
rz(-3.002758) q[0];
rz(0.88169634) q[1];
sx q[1];
rz(-1.9813184) q[1];
sx q[1];
rz(0.82463157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0966245) q[0];
sx q[0];
rz(-1.5486778) q[0];
sx q[0];
rz(1.478594) q[0];
x q[1];
rz(0.82197676) q[2];
sx q[2];
rz(-2.9993189) q[2];
sx q[2];
rz(1.7956487) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4887374) q[1];
sx q[1];
rz(-2.0664381) q[1];
sx q[1];
rz(-1.7486217) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10004754) q[3];
sx q[3];
rz(-2.2216019) q[3];
sx q[3];
rz(-1.5318961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7925988) q[2];
sx q[2];
rz(-0.94453347) q[2];
sx q[2];
rz(-2.3080431) q[2];
rz(-0.27093568) q[3];
sx q[3];
rz(-1.1214333) q[3];
sx q[3];
rz(0.84735316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9529097) q[0];
sx q[0];
rz(-1.2485349) q[0];
sx q[0];
rz(0.46440014) q[0];
rz(0.69350997) q[1];
sx q[1];
rz(-0.92718569) q[1];
sx q[1];
rz(0.69673353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7719846) q[0];
sx q[0];
rz(-1.4596997) q[0];
sx q[0];
rz(-0.055776061) q[0];
x q[1];
rz(1.3905754) q[2];
sx q[2];
rz(-1.6107133) q[2];
sx q[2];
rz(0.60933622) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2285959) q[1];
sx q[1];
rz(-1.3218398) q[1];
sx q[1];
rz(0.028214022) q[1];
rz(-pi) q[2];
rz(-0.64033033) q[3];
sx q[3];
rz(-1.5597831) q[3];
sx q[3];
rz(-2.2785694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7742179) q[2];
sx q[2];
rz(-3.0493224) q[2];
sx q[2];
rz(-0.83887678) q[2];
rz(-1.4539666) q[3];
sx q[3];
rz(-1.5435012) q[3];
sx q[3];
rz(-1.669223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.060318) q[0];
sx q[0];
rz(-1.3723137) q[0];
sx q[0];
rz(-1.6171612) q[0];
rz(0.35628191) q[1];
sx q[1];
rz(-1.4189015) q[1];
sx q[1];
rz(-1.3391395) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55521783) q[0];
sx q[0];
rz(-1.5065743) q[0];
sx q[0];
rz(-2.1703266) q[0];
rz(2.9288267) q[2];
sx q[2];
rz(-2.8147455) q[2];
sx q[2];
rz(-2.4642627) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0799651) q[1];
sx q[1];
rz(-1.4486169) q[1];
sx q[1];
rz(0.6811209) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4285092) q[3];
sx q[3];
rz(-1.5350966) q[3];
sx q[3];
rz(-1.3653099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7431006) q[2];
sx q[2];
rz(-2.3091381) q[2];
sx q[2];
rz(-1.8771646) q[2];
rz(-2.8049331) q[3];
sx q[3];
rz(-0.94069898) q[3];
sx q[3];
rz(-1.3338859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3867415) q[0];
sx q[0];
rz(-1.516153) q[0];
sx q[0];
rz(-1.0990912) q[0];
rz(-2.5418538) q[1];
sx q[1];
rz(-0.44282423) q[1];
sx q[1];
rz(2.5437358) q[1];
rz(-1.1588396) q[2];
sx q[2];
rz(-1.2745274) q[2];
sx q[2];
rz(1.6714255) q[2];
rz(1.5532495) q[3];
sx q[3];
rz(-2.1033136) q[3];
sx q[3];
rz(3.1284214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
