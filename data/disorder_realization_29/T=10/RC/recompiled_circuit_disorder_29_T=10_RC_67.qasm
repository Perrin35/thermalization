OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1155137) q[0];
sx q[0];
rz(-1.4839412) q[0];
sx q[0];
rz(-0.32615647) q[0];
rz(-1.1905319) q[1];
sx q[1];
rz(-1.3500554) q[1];
sx q[1];
rz(-1.5989369) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21539772) q[0];
sx q[0];
rz(-2.8370259) q[0];
sx q[0];
rz(-2.347441) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7803696) q[2];
sx q[2];
rz(-0.63399705) q[2];
sx q[2];
rz(2.0915973) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1118288) q[1];
sx q[1];
rz(-1.811869) q[1];
sx q[1];
rz(1.3328711) q[1];
rz(1.5845675) q[3];
sx q[3];
rz(-0.78409401) q[3];
sx q[3];
rz(-2.9053094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8618384) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(0.88511434) q[2];
rz(0.72201133) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(-3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9279813) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(2.0425178) q[0];
rz(2.4765769) q[1];
sx q[1];
rz(-1.4140833) q[1];
sx q[1];
rz(0.87759334) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7859902) q[0];
sx q[0];
rz(-2.2539833) q[0];
sx q[0];
rz(1.8049294) q[0];
rz(-pi) q[1];
rz(2.7769233) q[2];
sx q[2];
rz(-1.4450577) q[2];
sx q[2];
rz(2.113935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7922389) q[1];
sx q[1];
rz(-1.6041479) q[1];
sx q[1];
rz(-1.4637714) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33438501) q[3];
sx q[3];
rz(-1.9352203) q[3];
sx q[3];
rz(0.78727608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7426804) q[2];
sx q[2];
rz(-1.8360527) q[2];
sx q[2];
rz(2.0111283) q[2];
rz(1.8418664) q[3];
sx q[3];
rz(-1.2239417) q[3];
sx q[3];
rz(-1.6931504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14625064) q[0];
sx q[0];
rz(-1.8595707) q[0];
sx q[0];
rz(2.8515942) q[0];
rz(-2.4747804) q[1];
sx q[1];
rz(-1.0338444) q[1];
sx q[1];
rz(-3.0677632) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89798112) q[0];
sx q[0];
rz(-1.6323166) q[0];
sx q[0];
rz(0.089540066) q[0];
x q[1];
rz(-1.9129842) q[2];
sx q[2];
rz(-1.295919) q[2];
sx q[2];
rz(0.16155044) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0970043) q[1];
sx q[1];
rz(-1.5148331) q[1];
sx q[1];
rz(-1.7486649) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35145268) q[3];
sx q[3];
rz(-1.6822527) q[3];
sx q[3];
rz(-2.4748638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6510216) q[2];
sx q[2];
rz(-1.1940424) q[2];
sx q[2];
rz(-0.64669615) q[2];
rz(1.1086639) q[3];
sx q[3];
rz(-0.78287786) q[3];
sx q[3];
rz(1.1289319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.17764238) q[0];
sx q[0];
rz(-2.9681866) q[0];
sx q[0];
rz(1.9529163) q[0];
rz(1.0186609) q[1];
sx q[1];
rz(-0.97266346) q[1];
sx q[1];
rz(1.4368988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1543717) q[0];
sx q[0];
rz(-0.73723307) q[0];
sx q[0];
rz(0.1371951) q[0];
rz(-pi) q[1];
rz(1.2595348) q[2];
sx q[2];
rz(-1.2584104) q[2];
sx q[2];
rz(-0.54158467) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.77364591) q[1];
sx q[1];
rz(-1.2229496) q[1];
sx q[1];
rz(-2.0639973) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7014245) q[3];
sx q[3];
rz(-1.4454953) q[3];
sx q[3];
rz(3.0076722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.556095) q[2];
sx q[2];
rz(-1.4336339) q[2];
sx q[2];
rz(-2.0193224) q[2];
rz(2.1155817) q[3];
sx q[3];
rz(-2.3882073) q[3];
sx q[3];
rz(-2.1508353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.4595903) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(-2.7686152) q[0];
rz(-2.9176118) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(1.3164828) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7861159) q[0];
sx q[0];
rz(-1.1899723) q[0];
sx q[0];
rz(-2.064408) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3773642) q[2];
sx q[2];
rz(-0.35596213) q[2];
sx q[2];
rz(2.228235) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.74354467) q[1];
sx q[1];
rz(-3.0882356) q[1];
sx q[1];
rz(-1.3392901) q[1];
rz(0.54721197) q[3];
sx q[3];
rz(-1.1037165) q[3];
sx q[3];
rz(-2.3576749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.10107723) q[2];
sx q[2];
rz(-0.83156362) q[2];
sx q[2];
rz(-2.3357847) q[2];
rz(0.53330437) q[3];
sx q[3];
rz(-2.0093982) q[3];
sx q[3];
rz(1.0114975) q[3];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39078113) q[0];
sx q[0];
rz(-1.3178786) q[0];
sx q[0];
rz(0.090963013) q[0];
rz(-2.2816351) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(-1.3202753) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66619191) q[0];
sx q[0];
rz(-1.1328508) q[0];
sx q[0];
rz(1.9802718) q[0];
rz(1.7485511) q[2];
sx q[2];
rz(-2.1338935) q[2];
sx q[2];
rz(-0.69001889) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7399866) q[1];
sx q[1];
rz(-2.0410129) q[1];
sx q[1];
rz(-0.92958881) q[1];
rz(-pi) q[2];
rz(-2.8517013) q[3];
sx q[3];
rz(-2.206466) q[3];
sx q[3];
rz(-2.1260726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0992574) q[2];
sx q[2];
rz(-0.96747413) q[2];
sx q[2];
rz(-2.5406204) q[2];
rz(-2.6565334) q[3];
sx q[3];
rz(-2.9197013) q[3];
sx q[3];
rz(-1.4453567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27959529) q[0];
sx q[0];
rz(-1.1789362) q[0];
sx q[0];
rz(2.5860508) q[0];
rz(-0.034596054) q[1];
sx q[1];
rz(-2.3831773) q[1];
sx q[1];
rz(-1.3909891) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3334675) q[0];
sx q[0];
rz(-1.2398749) q[0];
sx q[0];
rz(0.56751859) q[0];
rz(-pi) q[1];
rz(-0.58629845) q[2];
sx q[2];
rz(-0.43160298) q[2];
sx q[2];
rz(-1.4177711) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2512867) q[1];
sx q[1];
rz(-1.444) q[1];
sx q[1];
rz(2.2433953) q[1];
x q[2];
rz(0.075918003) q[3];
sx q[3];
rz(-0.77709353) q[3];
sx q[3];
rz(2.0457552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3283219) q[2];
sx q[2];
rz(-2.6997824) q[2];
sx q[2];
rz(2.7461046) q[2];
rz(-1.288712) q[3];
sx q[3];
rz(-1.5356531) q[3];
sx q[3];
rz(-2.4718463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.678858) q[0];
sx q[0];
rz(-0.33927074) q[0];
sx q[0];
rz(-1.4920374) q[0];
rz(-2.18816) q[1];
sx q[1];
rz(-2.0326734) q[1];
sx q[1];
rz(1.4377726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4576365) q[0];
sx q[0];
rz(-0.55011612) q[0];
sx q[0];
rz(-1.9253299) q[0];
rz(-pi) q[1];
rz(1.4346801) q[2];
sx q[2];
rz(-1.4021177) q[2];
sx q[2];
rz(-2.8240734) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.86086035) q[1];
sx q[1];
rz(-2.3128465) q[1];
sx q[1];
rz(-0.6764722) q[1];
rz(0.020047232) q[3];
sx q[3];
rz(-2.0181977) q[3];
sx q[3];
rz(0.1964387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4961204) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(0.17871857) q[2];
rz(-0.86137613) q[3];
sx q[3];
rz(-1.9390315) q[3];
sx q[3];
rz(-2.7698959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9361967) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(2.0478915) q[0];
rz(0.73668346) q[1];
sx q[1];
rz(-1.8700347) q[1];
sx q[1];
rz(1.0587943) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7208784) q[0];
sx q[0];
rz(-1.2982561) q[0];
sx q[0];
rz(-1.0084624) q[0];
rz(-2.7339897) q[2];
sx q[2];
rz(-2.0098445) q[2];
sx q[2];
rz(1.1328732) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4150548) q[1];
sx q[1];
rz(-1.6931567) q[1];
sx q[1];
rz(1.9499669) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9731673) q[3];
sx q[3];
rz(-1.1306922) q[3];
sx q[3];
rz(2.6228867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8686691) q[2];
sx q[2];
rz(-0.69512853) q[2];
sx q[2];
rz(-1.6607364) q[2];
rz(2.7311834) q[3];
sx q[3];
rz(-1.7216262) q[3];
sx q[3];
rz(-0.15795344) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8695628) q[0];
sx q[0];
rz(-0.27150387) q[0];
sx q[0];
rz(-2.8503382) q[0];
rz(0.60925305) q[1];
sx q[1];
rz(-1.6758502) q[1];
sx q[1];
rz(1.7094918) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8490484) q[0];
sx q[0];
rz(-1.2546854) q[0];
sx q[0];
rz(-0.51878099) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3347689) q[2];
sx q[2];
rz(-2.0152976) q[2];
sx q[2];
rz(2.1399463) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6546302) q[1];
sx q[1];
rz(-1.6143867) q[1];
sx q[1];
rz(-0.55200465) q[1];
rz(-pi) q[2];
rz(1.3833952) q[3];
sx q[3];
rz(-2.7029576) q[3];
sx q[3];
rz(-2.8965829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46135205) q[2];
sx q[2];
rz(-2.0246918) q[2];
sx q[2];
rz(1.1414026) q[2];
rz(1.5348148) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(-0.46943584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5466945) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(2.7728511) q[1];
sx q[1];
rz(-1.2422961) q[1];
sx q[1];
rz(3.0098343) q[1];
rz(2.914631) q[2];
sx q[2];
rz(-1.3703176) q[2];
sx q[2];
rz(-0.3581518) q[2];
rz(-0.42320078) q[3];
sx q[3];
rz(-1.3907708) q[3];
sx q[3];
rz(1.6650865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];