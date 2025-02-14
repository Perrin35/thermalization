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
rz(-0.4975118) q[0];
sx q[0];
rz(-1.8026135) q[0];
sx q[0];
rz(-2.8259377) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(-0.38771954) q[1];
sx q[1];
rz(0.60400909) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9785288) q[0];
sx q[0];
rz(-1.8730358) q[0];
sx q[0];
rz(-3.0885484) q[0];
x q[1];
rz(-2.6409615) q[2];
sx q[2];
rz(-0.96304578) q[2];
sx q[2];
rz(2.1082102) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65434124) q[1];
sx q[1];
rz(-0.62955925) q[1];
sx q[1];
rz(-0.97626026) q[1];
rz(1.2617926) q[3];
sx q[3];
rz(-1.4674392) q[3];
sx q[3];
rz(0.37783937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9903367) q[2];
sx q[2];
rz(-1.9467111) q[2];
sx q[2];
rz(-1.7830431) q[2];
rz(1.5938866) q[3];
sx q[3];
rz(-2.2113776) q[3];
sx q[3];
rz(1.5096629) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58594054) q[0];
sx q[0];
rz(-2.5037615) q[0];
sx q[0];
rz(-1.9441388) q[0];
rz(2.6400631) q[1];
sx q[1];
rz(-1.5634368) q[1];
sx q[1];
rz(-1.7053568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7345463) q[0];
sx q[0];
rz(-2.1384031) q[0];
sx q[0];
rz(-0.27510402) q[0];
rz(-2.4386232) q[2];
sx q[2];
rz(-1.6466993) q[2];
sx q[2];
rz(3.0413306) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3489248) q[1];
sx q[1];
rz(-1.3745921) q[1];
sx q[1];
rz(1.6755107) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40879876) q[3];
sx q[3];
rz(-2.650947) q[3];
sx q[3];
rz(-1.7237494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88840914) q[2];
sx q[2];
rz(-1.2359572) q[2];
sx q[2];
rz(0.02296981) q[2];
rz(-0.90211558) q[3];
sx q[3];
rz(-1.918856) q[3];
sx q[3];
rz(0.42660108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4334634) q[0];
sx q[0];
rz(-1.7925649) q[0];
sx q[0];
rz(-0.9915114) q[0];
rz(2.1082711) q[1];
sx q[1];
rz(-2.8585377) q[1];
sx q[1];
rz(-1.2352157) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18444321) q[0];
sx q[0];
rz(-2.1547124) q[0];
sx q[0];
rz(-1.2725485) q[0];
rz(2.5052983) q[2];
sx q[2];
rz(-2.4264196) q[2];
sx q[2];
rz(-1.6805122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8950303) q[1];
sx q[1];
rz(-0.8324648) q[1];
sx q[1];
rz(-0.33729302) q[1];
x q[2];
rz(2.3091812) q[3];
sx q[3];
rz(-0.73736546) q[3];
sx q[3];
rz(-0.6399782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11423763) q[2];
sx q[2];
rz(-2.3894252) q[2];
sx q[2];
rz(0.99119622) q[2];
rz(1.2973971) q[3];
sx q[3];
rz(-1.8810561) q[3];
sx q[3];
rz(1.7245002) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39795136) q[0];
sx q[0];
rz(-0.93248168) q[0];
sx q[0];
rz(-2.3241296) q[0];
rz(-2.815333) q[1];
sx q[1];
rz(-1.9605109) q[1];
sx q[1];
rz(-0.78525966) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2366339) q[0];
sx q[0];
rz(-1.8371305) q[0];
sx q[0];
rz(1.8840709) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26008028) q[2];
sx q[2];
rz(-2.4818015) q[2];
sx q[2];
rz(1.545797) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7689159) q[1];
sx q[1];
rz(-2.0119889) q[1];
sx q[1];
rz(0.4695453) q[1];
x q[2];
rz(0.35492949) q[3];
sx q[3];
rz(-1.4961492) q[3];
sx q[3];
rz(-1.9982893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0205445) q[2];
sx q[2];
rz(-2.497017) q[2];
sx q[2];
rz(-0.56309593) q[2];
rz(-2.2480887) q[3];
sx q[3];
rz(-1.9681135) q[3];
sx q[3];
rz(1.2347429) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9017482) q[0];
sx q[0];
rz(-0.32308602) q[0];
sx q[0];
rz(-1.5455986) q[0];
rz(2.1196938) q[1];
sx q[1];
rz(-1.1232168) q[1];
sx q[1];
rz(-2.9603069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3485506) q[0];
sx q[0];
rz(-0.53213813) q[0];
sx q[0];
rz(1.1842968) q[0];
rz(-pi) q[1];
rz(-0.47563817) q[2];
sx q[2];
rz(-1.5283094) q[2];
sx q[2];
rz(-1.795639) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0376772) q[1];
sx q[1];
rz(-1.5427151) q[1];
sx q[1];
rz(1.7930536) q[1];
rz(1.7277754) q[3];
sx q[3];
rz(-2.008956) q[3];
sx q[3];
rz(-2.639132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.01866092) q[2];
sx q[2];
rz(-0.67492008) q[2];
sx q[2];
rz(1.8956511) q[2];
rz(2.5490226) q[3];
sx q[3];
rz(-1.4453245) q[3];
sx q[3];
rz(-0.84111253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5331921) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(-2.8327508) q[0];
rz(-0.32288512) q[1];
sx q[1];
rz(-2.7643118) q[1];
sx q[1];
rz(3.0016532) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4734264) q[0];
sx q[0];
rz(-0.65957171) q[0];
sx q[0];
rz(-1.1585328) q[0];
rz(-pi) q[1];
rz(-1.1148648) q[2];
sx q[2];
rz(-1.6559634) q[2];
sx q[2];
rz(-0.87907253) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9819239) q[1];
sx q[1];
rz(-2.2743011) q[1];
sx q[1];
rz(-2.7419006) q[1];
rz(-pi) q[2];
rz(1.1902963) q[3];
sx q[3];
rz(-0.69529136) q[3];
sx q[3];
rz(2.8075346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66190019) q[2];
sx q[2];
rz(-0.65325824) q[2];
sx q[2];
rz(-0.24197401) q[2];
rz(-0.74609977) q[3];
sx q[3];
rz(-1.7753764) q[3];
sx q[3];
rz(-1.7297176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.022843) q[0];
sx q[0];
rz(-1.043909) q[0];
sx q[0];
rz(1.2710849) q[0];
rz(-2.5785043) q[1];
sx q[1];
rz(-0.92910281) q[1];
sx q[1];
rz(-1.7005327) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25172397) q[0];
sx q[0];
rz(-1.2996724) q[0];
sx q[0];
rz(-0.95377484) q[0];
x q[1];
rz(-2.7602642) q[2];
sx q[2];
rz(-1.6994119) q[2];
sx q[2];
rz(-0.53024697) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5477834) q[1];
sx q[1];
rz(-1.154222) q[1];
sx q[1];
rz(-1.3273456) q[1];
x q[2];
rz(-1.3541) q[3];
sx q[3];
rz(-0.85734493) q[3];
sx q[3];
rz(0.59123033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.21961221) q[2];
sx q[2];
rz(-1.8461123) q[2];
sx q[2];
rz(1.8249576) q[2];
rz(0.5611788) q[3];
sx q[3];
rz(-2.1771274) q[3];
sx q[3];
rz(1.6130028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0372666) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(0.88687801) q[0];
rz(0.2001702) q[1];
sx q[1];
rz(-1.472241) q[1];
sx q[1];
rz(1.0221457) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8093531) q[0];
sx q[0];
rz(-0.13357559) q[0];
sx q[0];
rz(1.3788401) q[0];
rz(-pi) q[1];
rz(2.4046477) q[2];
sx q[2];
rz(-0.96074694) q[2];
sx q[2];
rz(1.7231154) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31412582) q[1];
sx q[1];
rz(-2.3955401) q[1];
sx q[1];
rz(1.5681727) q[1];
rz(-1.2194013) q[3];
sx q[3];
rz(-1.0236275) q[3];
sx q[3];
rz(-1.6148293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.48978051) q[2];
sx q[2];
rz(-0.45280364) q[2];
sx q[2];
rz(2.32302) q[2];
rz(-0.98322785) q[3];
sx q[3];
rz(-2.1849617) q[3];
sx q[3];
rz(2.6035068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035456903) q[0];
sx q[0];
rz(-1.9945194) q[0];
sx q[0];
rz(1.5863093) q[0];
rz(0.69681329) q[1];
sx q[1];
rz(-0.15334829) q[1];
sx q[1];
rz(-0.78479016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1479758) q[0];
sx q[0];
rz(-2.412556) q[0];
sx q[0];
rz(-2.4528798) q[0];
rz(-pi) q[1];
rz(-1.3680787) q[2];
sx q[2];
rz(-1.4215111) q[2];
sx q[2];
rz(-1.9307856) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.49336772) q[1];
sx q[1];
rz(-2.7052212) q[1];
sx q[1];
rz(-3.0619951) q[1];
rz(-0.611245) q[3];
sx q[3];
rz(-1.9211624) q[3];
sx q[3];
rz(0.65920748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67184225) q[2];
sx q[2];
rz(-1.8639114) q[2];
sx q[2];
rz(-2.3308241) q[2];
rz(1.9226711) q[3];
sx q[3];
rz(-2.6007077) q[3];
sx q[3];
rz(1.0651275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3588381) q[0];
sx q[0];
rz(-2.8420119) q[0];
sx q[0];
rz(1.4240356) q[0];
rz(-2.3639823) q[1];
sx q[1];
rz(-0.97336665) q[1];
sx q[1];
rz(-0.75540677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9841254) q[0];
sx q[0];
rz(-1.3608772) q[0];
sx q[0];
rz(-0.0046878417) q[0];
rz(-pi) q[1];
rz(0.39647409) q[2];
sx q[2];
rz(-1.586953) q[2];
sx q[2];
rz(2.0153869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.81138071) q[1];
sx q[1];
rz(-0.63618681) q[1];
sx q[1];
rz(2.5274171) q[1];
rz(-pi) q[2];
rz(-2.741119) q[3];
sx q[3];
rz(-2.4254834) q[3];
sx q[3];
rz(3.0587089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5648254) q[2];
sx q[2];
rz(-0.29594031) q[2];
sx q[2];
rz(2.9887065) q[2];
rz(-1.2711924) q[3];
sx q[3];
rz(-1.252424) q[3];
sx q[3];
rz(-0.66711867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85970238) q[0];
sx q[0];
rz(-2.6489881) q[0];
sx q[0];
rz(-0.76982605) q[0];
rz(1.9181171) q[1];
sx q[1];
rz(-1.2212831) q[1];
sx q[1];
rz(2.4881359) q[1];
rz(0.59719795) q[2];
sx q[2];
rz(-1.9382678) q[2];
sx q[2];
rz(1.1282327) q[2];
rz(1.1673234) q[3];
sx q[3];
rz(-1.2444166) q[3];
sx q[3];
rz(1.9593596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
