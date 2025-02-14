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
rz(3.8976964) q[1];
sx q[1];
rz(10.1367) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24589892) q[0];
sx q[0];
rz(-0.14764365) q[0];
sx q[0];
rz(2.9013322) q[0];
x q[1];
rz(-0.759077) q[2];
sx q[2];
rz(-1.402194) q[2];
sx q[2];
rz(2.1137266) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3040052) q[1];
sx q[1];
rz(-1.7978354) q[1];
sx q[1];
rz(0.63804896) q[1];
rz(-2.6308832) q[3];
sx q[3];
rz(-1.0282955) q[3];
sx q[3];
rz(1.5170908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11382515) q[2];
sx q[2];
rz(-2.1043468) q[2];
sx q[2];
rz(-2.2859331) q[2];
rz(1.8850108) q[3];
sx q[3];
rz(-2.515007) q[3];
sx q[3];
rz(1.9470866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7282495) q[0];
sx q[0];
rz(-3.0276868) q[0];
sx q[0];
rz(0.86694992) q[0];
rz(-1.4504245) q[1];
sx q[1];
rz(-1.9536641) q[1];
sx q[1];
rz(0.11122045) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0361476) q[0];
sx q[0];
rz(-0.67942028) q[0];
sx q[0];
rz(-2.5359175) q[0];
rz(0.48008474) q[2];
sx q[2];
rz(-1.089555) q[2];
sx q[2];
rz(0.54779886) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7780558) q[1];
sx q[1];
rz(-0.95603525) q[1];
sx q[1];
rz(2.3585783) q[1];
rz(2.2777249) q[3];
sx q[3];
rz(-1.6935395) q[3];
sx q[3];
rz(1.3161167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1078681) q[2];
sx q[2];
rz(-2.3794231) q[2];
sx q[2];
rz(0.49125683) q[2];
rz(1.8852425) q[3];
sx q[3];
rz(-1.4278853) q[3];
sx q[3];
rz(0.45891941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97378174) q[0];
sx q[0];
rz(-1.7901006) q[0];
sx q[0];
rz(-0.0034045086) q[0];
rz(-2.7963474) q[1];
sx q[1];
rz(-2.7216941) q[1];
sx q[1];
rz(2.2664864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7686979) q[0];
sx q[0];
rz(-0.56124306) q[0];
sx q[0];
rz(-0.13613693) q[0];
x q[1];
rz(-3.0217808) q[2];
sx q[2];
rz(-1.3917599) q[2];
sx q[2];
rz(-1.2546292) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5179123) q[1];
sx q[1];
rz(-2.0225683) q[1];
sx q[1];
rz(-2.6623515) q[1];
x q[2];
rz(-3.1140675) q[3];
sx q[3];
rz(-2.3701931) q[3];
sx q[3];
rz(-1.9302544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.7615937) q[2];
sx q[2];
rz(-0.21445175) q[2];
sx q[2];
rz(2.8144515) q[2];
rz(1.3867779) q[3];
sx q[3];
rz(-1.944647) q[3];
sx q[3];
rz(0.073971204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73579329) q[0];
sx q[0];
rz(-1.0013591) q[0];
sx q[0];
rz(1.557198) q[0];
rz(-2.7994075) q[1];
sx q[1];
rz(-1.0289959) q[1];
sx q[1];
rz(-0.4272961) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2262242) q[0];
sx q[0];
rz(-1.5378204) q[0];
sx q[0];
rz(-3.0874662) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37125606) q[2];
sx q[2];
rz(-1.0635738) q[2];
sx q[2];
rz(-1.660342) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8646331) q[1];
sx q[1];
rz(-1.0103419) q[1];
sx q[1];
rz(-2.2629281) q[1];
rz(-1.8713636) q[3];
sx q[3];
rz(-1.0247314) q[3];
sx q[3];
rz(2.6566243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4348609) q[2];
sx q[2];
rz(-2.2123983) q[2];
sx q[2];
rz(2.2959607) q[2];
rz(3.1137858) q[3];
sx q[3];
rz(-1.7866725) q[3];
sx q[3];
rz(-1.6644679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61331767) q[0];
sx q[0];
rz(-1.0266101) q[0];
sx q[0];
rz(-3.0897621) q[0];
rz(-1.9822281) q[1];
sx q[1];
rz(-0.66957384) q[1];
sx q[1];
rz(-2.5313098) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1185111) q[0];
sx q[0];
rz(-1.3661962) q[0];
sx q[0];
rz(-1.746063) q[0];
rz(-0.94172768) q[2];
sx q[2];
rz(-2.6607408) q[2];
sx q[2];
rz(0.22823315) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2817262) q[1];
sx q[1];
rz(-1.3944346) q[1];
sx q[1];
rz(-0.4875261) q[1];
rz(-pi) q[2];
rz(3.1018879) q[3];
sx q[3];
rz(-2.5330946) q[3];
sx q[3];
rz(2.41197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3940238) q[2];
sx q[2];
rz(-1.787786) q[2];
sx q[2];
rz(-0.58689153) q[2];
rz(-1.759257) q[3];
sx q[3];
rz(-2.0943677) q[3];
sx q[3];
rz(2.1658678) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1047644) q[0];
sx q[0];
rz(-0.47127518) q[0];
sx q[0];
rz(0.18071827) q[0];
rz(-1.7432927) q[1];
sx q[1];
rz(-1.2340052) q[1];
sx q[1];
rz(1.3999375) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64547951) q[0];
sx q[0];
rz(-2.3262657) q[0];
sx q[0];
rz(-0.3078356) q[0];
x q[1];
rz(-1.9701875) q[2];
sx q[2];
rz(-1.9353494) q[2];
sx q[2];
rz(-2.6973218) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8508262) q[1];
sx q[1];
rz(-1.4843478) q[1];
sx q[1];
rz(-0.87169164) q[1];
rz(-pi) q[2];
rz(-0.65989699) q[3];
sx q[3];
rz(-1.1142892) q[3];
sx q[3];
rz(2.5259301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3900628) q[2];
sx q[2];
rz(-2.2921102) q[2];
sx q[2];
rz(1.0706104) q[2];
rz(-1.4812329) q[3];
sx q[3];
rz(-0.79602066) q[3];
sx q[3];
rz(1.3720366) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2299131) q[0];
sx q[0];
rz(-2.3340618) q[0];
sx q[0];
rz(0.56021488) q[0];
rz(-2.920976) q[1];
sx q[1];
rz(-0.95181528) q[1];
sx q[1];
rz(0.87699786) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3426039) q[0];
sx q[0];
rz(-1.5734125) q[0];
sx q[0];
rz(-1.5883587) q[0];
x q[1];
rz(-0.58135245) q[2];
sx q[2];
rz(-2.675781) q[2];
sx q[2];
rz(1.9663481) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0045053) q[1];
sx q[1];
rz(-1.3045038) q[1];
sx q[1];
rz(0.10612635) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6991922) q[3];
sx q[3];
rz(-0.69046016) q[3];
sx q[3];
rz(-2.4336876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.47473869) q[2];
sx q[2];
rz(-0.89484221) q[2];
sx q[2];
rz(-1.8043) q[2];
rz(-2.9232591) q[3];
sx q[3];
rz(-1.0206914) q[3];
sx q[3];
rz(-1.7201503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6364994) q[0];
sx q[0];
rz(-1.2525083) q[0];
sx q[0];
rz(3.002758) q[0];
rz(2.2598963) q[1];
sx q[1];
rz(-1.9813184) q[1];
sx q[1];
rz(-0.82463157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0966245) q[0];
sx q[0];
rz(-1.5929149) q[0];
sx q[0];
rz(1.6629987) q[0];
x q[1];
rz(0.82197676) q[2];
sx q[2];
rz(-0.14227371) q[2];
sx q[2];
rz(1.345944) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6528553) q[1];
sx q[1];
rz(-1.0751546) q[1];
sx q[1];
rz(-1.3929709) q[1];
x q[2];
rz(-0.10004754) q[3];
sx q[3];
rz(-2.2216019) q[3];
sx q[3];
rz(1.5318961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7925988) q[2];
sx q[2];
rz(-0.94453347) q[2];
sx q[2];
rz(-2.3080431) q[2];
rz(-2.870657) q[3];
sx q[3];
rz(-2.0201594) q[3];
sx q[3];
rz(-2.2942395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18868294) q[0];
sx q[0];
rz(-1.8930577) q[0];
sx q[0];
rz(0.46440014) q[0];
rz(2.4480827) q[1];
sx q[1];
rz(-0.92718569) q[1];
sx q[1];
rz(2.4448591) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2385134) q[0];
sx q[0];
rz(-3.0173324) q[0];
sx q[0];
rz(-1.1073698) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7510173) q[2];
sx q[2];
rz(-1.6107133) q[2];
sx q[2];
rz(-0.60933622) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2285959) q[1];
sx q[1];
rz(-1.3218398) q[1];
sx q[1];
rz(-3.1133786) q[1];
rz(-pi) q[2];
rz(-3.1231606) q[3];
sx q[3];
rz(-2.5011809) q[3];
sx q[3];
rz(-0.7225534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36737475) q[2];
sx q[2];
rz(-0.092270277) q[2];
sx q[2];
rz(2.3027159) q[2];
rz(1.687626) q[3];
sx q[3];
rz(-1.5435012) q[3];
sx q[3];
rz(1.4723697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.060318) q[0];
sx q[0];
rz(-1.3723137) q[0];
sx q[0];
rz(-1.6171612) q[0];
rz(0.35628191) q[1];
sx q[1];
rz(-1.7226912) q[1];
sx q[1];
rz(1.3391395) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.08218) q[0];
sx q[0];
rz(-0.97267432) q[0];
sx q[0];
rz(3.0638543) q[0];
x q[1];
rz(0.21276591) q[2];
sx q[2];
rz(-2.8147455) q[2];
sx q[2];
rz(-0.67732993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.410706) q[1];
sx q[1];
rz(-2.2458898) q[1];
sx q[1];
rz(-1.7275586) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.036063866) q[3];
sx q[3];
rz(-1.7129922) q[3];
sx q[3];
rz(-2.9309931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39849207) q[2];
sx q[2];
rz(-2.3091381) q[2];
sx q[2];
rz(-1.2644281) q[2];
rz(0.33665952) q[3];
sx q[3];
rz(-2.2008937) q[3];
sx q[3];
rz(-1.8077067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.7548512) q[0];
sx q[0];
rz(-1.516153) q[0];
sx q[0];
rz(-1.0990912) q[0];
rz(-0.59973888) q[1];
sx q[1];
rz(-2.6987684) q[1];
sx q[1];
rz(-0.59785688) q[1];
rz(0.32156051) q[2];
sx q[2];
rz(-1.9637932) q[2];
sx q[2];
rz(3.1153352) q[2];
rz(3.1118274) q[3];
sx q[3];
rz(-0.53277848) q[3];
sx q[3];
rz(3.0938704) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
