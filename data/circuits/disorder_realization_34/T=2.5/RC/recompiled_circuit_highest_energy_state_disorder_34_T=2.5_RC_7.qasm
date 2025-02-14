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
rz(0.74855411) q[0];
sx q[0];
rz(-1.3286123) q[0];
sx q[0];
rz(0.42088977) q[0];
rz(-0.66840494) q[1];
sx q[1];
rz(4.6894046) q[1];
sx q[1];
rz(8.2735396) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6164857) q[0];
sx q[0];
rz(-2.6455204) q[0];
sx q[0];
rz(-2.708359) q[0];
rz(-pi) q[1];
rz(2.1333063) q[2];
sx q[2];
rz(-1.1739302) q[2];
sx q[2];
rz(-2.4728554) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0033773) q[1];
sx q[1];
rz(-2.1103519) q[1];
sx q[1];
rz(0.92745499) q[1];
x q[2];
rz(0.26845308) q[3];
sx q[3];
rz(-2.2055948) q[3];
sx q[3];
rz(-0.47273794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.067817299) q[2];
sx q[2];
rz(-0.99738085) q[2];
sx q[2];
rz(2.1201521) q[2];
rz(-1.3218309) q[3];
sx q[3];
rz(-0.50522155) q[3];
sx q[3];
rz(-0.57636133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8030871) q[0];
sx q[0];
rz(-1.0271238) q[0];
sx q[0];
rz(-2.8402253) q[0];
rz(-2.001568) q[1];
sx q[1];
rz(-0.81914425) q[1];
sx q[1];
rz(-2.0778621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6607644) q[0];
sx q[0];
rz(-1.6910292) q[0];
sx q[0];
rz(-1.8390435) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76164328) q[2];
sx q[2];
rz(-1.7178665) q[2];
sx q[2];
rz(2.0936793) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7980405) q[1];
sx q[1];
rz(-0.69188373) q[1];
sx q[1];
rz(0.50358285) q[1];
x q[2];
rz(3.0228457) q[3];
sx q[3];
rz(-2.1034965) q[3];
sx q[3];
rz(1.1632196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86338824) q[2];
sx q[2];
rz(-0.74625838) q[2];
sx q[2];
rz(0.14872742) q[2];
rz(1.9637828) q[3];
sx q[3];
rz(-1.1921459) q[3];
sx q[3];
rz(-1.8671794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8389559) q[0];
sx q[0];
rz(-1.1936854) q[0];
sx q[0];
rz(1.0021915) q[0];
rz(1.3305371) q[1];
sx q[1];
rz(-1.1794773) q[1];
sx q[1];
rz(-2.5340714) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9423381) q[0];
sx q[0];
rz(-0.94616854) q[0];
sx q[0];
rz(2.6074431) q[0];
rz(0.27670105) q[2];
sx q[2];
rz(-1.1778129) q[2];
sx q[2];
rz(0.0056841141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3866736) q[1];
sx q[1];
rz(-0.52820871) q[1];
sx q[1];
rz(-0.22294238) q[1];
rz(-2.5030977) q[3];
sx q[3];
rz(-2.3110664) q[3];
sx q[3];
rz(0.80017024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2398296) q[2];
sx q[2];
rz(-2.2679057) q[2];
sx q[2];
rz(-0.98881161) q[2];
rz(1.4766988) q[3];
sx q[3];
rz(-1.3380545) q[3];
sx q[3];
rz(0.57507676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7779509) q[0];
sx q[0];
rz(-0.87894428) q[0];
sx q[0];
rz(-1.2203891) q[0];
rz(1.3149423) q[1];
sx q[1];
rz(-1.5053791) q[1];
sx q[1];
rz(-0.13994089) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24880399) q[0];
sx q[0];
rz(-2.148845) q[0];
sx q[0];
rz(-2.0530353) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5342388) q[2];
sx q[2];
rz(-2.3343041) q[2];
sx q[2];
rz(-1.8497149) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.32711682) q[1];
sx q[1];
rz(-0.45683858) q[1];
sx q[1];
rz(-1.4903699) q[1];
rz(-pi) q[2];
rz(3.0698397) q[3];
sx q[3];
rz(-1.188031) q[3];
sx q[3];
rz(1.07384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4466897) q[2];
sx q[2];
rz(-1.0580772) q[2];
sx q[2];
rz(-3.0042082) q[2];
rz(1.3368227) q[3];
sx q[3];
rz(-1.5788014) q[3];
sx q[3];
rz(3.1409851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7430275) q[0];
sx q[0];
rz(-1.7337357) q[0];
sx q[0];
rz(0.16125691) q[0];
rz(-0.84363168) q[1];
sx q[1];
rz(-2.3944941) q[1];
sx q[1];
rz(1.3074494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3529328) q[0];
sx q[0];
rz(-0.82214117) q[0];
sx q[0];
rz(1.201512) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1990541) q[2];
sx q[2];
rz(-1.2110333) q[2];
sx q[2];
rz(-0.86866405) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2469663) q[1];
sx q[1];
rz(-1.6603396) q[1];
sx q[1];
rz(2.154661) q[1];
x q[2];
rz(-2.4668155) q[3];
sx q[3];
rz(-1.0193079) q[3];
sx q[3];
rz(1.2434208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.97189409) q[2];
sx q[2];
rz(-0.43115386) q[2];
sx q[2];
rz(-0.92740721) q[2];
rz(1.6124604) q[3];
sx q[3];
rz(-1.1011139) q[3];
sx q[3];
rz(0.35759887) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9277495) q[0];
sx q[0];
rz(-0.96927154) q[0];
sx q[0];
rz(-1.4056322) q[0];
rz(1.7270145) q[1];
sx q[1];
rz(-0.79726338) q[1];
sx q[1];
rz(2.3419211) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1522923) q[0];
sx q[0];
rz(-1.5591027) q[0];
sx q[0];
rz(-1.508827) q[0];
rz(-pi) q[1];
x q[1];
rz(1.437101) q[2];
sx q[2];
rz(-1.8779199) q[2];
sx q[2];
rz(1.4442867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.97669125) q[1];
sx q[1];
rz(-1.8215239) q[1];
sx q[1];
rz(1.2198971) q[1];
rz(-pi) q[2];
rz(-0.40414895) q[3];
sx q[3];
rz(-1.3938483) q[3];
sx q[3];
rz(1.1791942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8863525) q[2];
sx q[2];
rz(-0.73840529) q[2];
sx q[2];
rz(-2.3389471) q[2];
rz(0.32043996) q[3];
sx q[3];
rz(-1.3405864) q[3];
sx q[3];
rz(1.505544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36444148) q[0];
sx q[0];
rz(-0.027712263) q[0];
sx q[0];
rz(-3.1112772) q[0];
rz(-1.3069356) q[1];
sx q[1];
rz(-1.0825284) q[1];
sx q[1];
rz(-2.511715) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2455642) q[0];
sx q[0];
rz(-1.7587852) q[0];
sx q[0];
rz(2.7134368) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81923072) q[2];
sx q[2];
rz(-2.5360245) q[2];
sx q[2];
rz(-1.2450964) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3105433) q[1];
sx q[1];
rz(-1.2782974) q[1];
sx q[1];
rz(-1.7834375) q[1];
rz(-0.23655741) q[3];
sx q[3];
rz(-0.45522296) q[3];
sx q[3];
rz(0.61456028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.031875413) q[2];
sx q[2];
rz(-2.4745291) q[2];
sx q[2];
rz(1.9943705) q[2];
rz(2.4540497) q[3];
sx q[3];
rz(-2.5496428) q[3];
sx q[3];
rz(-1.8182925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3654093) q[0];
sx q[0];
rz(-0.23615806) q[0];
sx q[0];
rz(-0.48015204) q[0];
rz(-0.40223739) q[1];
sx q[1];
rz(-1.4996585) q[1];
sx q[1];
rz(0.38280907) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7744336) q[0];
sx q[0];
rz(-1.8978137) q[0];
sx q[0];
rz(2.9190382) q[0];
x q[1];
rz(1.2023964) q[2];
sx q[2];
rz(-2.830626) q[2];
sx q[2];
rz(0.63705618) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76012145) q[1];
sx q[1];
rz(-0.49182004) q[1];
sx q[1];
rz(-0.90729883) q[1];
x q[2];
rz(2.945369) q[3];
sx q[3];
rz(-1.5999402) q[3];
sx q[3];
rz(1.9118229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39168921) q[2];
sx q[2];
rz(-2.4112371) q[2];
sx q[2];
rz(-0.35823092) q[2];
rz(0.41424888) q[3];
sx q[3];
rz(-1.0284938) q[3];
sx q[3];
rz(-2.7481368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8527894) q[0];
sx q[0];
rz(-1.2816592) q[0];
sx q[0];
rz(2.7237256) q[0];
rz(1.6425543) q[1];
sx q[1];
rz(-0.42048979) q[1];
sx q[1];
rz(-2.5299759) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9238511) q[0];
sx q[0];
rz(-1.3425266) q[0];
sx q[0];
rz(-1.7817253) q[0];
rz(-pi) q[1];
rz(0.29199227) q[2];
sx q[2];
rz(-2.9270009) q[2];
sx q[2];
rz(0.24927441) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4222718) q[1];
sx q[1];
rz(-2.5818655) q[1];
sx q[1];
rz(-1.0258963) q[1];
rz(0.34682746) q[3];
sx q[3];
rz(-1.6478103) q[3];
sx q[3];
rz(-1.6517757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9548107) q[2];
sx q[2];
rz(-1.45767) q[2];
sx q[2];
rz(-2.7853454) q[2];
rz(2.6567843) q[3];
sx q[3];
rz(-1.7907413) q[3];
sx q[3];
rz(-0.029732186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.818882) q[0];
sx q[0];
rz(-2.0368545) q[0];
sx q[0];
rz(-1.4622965) q[0];
rz(-2.7541584) q[1];
sx q[1];
rz(-2.2144364) q[1];
sx q[1];
rz(-0.80518728) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0793251) q[0];
sx q[0];
rz(-2.1342417) q[0];
sx q[0];
rz(0.84765537) q[0];
rz(2.3024955) q[2];
sx q[2];
rz(-2.5658414) q[2];
sx q[2];
rz(2.8481399) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8987052) q[1];
sx q[1];
rz(-1.4243898) q[1];
sx q[1];
rz(2.8553748) q[1];
rz(1.8191387) q[3];
sx q[3];
rz(-1.9232009) q[3];
sx q[3];
rz(-0.067841522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33599535) q[2];
sx q[2];
rz(-2.5184641) q[2];
sx q[2];
rz(-1.222329) q[2];
rz(2.3313816) q[3];
sx q[3];
rz(-2.2170292) q[3];
sx q[3];
rz(-1.5626102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7280818) q[0];
sx q[0];
rz(-1.6376729) q[0];
sx q[0];
rz(-1.6117657) q[0];
rz(-1.8660846) q[1];
sx q[1];
rz(-2.7559912) q[1];
sx q[1];
rz(-1.4082946) q[1];
rz(-1.9553984) q[2];
sx q[2];
rz(-1.4438965) q[2];
sx q[2];
rz(0.039197103) q[2];
rz(2.1199216) q[3];
sx q[3];
rz(-0.74629489) q[3];
sx q[3];
rz(-1.3003579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
