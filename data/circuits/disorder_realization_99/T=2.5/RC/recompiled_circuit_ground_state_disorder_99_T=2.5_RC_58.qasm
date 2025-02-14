OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34243256) q[0];
sx q[0];
rz(-0.39781308) q[0];
sx q[0];
rz(-1.2678658) q[0];
rz(-0.59029382) q[1];
sx q[1];
rz(-0.78650147) q[1];
sx q[1];
rz(1.9222577) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0578227) q[0];
sx q[0];
rz(-1.738213) q[0];
sx q[0];
rz(0.98982248) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27094312) q[2];
sx q[2];
rz(-0.56519714) q[2];
sx q[2];
rz(-0.2172367) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.089069033) q[1];
sx q[1];
rz(-0.34862754) q[1];
sx q[1];
rz(2.2050956) q[1];
rz(-pi) q[2];
rz(-0.91359992) q[3];
sx q[3];
rz(-1.3845647) q[3];
sx q[3];
rz(-2.5600764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5775602) q[2];
sx q[2];
rz(-1.6200248) q[2];
sx q[2];
rz(-1.7731898) q[2];
rz(-2.2300301) q[3];
sx q[3];
rz(-1.3699968) q[3];
sx q[3];
rz(0.024356775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1034705) q[0];
sx q[0];
rz(-2.7439674) q[0];
sx q[0];
rz(0.11904112) q[0];
rz(2.0114404) q[1];
sx q[1];
rz(-0.96769133) q[1];
sx q[1];
rz(-1.7134604) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6551681) q[0];
sx q[0];
rz(-2.0858686) q[0];
sx q[0];
rz(1.0977655) q[0];
rz(-0.4006673) q[2];
sx q[2];
rz(-1.180449) q[2];
sx q[2];
rz(-2.0290749) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5806313) q[1];
sx q[1];
rz(-2.033618) q[1];
sx q[1];
rz(2.9918519) q[1];
x q[2];
rz(-0.83934966) q[3];
sx q[3];
rz(-2.7003059) q[3];
sx q[3];
rz(-1.4660783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52593645) q[2];
sx q[2];
rz(-1.6776513) q[2];
sx q[2];
rz(-0.76457912) q[2];
rz(0.48314759) q[3];
sx q[3];
rz(-1.8254447) q[3];
sx q[3];
rz(-2.29276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1044384) q[0];
sx q[0];
rz(-0.26532441) q[0];
sx q[0];
rz(-1.4720488) q[0];
rz(0.56023359) q[1];
sx q[1];
rz(-1.3872223) q[1];
sx q[1];
rz(1.4405506) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9870696) q[0];
sx q[0];
rz(-1.8049585) q[0];
sx q[0];
rz(2.9375141) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5650827) q[2];
sx q[2];
rz(-2.2457321) q[2];
sx q[2];
rz(2.142148) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6145103) q[1];
sx q[1];
rz(-2.0918814) q[1];
sx q[1];
rz(-0.39343843) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3769242) q[3];
sx q[3];
rz(-2.8199785) q[3];
sx q[3];
rz(-2.4177616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6560087) q[2];
sx q[2];
rz(-1.4914844) q[2];
sx q[2];
rz(-0.34240016) q[2];
rz(1.7229236) q[3];
sx q[3];
rz(-0.72903052) q[3];
sx q[3];
rz(0.5595783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2320084) q[0];
sx q[0];
rz(-2.7689458) q[0];
sx q[0];
rz(-1.5268071) q[0];
rz(-0.80742637) q[1];
sx q[1];
rz(-1.5734943) q[1];
sx q[1];
rz(1.2015013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1781834) q[0];
sx q[0];
rz(-0.82287517) q[0];
sx q[0];
rz(-1.9848518) q[0];
rz(2.3551201) q[2];
sx q[2];
rz(-0.82068372) q[2];
sx q[2];
rz(0.72696668) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0240942) q[1];
sx q[1];
rz(-1.9650241) q[1];
sx q[1];
rz(-2.3181163) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19530794) q[3];
sx q[3];
rz(-2.8805519) q[3];
sx q[3];
rz(1.1782427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80369192) q[2];
sx q[2];
rz(-0.98946977) q[2];
sx q[2];
rz(1.2376415) q[2];
rz(-1.8303653) q[3];
sx q[3];
rz(-1.6236191) q[3];
sx q[3];
rz(-1.9633912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.118498) q[0];
sx q[0];
rz(-1.3730405) q[0];
sx q[0];
rz(2.232724) q[0];
rz(0.88242775) q[1];
sx q[1];
rz(-0.71417037) q[1];
sx q[1];
rz(0.73208255) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013618795) q[0];
sx q[0];
rz(-2.2165856) q[0];
sx q[0];
rz(2.5396106) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3734748) q[2];
sx q[2];
rz(-1.7752194) q[2];
sx q[2];
rz(-2.8724175) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4914507) q[1];
sx q[1];
rz(-0.60381266) q[1];
sx q[1];
rz(0.73765386) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77712147) q[3];
sx q[3];
rz(-0.91360559) q[3];
sx q[3];
rz(-1.4025721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.070179209) q[2];
sx q[2];
rz(-2.3290403) q[2];
sx q[2];
rz(1.8522235) q[2];
rz(1.4687126) q[3];
sx q[3];
rz(-1.9332935) q[3];
sx q[3];
rz(-2.7220791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5656723) q[0];
sx q[0];
rz(-1.1812295) q[0];
sx q[0];
rz(-1.1015724) q[0];
rz(2.9167602) q[1];
sx q[1];
rz(-1.79554) q[1];
sx q[1];
rz(-0.98519957) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3168047) q[0];
sx q[0];
rz(-1.7508306) q[0];
sx q[0];
rz(0.92702528) q[0];
rz(-pi) q[1];
rz(-2.6654319) q[2];
sx q[2];
rz(-0.67103065) q[2];
sx q[2];
rz(-0.76045017) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0538865) q[1];
sx q[1];
rz(-2.3017677) q[1];
sx q[1];
rz(2.5165416) q[1];
rz(-pi) q[2];
rz(1.5456301) q[3];
sx q[3];
rz(-1.0428793) q[3];
sx q[3];
rz(0.86564579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78099403) q[2];
sx q[2];
rz(-2.4807319) q[2];
sx q[2];
rz(-1.194427) q[2];
rz(0.099695168) q[3];
sx q[3];
rz(-1.7710779) q[3];
sx q[3];
rz(-2.1626507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7128971) q[0];
sx q[0];
rz(-0.45658699) q[0];
sx q[0];
rz(-1.1140484) q[0];
rz(1.7440589) q[1];
sx q[1];
rz(-1.3797398) q[1];
sx q[1];
rz(2.8278415) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6614051) q[0];
sx q[0];
rz(-1.620694) q[0];
sx q[0];
rz(-1.2779092) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4391104) q[2];
sx q[2];
rz(-1.4073616) q[2];
sx q[2];
rz(1.2800467) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.75668979) q[1];
sx q[1];
rz(-0.80763615) q[1];
sx q[1];
rz(-0.4203504) q[1];
x q[2];
rz(2.8242495) q[3];
sx q[3];
rz(-2.4306524) q[3];
sx q[3];
rz(1.8985572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.067001192) q[2];
sx q[2];
rz(-0.40412298) q[2];
sx q[2];
rz(0.60603777) q[2];
rz(-2.0634985) q[3];
sx q[3];
rz(-0.86229101) q[3];
sx q[3];
rz(-1.7830474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17230497) q[0];
sx q[0];
rz(-2.2242039) q[0];
sx q[0];
rz(-2.0148) q[0];
rz(2.2276095) q[1];
sx q[1];
rz(-0.4387478) q[1];
sx q[1];
rz(2.9235358) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.954531) q[0];
sx q[0];
rz(-2.0045223) q[0];
sx q[0];
rz(-2.0824672) q[0];
x q[1];
rz(0.78812353) q[2];
sx q[2];
rz(-1.0255073) q[2];
sx q[2];
rz(-2.9719549) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40849028) q[1];
sx q[1];
rz(-0.6576076) q[1];
sx q[1];
rz(-1.3642743) q[1];
rz(-pi) q[2];
rz(3.0239437) q[3];
sx q[3];
rz(-1.7075286) q[3];
sx q[3];
rz(-2.4916388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5128936) q[2];
sx q[2];
rz(-2.4662374) q[2];
sx q[2];
rz(2.8524354) q[2];
rz(-2.1469927) q[3];
sx q[3];
rz(-1.1887487) q[3];
sx q[3];
rz(-2.6836256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4299803) q[0];
sx q[0];
rz(-2.0531605) q[0];
sx q[0];
rz(2.3751538) q[0];
rz(1.6038766) q[1];
sx q[1];
rz(-1.3130554) q[1];
sx q[1];
rz(2.2235353) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9452451) q[0];
sx q[0];
rz(-0.96789384) q[0];
sx q[0];
rz(-2.5640998) q[0];
x q[1];
rz(-0.70785825) q[2];
sx q[2];
rz(-0.22032693) q[2];
sx q[2];
rz(-2.3489621) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.854771) q[1];
sx q[1];
rz(-2.3769551) q[1];
sx q[1];
rz(-2.7809793) q[1];
rz(-0.10580801) q[3];
sx q[3];
rz(-1.322016) q[3];
sx q[3];
rz(2.5073568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90091577) q[2];
sx q[2];
rz(-0.46611163) q[2];
sx q[2];
rz(-0.051699836) q[2];
rz(-2.2895571) q[3];
sx q[3];
rz(-1.4308735) q[3];
sx q[3];
rz(-2.8813072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.32605115) q[0];
sx q[0];
rz(-1.4651848) q[0];
sx q[0];
rz(-0.091212243) q[0];
rz(0.7971898) q[1];
sx q[1];
rz(-1.3858831) q[1];
sx q[1];
rz(-0.43112722) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9190311) q[0];
sx q[0];
rz(-3.0662144) q[0];
sx q[0];
rz(-0.68392174) q[0];
rz(-pi) q[1];
rz(-1.4412142) q[2];
sx q[2];
rz(-0.65276546) q[2];
sx q[2];
rz(1.5127986) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3865652) q[1];
sx q[1];
rz(-1.3530827) q[1];
sx q[1];
rz(-2.6854807) q[1];
x q[2];
rz(0.98072106) q[3];
sx q[3];
rz(-2.4517165) q[3];
sx q[3];
rz(1.7074613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5690696) q[2];
sx q[2];
rz(-1.5745682) q[2];
sx q[2];
rz(1.265906) q[2];
rz(-1.2931394) q[3];
sx q[3];
rz(-2.0796516) q[3];
sx q[3];
rz(-2.7354447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1311998) q[0];
sx q[0];
rz(-2.4507903) q[0];
sx q[0];
rz(-2.3660085) q[0];
rz(-2.5776183) q[1];
sx q[1];
rz(-2.0717944) q[1];
sx q[1];
rz(2.0800128) q[1];
rz(1.1376913) q[2];
sx q[2];
rz(-1.4972996) q[2];
sx q[2];
rz(1.0232915) q[2];
rz(-0.98714491) q[3];
sx q[3];
rz(-2.5844283) q[3];
sx q[3];
rz(1.9682932) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
