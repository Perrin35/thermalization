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
rz(-1.0492078) q[0];
sx q[0];
rz(7.0424289) q[0];
sx q[0];
rz(9.0483604) q[0];
rz(1.5578101) q[1];
sx q[1];
rz(-1.0695142) q[1];
sx q[1];
rz(0.95626107) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1726748) q[0];
sx q[0];
rz(-0.71828523) q[0];
sx q[0];
rz(1.0802313) q[0];
rz(-1.3772493) q[2];
sx q[2];
rz(-1.1551394) q[2];
sx q[2];
rz(3.1061098) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8531283) q[1];
sx q[1];
rz(-0.83544105) q[1];
sx q[1];
rz(-2.9771027) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8703635) q[3];
sx q[3];
rz(-1.3406957) q[3];
sx q[3];
rz(2.7903737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.817953) q[2];
sx q[2];
rz(-2.2223667) q[2];
sx q[2];
rz(-1.2032571) q[2];
rz(2.1008927) q[3];
sx q[3];
rz(-2.2668656) q[3];
sx q[3];
rz(0.11999764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.1131209) q[0];
sx q[0];
rz(-2.5682243) q[0];
sx q[0];
rz(-2.0290802) q[0];
rz(-1.6587229) q[1];
sx q[1];
rz(-0.590938) q[1];
sx q[1];
rz(-0.15731752) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4388889) q[0];
sx q[0];
rz(-1.7458785) q[0];
sx q[0];
rz(-0.69467993) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21944616) q[2];
sx q[2];
rz(-2.7601961) q[2];
sx q[2];
rz(-1.4666605) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.29780218) q[1];
sx q[1];
rz(-0.71759598) q[1];
sx q[1];
rz(1.1949431) q[1];
rz(0.3213082) q[3];
sx q[3];
rz(-0.23443334) q[3];
sx q[3];
rz(-3.0911764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6436254) q[2];
sx q[2];
rz(-0.47933856) q[2];
sx q[2];
rz(0.41110006) q[2];
rz(-0.57605612) q[3];
sx q[3];
rz(-2.1274302) q[3];
sx q[3];
rz(-2.1435553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882947) q[0];
sx q[0];
rz(-1.0030712) q[0];
sx q[0];
rz(2.4004747) q[0];
rz(-1.6916212) q[1];
sx q[1];
rz(-1.8284109) q[1];
sx q[1];
rz(-0.57201874) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93942552) q[0];
sx q[0];
rz(-1.8644445) q[0];
sx q[0];
rz(2.9333326) q[0];
rz(2.8764804) q[2];
sx q[2];
rz(-1.371061) q[2];
sx q[2];
rz(2.0743362) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.95840824) q[1];
sx q[1];
rz(-2.5377877) q[1];
sx q[1];
rz(-2.3314657) q[1];
x q[2];
rz(-2.9111262) q[3];
sx q[3];
rz(-1.8982045) q[3];
sx q[3];
rz(-0.62124204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5570306) q[2];
sx q[2];
rz(-1.9494373) q[2];
sx q[2];
rz(-2.3438047) q[2];
rz(1.8320463) q[3];
sx q[3];
rz(-1.2582425) q[3];
sx q[3];
rz(-1.4057188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0107467) q[0];
sx q[0];
rz(-2.2444785) q[0];
sx q[0];
rz(-0.36706269) q[0];
rz(-1.2182073) q[1];
sx q[1];
rz(-1.7299078) q[1];
sx q[1];
rz(-0.22148111) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2838615) q[0];
sx q[0];
rz(-0.22131187) q[0];
sx q[0];
rz(-2.8279634) q[0];
rz(0.34299739) q[2];
sx q[2];
rz(-2.7532737) q[2];
sx q[2];
rz(-1.266154) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96326107) q[1];
sx q[1];
rz(-1.8823922) q[1];
sx q[1];
rz(-2.7072886) q[1];
rz(-pi) q[2];
rz(-0.4226004) q[3];
sx q[3];
rz(-0.56767332) q[3];
sx q[3];
rz(-1.467553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9941142) q[2];
sx q[2];
rz(-1.9693815) q[2];
sx q[2];
rz(-1.3383024) q[2];
rz(2.639468) q[3];
sx q[3];
rz(-0.10052557) q[3];
sx q[3];
rz(-0.24371915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5153656) q[0];
sx q[0];
rz(-0.60972649) q[0];
sx q[0];
rz(1.5355661) q[0];
rz(-3.0961127) q[1];
sx q[1];
rz(-1.7605503) q[1];
sx q[1];
rz(2.6754726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.202318) q[0];
sx q[0];
rz(-2.8554248) q[0];
sx q[0];
rz(-1.0658468) q[0];
x q[1];
rz(2.8117198) q[2];
sx q[2];
rz(-0.95392841) q[2];
sx q[2];
rz(-2.5170779) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.921195) q[1];
sx q[1];
rz(-0.092893727) q[1];
sx q[1];
rz(2.9542406) q[1];
rz(-pi) q[2];
rz(-0.65069549) q[3];
sx q[3];
rz(-1.078372) q[3];
sx q[3];
rz(0.7975815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8972299) q[2];
sx q[2];
rz(-1.676061) q[2];
sx q[2];
rz(-2.3587295) q[2];
rz(3.1252981) q[3];
sx q[3];
rz(-1.3085082) q[3];
sx q[3];
rz(2.079594) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8472854) q[0];
sx q[0];
rz(-2.6562302) q[0];
sx q[0];
rz(-0.37931994) q[0];
rz(2.8324221) q[1];
sx q[1];
rz(-1.199017) q[1];
sx q[1];
rz(-0.22383037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17344089) q[0];
sx q[0];
rz(-1.26014) q[0];
sx q[0];
rz(-2.2947427) q[0];
rz(3.0293944) q[2];
sx q[2];
rz(-1.7609481) q[2];
sx q[2];
rz(-0.79922215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1081738) q[1];
sx q[1];
rz(-1.6774607) q[1];
sx q[1];
rz(-1.4640385) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.623935) q[3];
sx q[3];
rz(-2.3316962) q[3];
sx q[3];
rz(-0.018190688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0763187) q[2];
sx q[2];
rz(-1.9741917) q[2];
sx q[2];
rz(-2.8889612) q[2];
rz(-0.97918716) q[3];
sx q[3];
rz(-0.88117176) q[3];
sx q[3];
rz(0.34669909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52794367) q[0];
sx q[0];
rz(-0.76512965) q[0];
sx q[0];
rz(1.2921523) q[0];
rz(-2.6076803) q[1];
sx q[1];
rz(-1.4521867) q[1];
sx q[1];
rz(-2.7814878) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6394516) q[0];
sx q[0];
rz(-1.2107673) q[0];
sx q[0];
rz(-3.0675824) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5582997) q[2];
sx q[2];
rz(-0.55634004) q[2];
sx q[2];
rz(-0.46890989) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95197612) q[1];
sx q[1];
rz(-1.8120753) q[1];
sx q[1];
rz(0.56778384) q[1];
rz(-3.1084314) q[3];
sx q[3];
rz(-1.6735184) q[3];
sx q[3];
rz(-2.9364613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1503633) q[2];
sx q[2];
rz(-1.6232619) q[2];
sx q[2];
rz(-0.72597996) q[2];
rz(0.43068543) q[3];
sx q[3];
rz(-0.78023282) q[3];
sx q[3];
rz(-0.45346692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78366572) q[0];
sx q[0];
rz(-1.6538606) q[0];
sx q[0];
rz(-0.7861535) q[0];
rz(-0.17732009) q[1];
sx q[1];
rz(-1.0847849) q[1];
sx q[1];
rz(1.0160944) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7854016) q[0];
sx q[0];
rz(-1.8520266) q[0];
sx q[0];
rz(-3.0321971) q[0];
x q[1];
rz(0.21206484) q[2];
sx q[2];
rz(-1.4432505) q[2];
sx q[2];
rz(-2.2812592) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7238771) q[1];
sx q[1];
rz(-2.2347365) q[1];
sx q[1];
rz(1.5491897) q[1];
rz(1.949259) q[3];
sx q[3];
rz(-2.2871823) q[3];
sx q[3];
rz(2.5797957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59755406) q[2];
sx q[2];
rz(-0.95557135) q[2];
sx q[2];
rz(2.3717144) q[2];
rz(0.17635135) q[3];
sx q[3];
rz(-1.4498815) q[3];
sx q[3];
rz(-0.31340733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65799323) q[0];
sx q[0];
rz(-3.0577116) q[0];
sx q[0];
rz(-1.0461079) q[0];
rz(-1.5484035) q[1];
sx q[1];
rz(-1.7240588) q[1];
sx q[1];
rz(1.2215325) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0451121) q[0];
sx q[0];
rz(-1.6754221) q[0];
sx q[0];
rz(-1.2136995) q[0];
x q[1];
rz(0.1510001) q[2];
sx q[2];
rz(-1.2556453) q[2];
sx q[2];
rz(2.0133108) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9323349) q[1];
sx q[1];
rz(-2.845721) q[1];
sx q[1];
rz(-2.3274755) q[1];
x q[2];
rz(-1.5658978) q[3];
sx q[3];
rz(-1.6143482) q[3];
sx q[3];
rz(2.5033308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.017642411) q[2];
sx q[2];
rz(-1.4175748) q[2];
sx q[2];
rz(1.9564015) q[2];
rz(-0.3398529) q[3];
sx q[3];
rz(-0.9797107) q[3];
sx q[3];
rz(0.11782304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3621984) q[0];
sx q[0];
rz(-1.3267936) q[0];
sx q[0];
rz(-2.4438044) q[0];
rz(-2.4367874) q[1];
sx q[1];
rz(-1.1230527) q[1];
sx q[1];
rz(0.25308457) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11624434) q[0];
sx q[0];
rz(-1.3379813) q[0];
sx q[0];
rz(3.065585) q[0];
rz(0.63941892) q[2];
sx q[2];
rz(-0.76335159) q[2];
sx q[2];
rz(-1.167041) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.53348161) q[1];
sx q[1];
rz(-0.93431707) q[1];
sx q[1];
rz(0.57884207) q[1];
rz(-1.8277824) q[3];
sx q[3];
rz(-0.8567613) q[3];
sx q[3];
rz(-2.11588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7207328) q[2];
sx q[2];
rz(-1.3497817) q[2];
sx q[2];
rz(0.68812686) q[2];
rz(0.94528919) q[3];
sx q[3];
rz(-1.9663845) q[3];
sx q[3];
rz(0.92760408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59104334) q[0];
sx q[0];
rz(-1.782438) q[0];
sx q[0];
rz(1.8989643) q[0];
rz(-0.91724829) q[1];
sx q[1];
rz(-1.6684253) q[1];
sx q[1];
rz(-2.565276) q[1];
rz(2.5653432) q[2];
sx q[2];
rz(-1.2628308) q[2];
sx q[2];
rz(-2.2028883) q[2];
rz(-1.7797021) q[3];
sx q[3];
rz(-2.7068797) q[3];
sx q[3];
rz(-2.5754365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
