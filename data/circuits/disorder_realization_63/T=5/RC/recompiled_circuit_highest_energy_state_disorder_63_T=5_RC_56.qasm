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
rz(2.0923848) q[0];
sx q[0];
rz(-0.75924358) q[0];
sx q[0];
rz(2.7651751) q[0];
rz(1.5578101) q[1];
sx q[1];
rz(-1.0695142) q[1];
sx q[1];
rz(0.95626107) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1726748) q[0];
sx q[0];
rz(-2.4233074) q[0];
sx q[0];
rz(-1.0802313) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7189288) q[2];
sx q[2];
rz(-1.7476805) q[2];
sx q[2];
rz(-1.4563349) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9701582) q[1];
sx q[1];
rz(-1.6925319) q[1];
sx q[1];
rz(2.3129169) q[1];
x q[2];
rz(0.27122916) q[3];
sx q[3];
rz(-1.3406957) q[3];
sx q[3];
rz(2.7903737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.817953) q[2];
sx q[2];
rz(-2.2223667) q[2];
sx q[2];
rz(1.2032571) q[2];
rz(-1.0407) q[3];
sx q[3];
rz(-2.2668656) q[3];
sx q[3];
rz(-3.021595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1131209) q[0];
sx q[0];
rz(-0.57336837) q[0];
sx q[0];
rz(1.1125125) q[0];
rz(1.6587229) q[1];
sx q[1];
rz(-2.5506546) q[1];
sx q[1];
rz(2.9842751) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1293662) q[0];
sx q[0];
rz(-2.2528045) q[0];
sx q[0];
rz(-1.3444882) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9221465) q[2];
sx q[2];
rz(-0.38139653) q[2];
sx q[2];
rz(-1.6749322) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.98401947) q[1];
sx q[1];
rz(-1.3270151) q[1];
sx q[1];
rz(-0.88884377) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4955161) q[3];
sx q[3];
rz(-1.7930248) q[3];
sx q[3];
rz(0.27929515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4979672) q[2];
sx q[2];
rz(-0.47933856) q[2];
sx q[2];
rz(-0.41110006) q[2];
rz(-2.5655365) q[3];
sx q[3];
rz(-2.1274302) q[3];
sx q[3];
rz(2.1435553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75329798) q[0];
sx q[0];
rz(-2.1385215) q[0];
sx q[0];
rz(0.74111795) q[0];
rz(1.4499715) q[1];
sx q[1];
rz(-1.8284109) q[1];
sx q[1];
rz(-0.57201874) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30878371) q[0];
sx q[0];
rz(-0.35825142) q[0];
sx q[0];
rz(-0.97106309) q[0];
rz(-1.3640312) q[2];
sx q[2];
rz(-1.3110779) q[2];
sx q[2];
rz(-2.5842359) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1831844) q[1];
sx q[1];
rz(-2.5377877) q[1];
sx q[1];
rz(-2.3314657) q[1];
x q[2];
rz(2.9111262) q[3];
sx q[3];
rz(-1.8982045) q[3];
sx q[3];
rz(0.62124204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.584562) q[2];
sx q[2];
rz(-1.1921554) q[2];
sx q[2];
rz(-2.3438047) q[2];
rz(1.3095464) q[3];
sx q[3];
rz(-1.2582425) q[3];
sx q[3];
rz(1.4057188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0107467) q[0];
sx q[0];
rz(-0.89711419) q[0];
sx q[0];
rz(2.77453) q[0];
rz(-1.2182073) q[1];
sx q[1];
rz(-1.4116849) q[1];
sx q[1];
rz(0.22148111) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7349827) q[0];
sx q[0];
rz(-1.5030229) q[0];
sx q[0];
rz(2.9307515) q[0];
rz(-pi) q[1];
rz(0.367737) q[2];
sx q[2];
rz(-1.6984816) q[2];
sx q[2];
rz(-2.5177296) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.023497009) q[1];
sx q[1];
rz(-0.52866908) q[1];
sx q[1];
rz(0.65331991) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52677299) q[3];
sx q[3];
rz(-1.793141) q[3];
sx q[3];
rz(-0.46565817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.14747846) q[2];
sx q[2];
rz(-1.1722112) q[2];
sx q[2];
rz(1.3383024) q[2];
rz(0.5021247) q[3];
sx q[3];
rz(-0.10052557) q[3];
sx q[3];
rz(0.24371915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5153656) q[0];
sx q[0];
rz(-2.5318662) q[0];
sx q[0];
rz(1.5355661) q[0];
rz(0.045479927) q[1];
sx q[1];
rz(-1.3810424) q[1];
sx q[1];
rz(0.46612003) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2855466) q[0];
sx q[0];
rz(-1.4338126) q[0];
sx q[0];
rz(-1.3187506) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32987288) q[2];
sx q[2];
rz(-2.1876642) q[2];
sx q[2];
rz(-0.62451476) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.032253232) q[1];
sx q[1];
rz(-1.6620599) q[1];
sx q[1];
rz(1.5534459) q[1];
x q[2];
rz(-0.97753559) q[3];
sx q[3];
rz(-2.1338507) q[3];
sx q[3];
rz(-2.7138591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8972299) q[2];
sx q[2];
rz(-1.676061) q[2];
sx q[2];
rz(2.3587295) q[2];
rz(-3.1252981) q[3];
sx q[3];
rz(-1.3085082) q[3];
sx q[3];
rz(1.0619987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29430729) q[0];
sx q[0];
rz(-0.48536244) q[0];
sx q[0];
rz(0.37931994) q[0];
rz(2.8324221) q[1];
sx q[1];
rz(-1.199017) q[1];
sx q[1];
rz(2.9177623) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4113385) q[0];
sx q[0];
rz(-0.77654949) q[0];
sx q[0];
rz(-1.1194487) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0293944) q[2];
sx q[2];
rz(-1.7609481) q[2];
sx q[2];
rz(-2.3423705) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47403112) q[1];
sx q[1];
rz(-1.6769451) q[1];
sx q[1];
rz(-0.1072704) q[1];
x q[2];
rz(-1.0204487) q[3];
sx q[3];
rz(-2.1990877) q[3];
sx q[3];
rz(2.3165143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0763187) q[2];
sx q[2];
rz(-1.9741917) q[2];
sx q[2];
rz(2.8889612) q[2];
rz(2.1624055) q[3];
sx q[3];
rz(-0.88117176) q[3];
sx q[3];
rz(0.34669909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.613649) q[0];
sx q[0];
rz(-0.76512965) q[0];
sx q[0];
rz(-1.2921523) q[0];
rz(-0.53391236) q[1];
sx q[1];
rz(-1.4521867) q[1];
sx q[1];
rz(2.7814878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6394516) q[0];
sx q[0];
rz(-1.9308254) q[0];
sx q[0];
rz(3.0675824) q[0];
rz(1.8890284) q[2];
sx q[2];
rz(-1.1063853) q[2];
sx q[2];
rz(-2.0385252) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77006631) q[1];
sx q[1];
rz(-2.1202104) q[1];
sx q[1];
rz(-1.2868164) q[1];
rz(1.259616) q[3];
sx q[3];
rz(-0.10792416) q[3];
sx q[3];
rz(-0.51801658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1503633) q[2];
sx q[2];
rz(-1.6232619) q[2];
sx q[2];
rz(2.4156127) q[2];
rz(2.7109072) q[3];
sx q[3];
rz(-0.78023282) q[3];
sx q[3];
rz(0.45346692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3579269) q[0];
sx q[0];
rz(-1.6538606) q[0];
sx q[0];
rz(-0.7861535) q[0];
rz(-2.9642726) q[1];
sx q[1];
rz(-2.0568078) q[1];
sx q[1];
rz(1.0160944) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1209512) q[0];
sx q[0];
rz(-0.30123152) q[0];
sx q[0];
rz(-1.2094638) q[0];
x q[1];
rz(2.9295278) q[2];
sx q[2];
rz(-1.6983422) q[2];
sx q[2];
rz(0.86033347) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1397649) q[1];
sx q[1];
rz(-1.5537801) q[1];
sx q[1];
rz(2.4775392) q[1];
rz(-pi) q[2];
x q[2];
rz(1.949259) q[3];
sx q[3];
rz(-0.85441033) q[3];
sx q[3];
rz(-2.5797957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59755406) q[2];
sx q[2];
rz(-0.95557135) q[2];
sx q[2];
rz(2.3717144) q[2];
rz(2.9652413) q[3];
sx q[3];
rz(-1.4498815) q[3];
sx q[3];
rz(-2.8281853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4835994) q[0];
sx q[0];
rz(-0.083881065) q[0];
sx q[0];
rz(1.0461079) q[0];
rz(-1.5931891) q[1];
sx q[1];
rz(-1.7240588) q[1];
sx q[1];
rz(1.9200602) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6283327) q[0];
sx q[0];
rz(-1.9258537) q[0];
sx q[0];
rz(3.0299788) q[0];
rz(-pi) q[1];
rz(2.0030973) q[2];
sx q[2];
rz(-0.34836787) q[2];
sx q[2];
rz(1.5846202) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0960779) q[1];
sx q[1];
rz(-1.3692665) q[1];
sx q[1];
rz(1.3526826) q[1];
rz(0.11193377) q[3];
sx q[3];
rz(-0.043826274) q[3];
sx q[3];
rz(-0.52622139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1239502) q[2];
sx q[2];
rz(-1.7240179) q[2];
sx q[2];
rz(-1.9564015) q[2];
rz(-0.3398529) q[3];
sx q[3];
rz(-2.161882) q[3];
sx q[3];
rz(-0.11782304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3621984) q[0];
sx q[0];
rz(-1.3267936) q[0];
sx q[0];
rz(-0.69778824) q[0];
rz(-2.4367874) q[1];
sx q[1];
rz(-2.0185399) q[1];
sx q[1];
rz(2.8885081) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9390209) q[0];
sx q[0];
rz(-0.24469412) q[0];
sx q[0];
rz(-1.2608725) q[0];
rz(-pi) q[1];
rz(-0.65480755) q[2];
sx q[2];
rz(-1.1455481) q[2];
sx q[2];
rz(0.08928334) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8423137) q[1];
sx q[1];
rz(-0.83228534) q[1];
sx q[1];
rz(0.93361093) q[1];
x q[2];
rz(2.4109957) q[3];
sx q[3];
rz(-1.7640778) q[3];
sx q[3];
rz(2.766942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42085984) q[2];
sx q[2];
rz(-1.791811) q[2];
sx q[2];
rz(2.4534658) q[2];
rz(-2.1963035) q[3];
sx q[3];
rz(-1.1752081) q[3];
sx q[3];
rz(-0.92760408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59104334) q[0];
sx q[0];
rz(-1.3591546) q[0];
sx q[0];
rz(-1.2426283) q[0];
rz(0.91724829) q[1];
sx q[1];
rz(-1.4731673) q[1];
sx q[1];
rz(0.57631667) q[1];
rz(2.6131931) q[2];
sx q[2];
rz(-0.64506275) q[2];
sx q[2];
rz(2.0730369) q[2];
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
