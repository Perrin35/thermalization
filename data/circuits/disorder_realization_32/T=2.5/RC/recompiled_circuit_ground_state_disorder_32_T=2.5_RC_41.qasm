OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6736421) q[0];
sx q[0];
rz(-1.8704432) q[0];
sx q[0];
rz(-0.76572642) q[0];
rz(2.7453121) q[1];
sx q[1];
rz(-3.0488465) q[1];
sx q[1];
rz(-0.5527817) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3512289) q[0];
sx q[0];
rz(-1.1811331) q[0];
sx q[0];
rz(-2.0084951) q[0];
rz(-0.88199204) q[2];
sx q[2];
rz(-2.246736) q[2];
sx q[2];
rz(-1.9641528) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.89239425) q[1];
sx q[1];
rz(-0.65649872) q[1];
sx q[1];
rz(1.0358443) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90867282) q[3];
sx q[3];
rz(-1.931802) q[3];
sx q[3];
rz(0.26546016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3679686) q[2];
sx q[2];
rz(-1.6398733) q[2];
sx q[2];
rz(0.088851301) q[2];
rz(-0.39696524) q[3];
sx q[3];
rz(-2.1019955) q[3];
sx q[3];
rz(1.6840434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6615768) q[0];
sx q[0];
rz(-2.6341697) q[0];
sx q[0];
rz(0.85451025) q[0];
rz(0.62141934) q[1];
sx q[1];
rz(-0.3365376) q[1];
sx q[1];
rz(-1.052676) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5993505) q[0];
sx q[0];
rz(-2.9711038) q[0];
sx q[0];
rz(0.64298274) q[0];
rz(-pi) q[1];
rz(-2.3143155) q[2];
sx q[2];
rz(-0.57663871) q[2];
sx q[2];
rz(2.1183543) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4667644) q[1];
sx q[1];
rz(-0.75777868) q[1];
sx q[1];
rz(-1.2762875) q[1];
x q[2];
rz(-1.5731811) q[3];
sx q[3];
rz(-2.5376476) q[3];
sx q[3];
rz(-1.2177474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66775995) q[2];
sx q[2];
rz(-2.2809873) q[2];
sx q[2];
rz(-0.94998002) q[2];
rz(-1.0085603) q[3];
sx q[3];
rz(-2.5086094) q[3];
sx q[3];
rz(2.8620201) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94473332) q[0];
sx q[0];
rz(-1.2161398) q[0];
sx q[0];
rz(1.3470294) q[0];
rz(-2.0836209) q[1];
sx q[1];
rz(-0.24669138) q[1];
sx q[1];
rz(-0.11014858) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5025217) q[0];
sx q[0];
rz(-2.5701231) q[0];
sx q[0];
rz(-2.1854464) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0579849) q[2];
sx q[2];
rz(-1.4232529) q[2];
sx q[2];
rz(1.7718441) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.55787841) q[1];
sx q[1];
rz(-1.6772575) q[1];
sx q[1];
rz(-0.6006247) q[1];
rz(-pi) q[2];
rz(2.3055656) q[3];
sx q[3];
rz(-0.6202094) q[3];
sx q[3];
rz(-0.42441747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9614253) q[2];
sx q[2];
rz(-1.578293) q[2];
sx q[2];
rz(3.0999198) q[2];
rz(-2.9336119) q[3];
sx q[3];
rz(-0.22170034) q[3];
sx q[3];
rz(2.8878133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5623077) q[0];
sx q[0];
rz(-1.8203745) q[0];
sx q[0];
rz(-0.98096171) q[0];
rz(2.7085069) q[1];
sx q[1];
rz(-1.667645) q[1];
sx q[1];
rz(-1.6281737) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53674066) q[0];
sx q[0];
rz(-1.3910595) q[0];
sx q[0];
rz(0.40412235) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5879405) q[2];
sx q[2];
rz(-0.671467) q[2];
sx q[2];
rz(1.4070828) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.95483741) q[1];
sx q[1];
rz(-1.4680956) q[1];
sx q[1];
rz(-2.0689194) q[1];
rz(-pi) q[2];
rz(0.93928503) q[3];
sx q[3];
rz(-1.8016677) q[3];
sx q[3];
rz(-2.758647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.607434) q[2];
sx q[2];
rz(-0.51924339) q[2];
sx q[2];
rz(-2.329211) q[2];
rz(1.8937998) q[3];
sx q[3];
rz(-1.8915853) q[3];
sx q[3];
rz(1.8947424) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8060551) q[0];
sx q[0];
rz(-2.1195109) q[0];
sx q[0];
rz(1.442765) q[0];
rz(0.45817786) q[1];
sx q[1];
rz(-1.5135601) q[1];
sx q[1];
rz(2.3853669) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6578015) q[0];
sx q[0];
rz(-1.7258125) q[0];
sx q[0];
rz(-2.20138) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5014072) q[2];
sx q[2];
rz(-2.019503) q[2];
sx q[2];
rz(1.9404836) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9971265) q[1];
sx q[1];
rz(-1.6580075) q[1];
sx q[1];
rz(0.13551305) q[1];
x q[2];
rz(-1.4706572) q[3];
sx q[3];
rz(-2.7960145) q[3];
sx q[3];
rz(2.341193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4592287) q[2];
sx q[2];
rz(-1.3882779) q[2];
sx q[2];
rz(2.1759822) q[2];
rz(-0.56973488) q[3];
sx q[3];
rz(-2.0163048) q[3];
sx q[3];
rz(-1.6857356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7392893) q[0];
sx q[0];
rz(-1.1981413) q[0];
sx q[0];
rz(-1.385561) q[0];
rz(-2.3784474) q[1];
sx q[1];
rz(-1.4475854) q[1];
sx q[1];
rz(-0.10173434) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83086936) q[0];
sx q[0];
rz(-2.1889157) q[0];
sx q[0];
rz(-3.028426) q[0];
rz(-2.7587682) q[2];
sx q[2];
rz(-1.5271795) q[2];
sx q[2];
rz(-0.62784615) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26120034) q[1];
sx q[1];
rz(-1.8414547) q[1];
sx q[1];
rz(-0.063505725) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2564628) q[3];
sx q[3];
rz(-2.1477774) q[3];
sx q[3];
rz(0.5761522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.10852854) q[2];
sx q[2];
rz(-1.6988924) q[2];
sx q[2];
rz(-0.70890439) q[2];
rz(-0.76062834) q[3];
sx q[3];
rz(-1.9899188) q[3];
sx q[3];
rz(0.93543783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075542299) q[0];
sx q[0];
rz(-1.9882555) q[0];
sx q[0];
rz(-0.72108889) q[0];
rz(-0.9779633) q[1];
sx q[1];
rz(-2.5826192) q[1];
sx q[1];
rz(1.579938) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52873351) q[0];
sx q[0];
rz(-1.4943143) q[0];
sx q[0];
rz(1.4974817) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97550895) q[2];
sx q[2];
rz(-1.3689318) q[2];
sx q[2];
rz(1.6109811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.60266337) q[1];
sx q[1];
rz(-2.437535) q[1];
sx q[1];
rz(2.5119147) q[1];
rz(2.2411602) q[3];
sx q[3];
rz(-1.2868902) q[3];
sx q[3];
rz(-2.8395124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8568023) q[2];
sx q[2];
rz(-2.6140116) q[2];
sx q[2];
rz(-1.52012) q[2];
rz(0.4367477) q[3];
sx q[3];
rz(-0.89478651) q[3];
sx q[3];
rz(2.6020218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.082551) q[0];
sx q[0];
rz(-2.3024547) q[0];
sx q[0];
rz(-0.82373291) q[0];
rz(-2.0027022) q[1];
sx q[1];
rz(-1.7887807) q[1];
sx q[1];
rz(-1.9653856) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3977073) q[0];
sx q[0];
rz(-1.526014) q[0];
sx q[0];
rz(2.0732047) q[0];
rz(-pi) q[1];
rz(1.9763821) q[2];
sx q[2];
rz(-2.1422804) q[2];
sx q[2];
rz(1.5290844) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.35299435) q[1];
sx q[1];
rz(-1.561015) q[1];
sx q[1];
rz(-0.23810975) q[1];
rz(-pi) q[2];
rz(1.6105936) q[3];
sx q[3];
rz(-2.1615513) q[3];
sx q[3];
rz(-1.9692957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99522432) q[2];
sx q[2];
rz(-0.58688846) q[2];
sx q[2];
rz(-0.5874908) q[2];
rz(-0.38582173) q[3];
sx q[3];
rz(-2.2930175) q[3];
sx q[3];
rz(2.2390656) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1189482) q[0];
sx q[0];
rz(-2.1653439) q[0];
sx q[0];
rz(2.5318085) q[0];
rz(1.7976286) q[1];
sx q[1];
rz(-1.983843) q[1];
sx q[1];
rz(-0.15300289) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1399779) q[0];
sx q[0];
rz(-1.5904332) q[0];
sx q[0];
rz(2.6756712) q[0];
rz(1.5605037) q[2];
sx q[2];
rz(-0.77559815) q[2];
sx q[2];
rz(-2.6960109) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9003332) q[1];
sx q[1];
rz(-1.5521084) q[1];
sx q[1];
rz(0.92872031) q[1];
rz(2.1243662) q[3];
sx q[3];
rz(-2.3690936) q[3];
sx q[3];
rz(-1.8541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0387705) q[2];
sx q[2];
rz(-1.5615347) q[2];
sx q[2];
rz(2.5193396) q[2];
rz(-1.9410939) q[3];
sx q[3];
rz(-1.7222722) q[3];
sx q[3];
rz(2.5672369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822561) q[0];
sx q[0];
rz(-0.065956235) q[0];
sx q[0];
rz(2.7863853) q[0];
rz(1.1025053) q[1];
sx q[1];
rz(-0.016977221) q[1];
sx q[1];
rz(0.65974081) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6924877) q[0];
sx q[0];
rz(-0.67506719) q[0];
sx q[0];
rz(2.8646742) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8546472) q[2];
sx q[2];
rz(-2.7128993) q[2];
sx q[2];
rz(2.5108166) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8568032) q[1];
sx q[1];
rz(-1.8823913) q[1];
sx q[1];
rz(-1.9143399) q[1];
rz(-1.6361314) q[3];
sx q[3];
rz(-2.4628277) q[3];
sx q[3];
rz(-2.9836054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.79992646) q[2];
sx q[2];
rz(-0.081192668) q[2];
sx q[2];
rz(-0.1989092) q[2];
rz(0.79814664) q[3];
sx q[3];
rz(-1.9559559) q[3];
sx q[3];
rz(1.2822436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2220919) q[0];
sx q[0];
rz(-1.6035447) q[0];
sx q[0];
rz(-1.0731687) q[0];
rz(0.78492289) q[1];
sx q[1];
rz(-0.85826086) q[1];
sx q[1];
rz(0.8716743) q[1];
rz(1.2047819) q[2];
sx q[2];
rz(-0.13199619) q[2];
sx q[2];
rz(0.97313626) q[2];
rz(-2.4724805) q[3];
sx q[3];
rz(-1.787686) q[3];
sx q[3];
rz(-0.91989582) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
