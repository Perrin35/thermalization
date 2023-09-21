OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7712819) q[0];
sx q[0];
rz(-0.9960649) q[0];
sx q[0];
rz(-0.87061849) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(-0.14970782) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71404845) q[0];
sx q[0];
rz(-2.2872426) q[0];
sx q[0];
rz(2.2358405) q[0];
x q[1];
rz(1.0295463) q[2];
sx q[2];
rz(-2.0143348) q[2];
sx q[2];
rz(0.29543791) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40587273) q[1];
sx q[1];
rz(-0.18438965) q[1];
sx q[1];
rz(-1.3609481) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1211795) q[3];
sx q[3];
rz(-1.4074416) q[3];
sx q[3];
rz(2.6640716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.25847882) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(2.1885833) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(0.63252226) q[0];
rz(0.44644341) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(-2.4893563) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3632293) q[0];
sx q[0];
rz(-0.031154545) q[0];
sx q[0];
rz(-0.40701436) q[0];
x q[1];
rz(2.2488238) q[2];
sx q[2];
rz(-0.36913482) q[2];
sx q[2];
rz(0.54112753) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3398509) q[1];
sx q[1];
rz(-2.7738214) q[1];
sx q[1];
rz(2.4778609) q[1];
x q[2];
rz(-2.2964301) q[3];
sx q[3];
rz(-2.3361428) q[3];
sx q[3];
rz(-1.0172539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5941045) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(1.9799505) q[2];
rz(-1.15796) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(-1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050425477) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(0.85025775) q[0];
rz(0.49750528) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(1.3495548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7499381) q[0];
sx q[0];
rz(-2.2313801) q[0];
sx q[0];
rz(-1.8750989) q[0];
rz(2.0850052) q[2];
sx q[2];
rz(-0.075767013) q[2];
sx q[2];
rz(2.5839992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9393443) q[1];
sx q[1];
rz(-1.0679686) q[1];
sx q[1];
rz(1.8400251) q[1];
rz(3.1317741) q[3];
sx q[3];
rz(-1.9968642) q[3];
sx q[3];
rz(0.94235086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0321908) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6501453) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(1.7310463) q[0];
rz(0.63181216) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(-3.1052123) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2934389) q[0];
sx q[0];
rz(-1.358195) q[0];
sx q[0];
rz(-0.92377499) q[0];
rz(0.58156275) q[2];
sx q[2];
rz(-1.8459324) q[2];
sx q[2];
rz(-0.046422596) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42602793) q[1];
sx q[1];
rz(-1.702311) q[1];
sx q[1];
rz(0.87042602) q[1];
rz(-pi) q[2];
rz(-2.9070204) q[3];
sx q[3];
rz(-1.4751225) q[3];
sx q[3];
rz(0.58805874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2146384) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(1.1882163) q[2];
rz(-0.67048091) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(-2.5454583) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(2.9751076) q[0];
rz(0.75603756) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(-0.23434815) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4011824) q[0];
sx q[0];
rz(-2.4966842) q[0];
sx q[0];
rz(0.58437225) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7469823) q[2];
sx q[2];
rz(-2.641045) q[2];
sx q[2];
rz(1.8793775) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.036451642) q[1];
sx q[1];
rz(-0.76117939) q[1];
sx q[1];
rz(-2.9245604) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0626416) q[3];
sx q[3];
rz(-1.2478932) q[3];
sx q[3];
rz(2.3372646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70871893) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(0.22053545) q[2];
rz(-2.7045414) q[3];
sx q[3];
rz(-2.1190937) q[3];
sx q[3];
rz(2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41546145) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(-0.81714001) q[0];
rz(2.5754886) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(1.9979427) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8457984) q[0];
sx q[0];
rz(-0.37308559) q[0];
sx q[0];
rz(-3.0315115) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32240378) q[2];
sx q[2];
rz(-0.69903261) q[2];
sx q[2];
rz(0.93271819) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3559349) q[1];
sx q[1];
rz(-2.6895084) q[1];
sx q[1];
rz(1.362734) q[1];
rz(-pi) q[2];
rz(0.46623047) q[3];
sx q[3];
rz(-0.56785184) q[3];
sx q[3];
rz(-3.1371491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.75366655) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(0.4894408) q[2];
rz(0.22805452) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(-0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.56931) q[0];
sx q[0];
rz(-2.4991878) q[0];
sx q[0];
rz(-1.2868767) q[0];
rz(-2.4781748) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(1.2333262) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39499261) q[0];
sx q[0];
rz(-1.9548423) q[0];
sx q[0];
rz(0.011944255) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2912824) q[2];
sx q[2];
rz(-1.0734953) q[2];
sx q[2];
rz(-0.26956272) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8334956) q[1];
sx q[1];
rz(-1.4596246) q[1];
sx q[1];
rz(1.584946) q[1];
x q[2];
rz(0.4864278) q[3];
sx q[3];
rz(-1.8842116) q[3];
sx q[3];
rz(0.13669554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.104091) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(2.1255778) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(-2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0163517) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(1.6960779) q[0];
rz(2.9267172) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(-1.8833556) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2381666) q[0];
sx q[0];
rz(-1.705372) q[0];
sx q[0];
rz(-2.0057136) q[0];
rz(-pi) q[1];
rz(-0.6446722) q[2];
sx q[2];
rz(-1.03627) q[2];
sx q[2];
rz(1.0169741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5453323) q[1];
sx q[1];
rz(-2.5858871) q[1];
sx q[1];
rz(-0.19821367) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97790896) q[3];
sx q[3];
rz(-2.3330354) q[3];
sx q[3];
rz(2.9582994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4593279) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(2.5788467) q[2];
rz(0.19966666) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(-1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257618) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(0.2510221) q[0];
rz(0.42731467) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(0.02773157) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4902892) q[0];
sx q[0];
rz(-0.9285183) q[0];
sx q[0];
rz(-0.62255967) q[0];
x q[1];
rz(0.8717732) q[2];
sx q[2];
rz(-1.7773526) q[2];
sx q[2];
rz(2.2733462) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8393644) q[1];
sx q[1];
rz(-1.6703969) q[1];
sx q[1];
rz(-2.539413) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6373709) q[3];
sx q[3];
rz(-1.596631) q[3];
sx q[3];
rz(-2.0500101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.015908265) q[2];
sx q[2];
rz(-0.97970825) q[2];
sx q[2];
rz(-1.998418) q[2];
rz(2.9987191) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(2.2843602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7509572) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(-0.45387682) q[0];
rz(-0.67165309) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(-2.8840816) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2902381) q[0];
sx q[0];
rz(-0.60831735) q[0];
sx q[0];
rz(-0.71382199) q[0];
x q[1];
rz(-2.4117878) q[2];
sx q[2];
rz(-2.4591755) q[2];
sx q[2];
rz(-0.44621106) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7911437) q[1];
sx q[1];
rz(-0.87462438) q[1];
sx q[1];
rz(1.5894366) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0912283) q[3];
sx q[3];
rz(-1.3551095) q[3];
sx q[3];
rz(-1.2558162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3283078) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(-0.60662398) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(-2.301208) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3474779) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(-2.9150302) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(-1.9526341) q[2];
sx q[2];
rz(-1.1321862) q[2];
sx q[2];
rz(1.95375) q[2];
rz(2.0541035) q[3];
sx q[3];
rz(-1.3369505) q[3];
sx q[3];
rz(-1.7575775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];