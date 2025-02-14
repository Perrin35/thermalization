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
rz(1.3311812) q[0];
sx q[0];
rz(-1.5067195) q[0];
sx q[0];
rz(0.17568406) q[0];
rz(1.0653347) q[1];
sx q[1];
rz(4.4720286) q[1];
sx q[1];
rz(8.7695697) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5979008) q[0];
sx q[0];
rz(-0.59698623) q[0];
sx q[0];
rz(-2.1587121) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5152446) q[2];
sx q[2];
rz(-1.6695938) q[2];
sx q[2];
rz(2.9227119) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0111078) q[1];
sx q[1];
rz(-1.4629416) q[1];
sx q[1];
rz(1.4977096) q[1];
x q[2];
rz(0.96027042) q[3];
sx q[3];
rz(-1.8132079) q[3];
sx q[3];
rz(-1.2038972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6931927) q[2];
sx q[2];
rz(-3.0312263) q[2];
sx q[2];
rz(2.8006862) q[2];
rz(2.36813) q[3];
sx q[3];
rz(-1.0186467) q[3];
sx q[3];
rz(0.10507467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.676749) q[0];
sx q[0];
rz(-0.28443795) q[0];
sx q[0];
rz(-2.8819717) q[0];
rz(0.78850857) q[1];
sx q[1];
rz(-2.2854243) q[1];
sx q[1];
rz(1.2451008) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21331025) q[0];
sx q[0];
rz(-1.4699555) q[0];
sx q[0];
rz(-2.1891602) q[0];
rz(-pi) q[1];
rz(0.2971895) q[2];
sx q[2];
rz(-2.8126394) q[2];
sx q[2];
rz(2.7847814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40971995) q[1];
sx q[1];
rz(-1.9202246) q[1];
sx q[1];
rz(-0.81373377) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29446843) q[3];
sx q[3];
rz(-0.52669062) q[3];
sx q[3];
rz(1.1931488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9022687) q[2];
sx q[2];
rz(-2.6332492) q[2];
sx q[2];
rz(-0.90947214) q[2];
rz(-0.66696143) q[3];
sx q[3];
rz(-1.704155) q[3];
sx q[3];
rz(3.0392569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42720902) q[0];
sx q[0];
rz(-2.7294071) q[0];
sx q[0];
rz(1.7538196) q[0];
rz(-0.99675238) q[1];
sx q[1];
rz(-2.4501188) q[1];
sx q[1];
rz(-2.6549285) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61302859) q[0];
sx q[0];
rz(-1.4183622) q[0];
sx q[0];
rz(-1.9804762) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85993729) q[2];
sx q[2];
rz(-0.87274466) q[2];
sx q[2];
rz(2.6280478) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4506766) q[1];
sx q[1];
rz(-1.449739) q[1];
sx q[1];
rz(0.52074213) q[1];
rz(2.9938042) q[3];
sx q[3];
rz(-0.94919616) q[3];
sx q[3];
rz(-3.0382267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.77728689) q[2];
sx q[2];
rz(-1.8824258) q[2];
sx q[2];
rz(2.6510009) q[2];
rz(-2.3250438) q[3];
sx q[3];
rz(-2.3994763) q[3];
sx q[3];
rz(-2.1013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33646026) q[0];
sx q[0];
rz(-1.8959624) q[0];
sx q[0];
rz(-2.1601987) q[0];
rz(-2.7384752) q[1];
sx q[1];
rz(-0.50794452) q[1];
sx q[1];
rz(0.53781646) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92141672) q[0];
sx q[0];
rz(-0.28982899) q[0];
sx q[0];
rz(-1.0224065) q[0];
rz(-0.85015202) q[2];
sx q[2];
rz(-1.6714897) q[2];
sx q[2];
rz(0.062391524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.06925) q[1];
sx q[1];
rz(-1.1872655) q[1];
sx q[1];
rz(1.0479397) q[1];
rz(-pi) q[2];
rz(-1.1868401) q[3];
sx q[3];
rz(-0.71184413) q[3];
sx q[3];
rz(-1.2983589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2526523) q[2];
sx q[2];
rz(-2.6344968) q[2];
sx q[2];
rz(-1.0520774) q[2];
rz(2.5893411) q[3];
sx q[3];
rz(-1.4058607) q[3];
sx q[3];
rz(1.9210057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1110765) q[0];
sx q[0];
rz(-1.3330326) q[0];
sx q[0];
rz(-0.037574969) q[0];
rz(0.75282085) q[1];
sx q[1];
rz(-1.5767187) q[1];
sx q[1];
rz(-3.0573696) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0287474) q[0];
sx q[0];
rz(-2.0984834) q[0];
sx q[0];
rz(-0.28808388) q[0];
rz(-1.8657612) q[2];
sx q[2];
rz(-1.2216481) q[2];
sx q[2];
rz(-2.5087506) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4956717) q[1];
sx q[1];
rz(-2.715647) q[1];
sx q[1];
rz(-1.6712214) q[1];
rz(-1.7315167) q[3];
sx q[3];
rz(-1.525021) q[3];
sx q[3];
rz(1.20884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5273744) q[2];
sx q[2];
rz(-2.1671593) q[2];
sx q[2];
rz(2.7091889) q[2];
rz(-1.6480986) q[3];
sx q[3];
rz(-1.3727539) q[3];
sx q[3];
rz(-2.4026332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37345165) q[0];
sx q[0];
rz(-1.6141163) q[0];
sx q[0];
rz(-0.76171184) q[0];
rz(-2.1181882) q[1];
sx q[1];
rz(-0.49097148) q[1];
sx q[1];
rz(-0.34058079) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7135431) q[0];
sx q[0];
rz(-0.78589702) q[0];
sx q[0];
rz(-3.0160496) q[0];
rz(1.7271785) q[2];
sx q[2];
rz(-2.1264653) q[2];
sx q[2];
rz(-1.2909691) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2557166) q[1];
sx q[1];
rz(-1.2310561) q[1];
sx q[1];
rz(2.8541616) q[1];
x q[2];
rz(-1.4177001) q[3];
sx q[3];
rz(-2.3153437) q[3];
sx q[3];
rz(1.5797256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9684101) q[2];
sx q[2];
rz(-2.0081655) q[2];
sx q[2];
rz(-0.75562149) q[2];
rz(-0.38718811) q[3];
sx q[3];
rz(-1.7500992) q[3];
sx q[3];
rz(-2.3473327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25283915) q[0];
sx q[0];
rz(-0.08789739) q[0];
sx q[0];
rz(-1.2095691) q[0];
rz(-1.5225438) q[1];
sx q[1];
rz(-2.3698898) q[1];
sx q[1];
rz(-1.3873772) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1529941) q[0];
sx q[0];
rz(-1.261088) q[0];
sx q[0];
rz(-1.9266846) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0389236) q[2];
sx q[2];
rz(-2.515677) q[2];
sx q[2];
rz(-2.4374838) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.61664192) q[1];
sx q[1];
rz(-1.157028) q[1];
sx q[1];
rz(-2.3828808) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8903743) q[3];
sx q[3];
rz(-1.1261903) q[3];
sx q[3];
rz(-0.79892677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1986971) q[2];
sx q[2];
rz(-0.2936475) q[2];
sx q[2];
rz(2.7315308) q[2];
rz(-0.53747082) q[3];
sx q[3];
rz(-1.042807) q[3];
sx q[3];
rz(1.0679831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10802565) q[0];
sx q[0];
rz(-1.2355868) q[0];
sx q[0];
rz(-2.1349452) q[0];
rz(1.135723) q[1];
sx q[1];
rz(-0.4907116) q[1];
sx q[1];
rz(-2.8699285) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0245117) q[0];
sx q[0];
rz(-0.70444626) q[0];
sx q[0];
rz(-0.68450494) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3799963) q[2];
sx q[2];
rz(-0.83405364) q[2];
sx q[2];
rz(2.746838) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2073839) q[1];
sx q[1];
rz(-1.0107688) q[1];
sx q[1];
rz(-2.4337342) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8992796) q[3];
sx q[3];
rz(-0.75868443) q[3];
sx q[3];
rz(-3.098947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9378308) q[2];
sx q[2];
rz(-0.93272847) q[2];
sx q[2];
rz(2.0358098) q[2];
rz(-1.743099) q[3];
sx q[3];
rz(-0.98413697) q[3];
sx q[3];
rz(-2.2891146) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49649134) q[0];
sx q[0];
rz(-2.190525) q[0];
sx q[0];
rz(-2.7880461) q[0];
rz(-1.722466) q[1];
sx q[1];
rz(-1.5357693) q[1];
sx q[1];
rz(0.062573418) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.94189) q[0];
sx q[0];
rz(-2.1158005) q[0];
sx q[0];
rz(0.99129321) q[0];
rz(0.57762164) q[2];
sx q[2];
rz(-2.2934539) q[2];
sx q[2];
rz(1.3179145) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.313301) q[1];
sx q[1];
rz(-0.57692674) q[1];
sx q[1];
rz(-2.5531656) q[1];
rz(-pi) q[2];
rz(-3.022404) q[3];
sx q[3];
rz(-2.475563) q[3];
sx q[3];
rz(-2.5882849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.320437) q[2];
sx q[2];
rz(-1.374036) q[2];
sx q[2];
rz(2.721411) q[2];
rz(-0.34475103) q[3];
sx q[3];
rz(-0.77778608) q[3];
sx q[3];
rz(-2.2043118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6688113) q[0];
sx q[0];
rz(-2.9720699) q[0];
sx q[0];
rz(1.2813168) q[0];
rz(-1.5538813) q[1];
sx q[1];
rz(-2.3322767) q[1];
sx q[1];
rz(-0.70714998) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90072322) q[0];
sx q[0];
rz(-1.5587731) q[0];
sx q[0];
rz(0.94038402) q[0];
x q[1];
rz(-2.1538582) q[2];
sx q[2];
rz(-1.9473922) q[2];
sx q[2];
rz(2.0104806) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8445804) q[1];
sx q[1];
rz(-1.2606916) q[1];
sx q[1];
rz(-1.8052141) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1955802) q[3];
sx q[3];
rz(-2.1008228) q[3];
sx q[3];
rz(1.8998787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3215434) q[2];
sx q[2];
rz(-1.1300491) q[2];
sx q[2];
rz(0.496544) q[2];
rz(-1.8494362) q[3];
sx q[3];
rz(-0.29785952) q[3];
sx q[3];
rz(-0.37798247) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9847066) q[0];
sx q[0];
rz(-0.26837415) q[0];
sx q[0];
rz(-1.0347086) q[0];
rz(2.5490419) q[1];
sx q[1];
rz(-0.68065803) q[1];
sx q[1];
rz(-1.5417644) q[1];
rz(-1.401027) q[2];
sx q[2];
rz(-0.23525722) q[2];
sx q[2];
rz(2.7973882) q[2];
rz(1.7962528) q[3];
sx q[3];
rz(-0.97889789) q[3];
sx q[3];
rz(-1.2759664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
