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
rz(-3.0814085) q[0];
sx q[0];
rz(-1.0589851) q[0];
sx q[0];
rz(-2.1314148) q[0];
rz(2.3330359) q[1];
sx q[1];
rz(-2.8454236) q[1];
sx q[1];
rz(-0.30997601) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3637562) q[0];
sx q[0];
rz(-1.68206) q[0];
sx q[0];
rz(-2.7194752) q[0];
rz(-pi) q[1];
rz(-1.593179) q[2];
sx q[2];
rz(-1.7760008) q[2];
sx q[2];
rz(0.6485282) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3821564) q[1];
sx q[1];
rz(-1.2431989) q[1];
sx q[1];
rz(-2.9345678) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3100389) q[3];
sx q[3];
rz(-1.420701) q[3];
sx q[3];
rz(0.41285601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4787204) q[2];
sx q[2];
rz(-2.8933849) q[2];
sx q[2];
rz(-0.09566801) q[2];
rz(3.0293448) q[3];
sx q[3];
rz(-0.87533689) q[3];
sx q[3];
rz(-1.4314502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9005168) q[0];
sx q[0];
rz(-2.2523585) q[0];
sx q[0];
rz(0.61479968) q[0];
rz(2.2531033) q[1];
sx q[1];
rz(-1.5526086) q[1];
sx q[1];
rz(0.49957553) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13321484) q[0];
sx q[0];
rz(-1.7724123) q[0];
sx q[0];
rz(-1.2382384) q[0];
rz(-pi) q[1];
rz(1.7287298) q[2];
sx q[2];
rz(-2.0896032) q[2];
sx q[2];
rz(0.11473303) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0608637) q[1];
sx q[1];
rz(-1.7434967) q[1];
sx q[1];
rz(0.91367803) q[1];
rz(2.8039475) q[3];
sx q[3];
rz(-1.0146552) q[3];
sx q[3];
rz(0.82586702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8932314) q[2];
sx q[2];
rz(-1.9042559) q[2];
sx q[2];
rz(0.020817967) q[2];
rz(1.9013532) q[3];
sx q[3];
rz(-1.6720684) q[3];
sx q[3];
rz(-1.9564691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6388539) q[0];
sx q[0];
rz(-0.26832142) q[0];
sx q[0];
rz(0.33367208) q[0];
rz(2.4036713) q[1];
sx q[1];
rz(-1.2260022) q[1];
sx q[1];
rz(0.12942448) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9343468) q[0];
sx q[0];
rz(-0.6967623) q[0];
sx q[0];
rz(1.3214701) q[0];
rz(-0.67141165) q[2];
sx q[2];
rz(-2.0575626) q[2];
sx q[2];
rz(1.1306131) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30363917) q[1];
sx q[1];
rz(-0.15175125) q[1];
sx q[1];
rz(0.93317731) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.09472) q[3];
sx q[3];
rz(-0.84730803) q[3];
sx q[3];
rz(1.0915966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0176598) q[2];
sx q[2];
rz(-1.61597) q[2];
sx q[2];
rz(-1.7714436) q[2];
rz(0.74639368) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(0.082775041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0081886) q[0];
sx q[0];
rz(-1.2097825) q[0];
sx q[0];
rz(-2.5764537) q[0];
rz(2.7124229) q[1];
sx q[1];
rz(-2.4950835) q[1];
sx q[1];
rz(0.82829222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1278616) q[0];
sx q[0];
rz(-0.7452226) q[0];
sx q[0];
rz(0.11786129) q[0];
x q[1];
rz(1.0617375) q[2];
sx q[2];
rz(-2.5974496) q[2];
sx q[2];
rz(0.5296112) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1431378) q[1];
sx q[1];
rz(-1.8363442) q[1];
sx q[1];
rz(2.6396855) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0542442) q[3];
sx q[3];
rz(-1.3540467) q[3];
sx q[3];
rz(-0.79424373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68507489) q[2];
sx q[2];
rz(-1.061331) q[2];
sx q[2];
rz(-1.431541) q[2];
rz(0.93959129) q[3];
sx q[3];
rz(-2.6288433) q[3];
sx q[3];
rz(-0.00069869839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13570304) q[0];
sx q[0];
rz(-2.8764184) q[0];
sx q[0];
rz(-1.3294719) q[0];
rz(1.1520518) q[1];
sx q[1];
rz(-1.0877129) q[1];
sx q[1];
rz(-0.00024814127) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6085998) q[0];
sx q[0];
rz(-1.8155351) q[0];
sx q[0];
rz(-0.18969638) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0543112) q[2];
sx q[2];
rz(-1.0382663) q[2];
sx q[2];
rz(-2.0244348) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6387512) q[1];
sx q[1];
rz(-1.0783245) q[1];
sx q[1];
rz(2.7472463) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4983921) q[3];
sx q[3];
rz(-2.3010217) q[3];
sx q[3];
rz(0.33354586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1382711) q[2];
sx q[2];
rz(-1.891581) q[2];
sx q[2];
rz(0.0058343466) q[2];
rz(-2.2516069) q[3];
sx q[3];
rz(-0.68166387) q[3];
sx q[3];
rz(2.0148923) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24790813) q[0];
sx q[0];
rz(-1.6981145) q[0];
sx q[0];
rz(-1.4877315) q[0];
rz(-0.030390175) q[1];
sx q[1];
rz(-1.1266212) q[1];
sx q[1];
rz(0.94246513) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93118942) q[0];
sx q[0];
rz(-1.7587587) q[0];
sx q[0];
rz(-0.24954777) q[0];
rz(-0.70955388) q[2];
sx q[2];
rz(-0.38564577) q[2];
sx q[2];
rz(0.31277671) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8928463) q[1];
sx q[1];
rz(-2.5809118) q[1];
sx q[1];
rz(1.0113202) q[1];
rz(-pi) q[2];
rz(1.0932572) q[3];
sx q[3];
rz(-2.6160588) q[3];
sx q[3];
rz(-2.9346187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5421062) q[2];
sx q[2];
rz(-2.3607871) q[2];
sx q[2];
rz(1.3055118) q[2];
rz(-2.5944338) q[3];
sx q[3];
rz(-1.9360417) q[3];
sx q[3];
rz(-2.3356596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18043537) q[0];
sx q[0];
rz(-2.1077709) q[0];
sx q[0];
rz(-0.10051522) q[0];
rz(0.48249498) q[1];
sx q[1];
rz(-1.9827739) q[1];
sx q[1];
rz(-0.85711342) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094497546) q[0];
sx q[0];
rz(-2.408354) q[0];
sx q[0];
rz(2.6493957) q[0];
rz(-3.0761511) q[2];
sx q[2];
rz(-2.4937154) q[2];
sx q[2];
rz(2.6104948) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7776612) q[1];
sx q[1];
rz(-2.6865733) q[1];
sx q[1];
rz(-2.6520686) q[1];
rz(-1.9590274) q[3];
sx q[3];
rz(-0.60743466) q[3];
sx q[3];
rz(-1.5608112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7414005) q[2];
sx q[2];
rz(-2.1705748) q[2];
sx q[2];
rz(2.8391489) q[2];
rz(-0.97964573) q[3];
sx q[3];
rz(-1.0023578) q[3];
sx q[3];
rz(-1.8625331) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.361146) q[0];
sx q[0];
rz(-0.13872153) q[0];
sx q[0];
rz(-0.030990344) q[0];
rz(-2.567645) q[1];
sx q[1];
rz(-1.4731044) q[1];
sx q[1];
rz(1.8642289) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81640667) q[0];
sx q[0];
rz(-1.2379097) q[0];
sx q[0];
rz(-1.216796) q[0];
rz(-1.1436815) q[2];
sx q[2];
rz(-0.53887212) q[2];
sx q[2];
rz(0.32921916) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9121418) q[1];
sx q[1];
rz(-2.0566642) q[1];
sx q[1];
rz(3.058601) q[1];
rz(-pi) q[2];
rz(-1.8559009) q[3];
sx q[3];
rz(-1.630097) q[3];
sx q[3];
rz(-1.3467195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0370827) q[2];
sx q[2];
rz(-3.0057378) q[2];
sx q[2];
rz(0.48745298) q[2];
rz(2.4449352) q[3];
sx q[3];
rz(-2.2127377) q[3];
sx q[3];
rz(0.063974403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0697407) q[0];
sx q[0];
rz(-0.47436473) q[0];
sx q[0];
rz(-0.35097861) q[0];
rz(0.17414302) q[1];
sx q[1];
rz(-1.5085647) q[1];
sx q[1];
rz(1.4016271) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021371776) q[0];
sx q[0];
rz(-1.8184156) q[0];
sx q[0];
rz(-1.6970474) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6675111) q[2];
sx q[2];
rz(-2.3326121) q[2];
sx q[2];
rz(-0.37341973) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.40608866) q[1];
sx q[1];
rz(-1.5998518) q[1];
sx q[1];
rz(1.3872223) q[1];
x q[2];
rz(-0.63302682) q[3];
sx q[3];
rz(-2.2928975) q[3];
sx q[3];
rz(-0.71114572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1104687) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(0.6366716) q[2];
rz(-1.1104256) q[3];
sx q[3];
rz(-1.597155) q[3];
sx q[3];
rz(-0.5184263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7745895) q[0];
sx q[0];
rz(-0.29357266) q[0];
sx q[0];
rz(2.0830182) q[0];
rz(-3.1029347) q[1];
sx q[1];
rz(-1.5267173) q[1];
sx q[1];
rz(-1.0640594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2702613) q[0];
sx q[0];
rz(-0.0055905213) q[0];
sx q[0];
rz(2.5020775) q[0];
x q[1];
rz(0.80697717) q[2];
sx q[2];
rz(-2.7514571) q[2];
sx q[2];
rz(-1.5240508) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1326931) q[1];
sx q[1];
rz(-1.510396) q[1];
sx q[1];
rz(1.6113043) q[1];
rz(0.070399447) q[3];
sx q[3];
rz(-2.681085) q[3];
sx q[3];
rz(-1.2706336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71558636) q[2];
sx q[2];
rz(-0.49804372) q[2];
sx q[2];
rz(-0.074020298) q[2];
rz(0.38153875) q[3];
sx q[3];
rz(-1.8266725) q[3];
sx q[3];
rz(-0.35416245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.030180177) q[0];
sx q[0];
rz(-1.8288061) q[0];
sx q[0];
rz(-0.70963138) q[0];
rz(-2.5297655) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(-0.62803531) q[2];
sx q[2];
rz(-2.5428094) q[2];
sx q[2];
rz(1.6369292) q[2];
rz(2.6388219) q[3];
sx q[3];
rz(-2.3270861) q[3];
sx q[3];
rz(2.6181639) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
