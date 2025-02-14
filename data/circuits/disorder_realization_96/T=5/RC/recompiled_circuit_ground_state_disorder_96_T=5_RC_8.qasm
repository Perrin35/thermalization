OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.191303) q[0];
sx q[0];
rz(-2.8714955) q[0];
sx q[0];
rz(2.252993) q[0];
rz(1.8285881) q[1];
sx q[1];
rz(4.7409952) q[1];
sx q[1];
rz(10.805605) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61533538) q[0];
sx q[0];
rz(-1.7563631) q[0];
sx q[0];
rz(0.19677563) q[0];
rz(2.3418597) q[2];
sx q[2];
rz(-1.1812484) q[2];
sx q[2];
rz(-2.3461297) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9819698) q[1];
sx q[1];
rz(-1.2352984) q[1];
sx q[1];
rz(1.5606176) q[1];
rz(-pi) q[2];
rz(1.2065776) q[3];
sx q[3];
rz(-1.3602339) q[3];
sx q[3];
rz(2.5488146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.31972739) q[2];
sx q[2];
rz(-1.3250019) q[2];
sx q[2];
rz(2.3925609) q[2];
rz(-2.9876515) q[3];
sx q[3];
rz(-1.0356244) q[3];
sx q[3];
rz(-0.099893959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.446949) q[0];
sx q[0];
rz(-1.6634989) q[0];
sx q[0];
rz(1.2167759) q[0];
rz(-1.0617537) q[1];
sx q[1];
rz(-0.80445015) q[1];
sx q[1];
rz(0.42713508) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1841284) q[0];
sx q[0];
rz(-1.2250568) q[0];
sx q[0];
rz(-0.60681245) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2522302) q[2];
sx q[2];
rz(-0.67020352) q[2];
sx q[2];
rz(2.2549652) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1698902) q[1];
sx q[1];
rz(-2.8198543) q[1];
sx q[1];
rz(1.1119026) q[1];
rz(-1.7747545) q[3];
sx q[3];
rz(-1.5195623) q[3];
sx q[3];
rz(0.27328396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.079166807) q[2];
sx q[2];
rz(-1.4126974) q[2];
sx q[2];
rz(1.742935) q[2];
rz(2.9178197) q[3];
sx q[3];
rz(-2.2327435) q[3];
sx q[3];
rz(1.5435425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.021521213) q[0];
sx q[0];
rz(-1.3503617) q[0];
sx q[0];
rz(-0.045510005) q[0];
rz(1.1687763) q[1];
sx q[1];
rz(-1.2380506) q[1];
sx q[1];
rz(-2.5659836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63059124) q[0];
sx q[0];
rz(-1.5596034) q[0];
sx q[0];
rz(-0.63909792) q[0];
rz(-2.2314638) q[2];
sx q[2];
rz(-1.5843387) q[2];
sx q[2];
rz(0.34156583) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8519046) q[1];
sx q[1];
rz(-0.59483268) q[1];
sx q[1];
rz(-1.1762876) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6034254) q[3];
sx q[3];
rz(-0.65308076) q[3];
sx q[3];
rz(0.075270502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1316954) q[2];
sx q[2];
rz(-0.74247777) q[2];
sx q[2];
rz(2.1843145) q[2];
rz(-0.57473985) q[3];
sx q[3];
rz(-1.6114019) q[3];
sx q[3];
rz(1.8360229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9855373) q[0];
sx q[0];
rz(-0.95273459) q[0];
sx q[0];
rz(1.2611457) q[0];
rz(3.0592697) q[1];
sx q[1];
rz(-1.0876834) q[1];
sx q[1];
rz(-1.3538768) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9808423) q[0];
sx q[0];
rz(-0.35355648) q[0];
sx q[0];
rz(2.4048231) q[0];
x q[1];
rz(0.48575966) q[2];
sx q[2];
rz(-1.4076715) q[2];
sx q[2];
rz(2.8201134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8396047) q[1];
sx q[1];
rz(-2.0652373) q[1];
sx q[1];
rz(-1.6322664) q[1];
rz(-0.98815496) q[3];
sx q[3];
rz(-2.372962) q[3];
sx q[3];
rz(-2.5386794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42133078) q[2];
sx q[2];
rz(-2.223184) q[2];
sx q[2];
rz(-1.8850231) q[2];
rz(2.4300857) q[3];
sx q[3];
rz(-1.3308176) q[3];
sx q[3];
rz(2.8594657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8307777) q[0];
sx q[0];
rz(-2.2956235) q[0];
sx q[0];
rz(-0.34314439) q[0];
rz(0.061773069) q[1];
sx q[1];
rz(-0.97007483) q[1];
sx q[1];
rz(-1.7247346) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23487906) q[0];
sx q[0];
rz(-0.32787927) q[0];
sx q[0];
rz(2.1239141) q[0];
rz(-pi) q[1];
rz(-0.38133905) q[2];
sx q[2];
rz(-0.73753192) q[2];
sx q[2];
rz(-2.2708837) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0793905) q[1];
sx q[1];
rz(-0.50143999) q[1];
sx q[1];
rz(-1.8556684) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9067578) q[3];
sx q[3];
rz(-1.9855301) q[3];
sx q[3];
rz(2.107634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5477649) q[2];
sx q[2];
rz(-1.6739028) q[2];
sx q[2];
rz(-1.0859547) q[2];
rz(0.91935277) q[3];
sx q[3];
rz(-1.7707526) q[3];
sx q[3];
rz(-0.61387387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3980961) q[0];
sx q[0];
rz(-1.9323823) q[0];
sx q[0];
rz(-0.79291517) q[0];
rz(-2.4010557) q[1];
sx q[1];
rz(-0.99895993) q[1];
sx q[1];
rz(-2.3103796) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8342469) q[0];
sx q[0];
rz(-1.9458658) q[0];
sx q[0];
rz(0.22386472) q[0];
rz(-pi) q[1];
rz(1.8478644) q[2];
sx q[2];
rz(-1.5779621) q[2];
sx q[2];
rz(2.9286419) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.97396353) q[1];
sx q[1];
rz(-2.2669753) q[1];
sx q[1];
rz(-2.7401398) q[1];
rz(1.6115723) q[3];
sx q[3];
rz(-0.47104657) q[3];
sx q[3];
rz(-2.149793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.071216019) q[2];
sx q[2];
rz(-1.8698317) q[2];
sx q[2];
rz(1.1191204) q[2];
rz(-1.8708771) q[3];
sx q[3];
rz(-1.1071353) q[3];
sx q[3];
rz(-1.7601815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0041466) q[0];
sx q[0];
rz(-2.9747712) q[0];
sx q[0];
rz(1.5420089) q[0];
rz(-2.1265325) q[1];
sx q[1];
rz(-1.5856182) q[1];
sx q[1];
rz(0.17280811) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27903444) q[0];
sx q[0];
rz(-2.1365215) q[0];
sx q[0];
rz(0.59086694) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68384275) q[2];
sx q[2];
rz(-1.8617861) q[2];
sx q[2];
rz(-3.0802205) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0902675) q[1];
sx q[1];
rz(-0.50053144) q[1];
sx q[1];
rz(-2.1085018) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0181581) q[3];
sx q[3];
rz(-2.2595836) q[3];
sx q[3];
rz(1.3017201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66199866) q[2];
sx q[2];
rz(-0.84344232) q[2];
sx q[2];
rz(2.1638828) q[2];
rz(0.57502037) q[3];
sx q[3];
rz(-1.2049371) q[3];
sx q[3];
rz(1.9112816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32790312) q[0];
sx q[0];
rz(-0.29569018) q[0];
sx q[0];
rz(3.1258702) q[0];
rz(-1.9533336) q[1];
sx q[1];
rz(-2.5091722) q[1];
sx q[1];
rz(-1.9421008) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2000113) q[0];
sx q[0];
rz(-2.0111548) q[0];
sx q[0];
rz(2.8811127) q[0];
rz(-2.0616777) q[2];
sx q[2];
rz(-1.5245617) q[2];
sx q[2];
rz(2.2904928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6024692) q[1];
sx q[1];
rz(-2.0401461) q[1];
sx q[1];
rz(-2.5825809) q[1];
rz(-pi) q[2];
rz(0.30808456) q[3];
sx q[3];
rz(-1.020806) q[3];
sx q[3];
rz(0.96984449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1477995) q[2];
sx q[2];
rz(-1.2387929) q[2];
sx q[2];
rz(-1.3851059) q[2];
rz(-2.5177453) q[3];
sx q[3];
rz(-1.8892989) q[3];
sx q[3];
rz(2.0417716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95595908) q[0];
sx q[0];
rz(-0.82056844) q[0];
sx q[0];
rz(-0.69325915) q[0];
rz(-2.8765826) q[1];
sx q[1];
rz(-0.82844228) q[1];
sx q[1];
rz(-1.5230491) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1177445) q[0];
sx q[0];
rz(-1.5818412) q[0];
sx q[0];
rz(-2.3680229) q[0];
x q[1];
rz(-0.74769865) q[2];
sx q[2];
rz(-0.77538632) q[2];
sx q[2];
rz(0.86624399) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56227389) q[1];
sx q[1];
rz(-2.8011311) q[1];
sx q[1];
rz(0.071694362) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6783488) q[3];
sx q[3];
rz(-2.5578024) q[3];
sx q[3];
rz(-1.8369758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.30274621) q[2];
sx q[2];
rz(-1.6281717) q[2];
sx q[2];
rz(-1.6839074) q[2];
rz(-0.95528209) q[3];
sx q[3];
rz(-2.6421319) q[3];
sx q[3];
rz(-2.3063851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38462287) q[0];
sx q[0];
rz(-2.5769825) q[0];
sx q[0];
rz(-2.6589174) q[0];
rz(2.2118498) q[1];
sx q[1];
rz(-2.1371806) q[1];
sx q[1];
rz(-2.7453056) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1866465) q[0];
sx q[0];
rz(-2.5781693) q[0];
sx q[0];
rz(-2.4231829) q[0];
rz(0.17896374) q[2];
sx q[2];
rz(-1.4216627) q[2];
sx q[2];
rz(-0.75527945) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0490108) q[1];
sx q[1];
rz(-1.3900458) q[1];
sx q[1];
rz(-2.3011219) q[1];
rz(1.2108874) q[3];
sx q[3];
rz(-1.1904089) q[3];
sx q[3];
rz(1.4402657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3489909) q[2];
sx q[2];
rz(-1.4395809) q[2];
sx q[2];
rz(-2.8144042) q[2];
rz(-0.78091019) q[3];
sx q[3];
rz(-1.1612929) q[3];
sx q[3];
rz(1.4769295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0384211) q[0];
sx q[0];
rz(-2.1506943) q[0];
sx q[0];
rz(-2.8839169) q[0];
rz(0.36956638) q[1];
sx q[1];
rz(-2.259544) q[1];
sx q[1];
rz(2.4824711) q[1];
rz(0.054389537) q[2];
sx q[2];
rz(-2.4475606) q[2];
sx q[2];
rz(-2.9771752) q[2];
rz(-2.7575708) q[3];
sx q[3];
rz(-0.99746084) q[3];
sx q[3];
rz(1.8538047) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
