OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9935432) q[0];
sx q[0];
rz(-2.7033959) q[0];
sx q[0];
rz(2.4105657) q[0];
rz(-2.464715) q[1];
sx q[1];
rz(-0.69898611) q[1];
sx q[1];
rz(0.55404034) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1822788) q[0];
sx q[0];
rz(-2.4327299) q[0];
sx q[0];
rz(-0.76226632) q[0];
rz(1.117716) q[2];
sx q[2];
rz(-0.52342285) q[2];
sx q[2];
rz(-0.40016178) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4198534) q[1];
sx q[1];
rz(-1.6045147) q[1];
sx q[1];
rz(0.88326247) q[1];
x q[2];
rz(-2.2525134) q[3];
sx q[3];
rz(-1.9535716) q[3];
sx q[3];
rz(2.0815522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76250166) q[2];
sx q[2];
rz(-1.5753626) q[2];
sx q[2];
rz(0.16538922) q[2];
rz(-0.23073828) q[3];
sx q[3];
rz(-2.8985891) q[3];
sx q[3];
rz(2.2802584) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47925258) q[0];
sx q[0];
rz(-0.32821822) q[0];
sx q[0];
rz(-1.7829371) q[0];
rz(-1.6183629) q[1];
sx q[1];
rz(-2.3871469) q[1];
sx q[1];
rz(-0.84017909) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019842783) q[0];
sx q[0];
rz(-0.54980924) q[0];
sx q[0];
rz(-1.0417095) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7457623) q[2];
sx q[2];
rz(-2.2326222) q[2];
sx q[2];
rz(-1.0368376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7666553) q[1];
sx q[1];
rz(-1.8647653) q[1];
sx q[1];
rz(2.4092182) q[1];
rz(-pi) q[2];
rz(-2.9884913) q[3];
sx q[3];
rz(-2.370122) q[3];
sx q[3];
rz(0.16530802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3885865) q[2];
sx q[2];
rz(-1.2331839) q[2];
sx q[2];
rz(-2.2071655) q[2];
rz(-1.8560483) q[3];
sx q[3];
rz(-2.3617187) q[3];
sx q[3];
rz(2.2540895) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0282106) q[0];
sx q[0];
rz(-2.1206355) q[0];
sx q[0];
rz(1.5519979) q[0];
rz(1.3681083) q[1];
sx q[1];
rz(-1.869447) q[1];
sx q[1];
rz(2.0210463) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5277953) q[0];
sx q[0];
rz(-0.58301914) q[0];
sx q[0];
rz(2.1979419) q[0];
x q[1];
rz(2.8875966) q[2];
sx q[2];
rz(-1.996737) q[2];
sx q[2];
rz(-2.6563702) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9160794) q[1];
sx q[1];
rz(-2.597993) q[1];
sx q[1];
rz(-0.17982843) q[1];
rz(-pi) q[2];
rz(0.12579671) q[3];
sx q[3];
rz(-1.9748944) q[3];
sx q[3];
rz(0.27751291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0418732) q[2];
sx q[2];
rz(-2.9283044) q[2];
sx q[2];
rz(-0.29207692) q[2];
rz(1.2415576) q[3];
sx q[3];
rz(-1.440666) q[3];
sx q[3];
rz(0.45274538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9860155) q[0];
sx q[0];
rz(-0.41322511) q[0];
sx q[0];
rz(0.43169942) q[0];
rz(2.9875634) q[1];
sx q[1];
rz(-1.5700211) q[1];
sx q[1];
rz(-1.0607176) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78846473) q[0];
sx q[0];
rz(-1.2305088) q[0];
sx q[0];
rz(-2.8022604) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3260822) q[2];
sx q[2];
rz(-1.3977461) q[2];
sx q[2];
rz(1.3428549) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0077483245) q[1];
sx q[1];
rz(-0.40002003) q[1];
sx q[1];
rz(-0.96537867) q[1];
rz(-pi) q[2];
rz(-0.19017724) q[3];
sx q[3];
rz(-1.7253644) q[3];
sx q[3];
rz(3.0319253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.74956191) q[2];
sx q[2];
rz(-1.6101086) q[2];
sx q[2];
rz(2.5677666) q[2];
rz(3.0841893) q[3];
sx q[3];
rz(-1.6618238) q[3];
sx q[3];
rz(-0.81006947) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60657984) q[0];
sx q[0];
rz(-2.8856394) q[0];
sx q[0];
rz(-2.8888597) q[0];
rz(-0.17929721) q[1];
sx q[1];
rz(-1.7959692) q[1];
sx q[1];
rz(1.4089233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63114595) q[0];
sx q[0];
rz(-2.8928693) q[0];
sx q[0];
rz(-1.2653217) q[0];
rz(-0.42667146) q[2];
sx q[2];
rz(-1.2152078) q[2];
sx q[2];
rz(-1.5384962) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2329916) q[1];
sx q[1];
rz(-1.4484693) q[1];
sx q[1];
rz(-0.78559383) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6176164) q[3];
sx q[3];
rz(-1.8856388) q[3];
sx q[3];
rz(3.063692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0130284) q[2];
sx q[2];
rz(-1.3821673) q[2];
sx q[2];
rz(-1.9405091) q[2];
rz(2.3342093) q[3];
sx q[3];
rz(-1.3599334) q[3];
sx q[3];
rz(2.0090296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078868911) q[0];
sx q[0];
rz(-1.4987334) q[0];
sx q[0];
rz(0.8412745) q[0];
rz(-1.9367283) q[1];
sx q[1];
rz(-1.8236022) q[1];
sx q[1];
rz(2.1119609) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1618274) q[0];
sx q[0];
rz(-2.0812391) q[0];
sx q[0];
rz(-0.011988601) q[0];
rz(-pi) q[1];
rz(-0.51635833) q[2];
sx q[2];
rz(-1.8022924) q[2];
sx q[2];
rz(1.6891434) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8334044) q[1];
sx q[1];
rz(-0.17886111) q[1];
sx q[1];
rz(3.0977644) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7577517) q[3];
sx q[3];
rz(-1.2724981) q[3];
sx q[3];
rz(-1.3307856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.131375) q[2];
sx q[2];
rz(-2.3848332) q[2];
sx q[2];
rz(-0.49794623) q[2];
rz(0.52305269) q[3];
sx q[3];
rz(-1.7224576) q[3];
sx q[3];
rz(-0.12652346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.2987357) q[0];
sx q[0];
rz(-3.0099478) q[0];
sx q[0];
rz(0.81197062) q[0];
rz(0.89093351) q[1];
sx q[1];
rz(-1.1633326) q[1];
sx q[1];
rz(2.1422211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10103664) q[0];
sx q[0];
rz(-2.479739) q[0];
sx q[0];
rz(-2.3836813) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81951253) q[2];
sx q[2];
rz(-0.14061059) q[2];
sx q[2];
rz(-2.484172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2216894) q[1];
sx q[1];
rz(-1.0070325) q[1];
sx q[1];
rz(-0.99181045) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29045297) q[3];
sx q[3];
rz(-0.92602611) q[3];
sx q[3];
rz(1.5034663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31819764) q[2];
sx q[2];
rz(-2.2237033) q[2];
sx q[2];
rz(-1.2981752) q[2];
rz(3.1031109) q[3];
sx q[3];
rz(-1.1063856) q[3];
sx q[3];
rz(-1.5255671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0865974) q[0];
sx q[0];
rz(-1.1308068) q[0];
sx q[0];
rz(0.31563345) q[0];
rz(-1.9075958) q[1];
sx q[1];
rz(-0.71593586) q[1];
sx q[1];
rz(-1.2589781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6427073) q[0];
sx q[0];
rz(-0.85432893) q[0];
sx q[0];
rz(1.1906719) q[0];
x q[1];
rz(0.18017144) q[2];
sx q[2];
rz(-2.368026) q[2];
sx q[2];
rz(-0.13828466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.054112) q[1];
sx q[1];
rz(-1.1743944) q[1];
sx q[1];
rz(1.1328843) q[1];
rz(-pi) q[2];
rz(-3.0397814) q[3];
sx q[3];
rz(-1.57668) q[3];
sx q[3];
rz(2.7525097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.2235609) q[2];
sx q[2];
rz(-0.74350244) q[2];
sx q[2];
rz(-0.027916748) q[2];
rz(-0.74808407) q[3];
sx q[3];
rz(-1.4192162) q[3];
sx q[3];
rz(-0.15302756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41392031) q[0];
sx q[0];
rz(-2.1509009) q[0];
sx q[0];
rz(-0.45898166) q[0];
rz(2.0052295) q[1];
sx q[1];
rz(-1.4308948) q[1];
sx q[1];
rz(-0.21839011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1735575) q[0];
sx q[0];
rz(-1.5941248) q[0];
sx q[0];
rz(-1.5541881) q[0];
x q[1];
rz(1.9543477) q[2];
sx q[2];
rz(-1.4462894) q[2];
sx q[2];
rz(2.2030168) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.23109197) q[1];
sx q[1];
rz(-2.0606961) q[1];
sx q[1];
rz(1.4040158) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1827385) q[3];
sx q[3];
rz(-1.4493296) q[3];
sx q[3];
rz(-1.1282008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8454664) q[2];
sx q[2];
rz(-0.25456905) q[2];
sx q[2];
rz(-1.0477585) q[2];
rz(2.3962077) q[3];
sx q[3];
rz(-1.5630009) q[3];
sx q[3];
rz(-1.6489702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1202241) q[0];
sx q[0];
rz(-0.57149082) q[0];
sx q[0];
rz(-0.34287232) q[0];
rz(2.951237) q[1];
sx q[1];
rz(-1.0089259) q[1];
sx q[1];
rz(0.51317936) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8489055) q[0];
sx q[0];
rz(-0.93528895) q[0];
sx q[0];
rz(0.19750316) q[0];
x q[1];
rz(2.076976) q[2];
sx q[2];
rz(-1.8864041) q[2];
sx q[2];
rz(-3.0402355) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.045195508) q[1];
sx q[1];
rz(-1.1151114) q[1];
sx q[1];
rz(-1.1229188) q[1];
x q[2];
rz(-2.9997708) q[3];
sx q[3];
rz(-1.98787) q[3];
sx q[3];
rz(2.6399222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43789431) q[2];
sx q[2];
rz(-2.5056705) q[2];
sx q[2];
rz(-0.61100125) q[2];
rz(0.25887394) q[3];
sx q[3];
rz(-1.2114108) q[3];
sx q[3];
rz(-1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7543058) q[0];
sx q[0];
rz(-1.5911234) q[0];
sx q[0];
rz(1.57244) q[0];
rz(-2.1517131) q[1];
sx q[1];
rz(-1.2073333) q[1];
sx q[1];
rz(-2.5055199) q[1];
rz(1.5237332) q[2];
sx q[2];
rz(-2.9527391) q[2];
sx q[2];
rz(-2.7486026) q[2];
rz(2.4425735) q[3];
sx q[3];
rz(-1.5200281) q[3];
sx q[3];
rz(2.716223) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
