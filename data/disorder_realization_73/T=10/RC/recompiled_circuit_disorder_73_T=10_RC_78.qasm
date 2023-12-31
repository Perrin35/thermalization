OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9956545) q[0];
sx q[0];
rz(-0.50322682) q[0];
sx q[0];
rz(-0.72416645) q[0];
rz(-2.5016298) q[1];
sx q[1];
rz(-2.6115186) q[1];
sx q[1];
rz(-2.35676) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7895176) q[0];
sx q[0];
rz(-1.8624458) q[0];
sx q[0];
rz(-1.350499) q[0];
x q[1];
rz(-2.7230524) q[2];
sx q[2];
rz(-1.6633908) q[2];
sx q[2];
rz(1.6469524) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6341056) q[1];
sx q[1];
rz(-1.7934985) q[1];
sx q[1];
rz(2.1102064) q[1];
rz(-1.647244) q[3];
sx q[3];
rz(-2.7259698) q[3];
sx q[3];
rz(1.3941744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.589754) q[2];
sx q[2];
rz(-1.4171615) q[2];
sx q[2];
rz(0.067967728) q[2];
rz(0.12456482) q[3];
sx q[3];
rz(-2.8187276) q[3];
sx q[3];
rz(1.7547866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92000604) q[0];
sx q[0];
rz(-3.0060372) q[0];
sx q[0];
rz(2.8979229) q[0];
rz(0.63175732) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(-1.7858645) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1649363) q[0];
sx q[0];
rz(-1.5436085) q[0];
sx q[0];
rz(-1.3154161) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83265702) q[2];
sx q[2];
rz(-0.50561935) q[2];
sx q[2];
rz(2.6308699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.020097453) q[1];
sx q[1];
rz(-2.0165682) q[1];
sx q[1];
rz(-2.9989468) q[1];
rz(-0.59216604) q[3];
sx q[3];
rz(-1.5321931) q[3];
sx q[3];
rz(1.6174699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0791066) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(2.8919354) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(-0.33199582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-2.8963985) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(2.2431592) q[0];
rz(-1.8067182) q[1];
sx q[1];
rz(-1.9060262) q[1];
sx q[1];
rz(1.2737087) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1153591) q[0];
sx q[0];
rz(-0.44239487) q[0];
sx q[0];
rz(-0.94985234) q[0];
rz(-pi) q[1];
rz(1.4726228) q[2];
sx q[2];
rz(-1.2766826) q[2];
sx q[2];
rz(-1.836118) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85325235) q[1];
sx q[1];
rz(-0.77008343) q[1];
sx q[1];
rz(-0.5132765) q[1];
x q[2];
rz(2.5898315) q[3];
sx q[3];
rz(-2.3380087) q[3];
sx q[3];
rz(1.4078275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0818103) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(-0.95345062) q[2];
rz(-0.034514286) q[3];
sx q[3];
rz(-2.3551066) q[3];
sx q[3];
rz(0.22687337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8614486) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(2.8934073) q[0];
rz(-1.0379627) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(0.074137069) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.603133) q[0];
sx q[0];
rz(-0.54766253) q[0];
sx q[0];
rz(1.0663701) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2494227) q[2];
sx q[2];
rz(-1.2954419) q[2];
sx q[2];
rz(-2.8868669) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76064199) q[1];
sx q[1];
rz(-1.4529072) q[1];
sx q[1];
rz(1.8878493) q[1];
x q[2];
rz(0.24012633) q[3];
sx q[3];
rz(-1.4014981) q[3];
sx q[3];
rz(1.5508159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2512102) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(-0.038643535) q[2];
rz(2.1679227) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53428179) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(-1.3624396) q[0];
rz(-2.3249987) q[1];
sx q[1];
rz(-1.8530308) q[1];
sx q[1];
rz(1.978925) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8487932) q[0];
sx q[0];
rz(-1.6822364) q[0];
sx q[0];
rz(-0.029990002) q[0];
rz(-0.16935279) q[2];
sx q[2];
rz(-1.8893818) q[2];
sx q[2];
rz(0.40118518) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.34449023) q[1];
sx q[1];
rz(-1.3791729) q[1];
sx q[1];
rz(-1.9413129) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88697042) q[3];
sx q[3];
rz(-0.63347048) q[3];
sx q[3];
rz(2.5630643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5148619) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(2.999372) q[2];
rz(0.90406117) q[3];
sx q[3];
rz(-1.8217434) q[3];
sx q[3];
rz(-2.952125) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11480039) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(1.6739155) q[0];
rz(0.57178512) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(-2.8335559) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338617) q[0];
sx q[0];
rz(-2.9836285) q[0];
sx q[0];
rz(-2.4979742) q[0];
x q[1];
rz(0.96254827) q[2];
sx q[2];
rz(-1.0481917) q[2];
sx q[2];
rz(0.42524291) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.088965485) q[1];
sx q[1];
rz(-1.7672156) q[1];
sx q[1];
rz(0.88167015) q[1];
rz(2.7932348) q[3];
sx q[3];
rz(-1.6568686) q[3];
sx q[3];
rz(2.0406046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3623111) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(-1.1479088) q[2];
rz(2.4273196) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(0.64546293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66184735) q[0];
sx q[0];
rz(-0.84091887) q[0];
sx q[0];
rz(-3.0116144) q[0];
rz(0.030844363) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(-2.470509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13043159) q[0];
sx q[0];
rz(-1.6058308) q[0];
sx q[0];
rz(-0.053671562) q[0];
rz(-1.7550049) q[2];
sx q[2];
rz(-1.6008018) q[2];
sx q[2];
rz(-1.9629994) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5741201) q[1];
sx q[1];
rz(-1.4680871) q[1];
sx q[1];
rz(1.1979539) q[1];
rz(-pi) q[2];
rz(0.89550771) q[3];
sx q[3];
rz(-1.5338147) q[3];
sx q[3];
rz(-1.0469588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0188296) q[2];
sx q[2];
rz(-1.5315703) q[2];
sx q[2];
rz(-2.5637131) q[2];
rz(0.028586483) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(-1.2602497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9776483) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(-0.41241616) q[0];
rz(-1.6917797) q[1];
sx q[1];
rz(-1.342536) q[1];
sx q[1];
rz(-1.1669881) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3229423) q[0];
sx q[0];
rz(-2.4837821) q[0];
sx q[0];
rz(-1.2176745) q[0];
rz(-pi) q[1];
rz(-2.8655878) q[2];
sx q[2];
rz(-0.88062421) q[2];
sx q[2];
rz(1.3921757) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7653212) q[1];
sx q[1];
rz(-1.6011392) q[1];
sx q[1];
rz(-0.18860753) q[1];
rz(-2.1006881) q[3];
sx q[3];
rz(-1.4412672) q[3];
sx q[3];
rz(-1.1216175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.940544) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(2.9525625) q[2];
rz(-0.14686251) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(1.3930901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1259595) q[0];
sx q[0];
rz(-1.3110315) q[0];
sx q[0];
rz(0.92700672) q[0];
rz(-1.3828297) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(1.6519201) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531567) q[0];
sx q[0];
rz(-2.2373767) q[0];
sx q[0];
rz(0.96320926) q[0];
rz(-1.6696879) q[2];
sx q[2];
rz(-2.6675468) q[2];
sx q[2];
rz(-2.3941819) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5628964) q[1];
sx q[1];
rz(-1.3778731) q[1];
sx q[1];
rz(1.3955411) q[1];
rz(-2.1065815) q[3];
sx q[3];
rz(-2.0097369) q[3];
sx q[3];
rz(0.5700232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5902517) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(-1.7112973) q[2];
rz(2.5643505) q[3];
sx q[3];
rz(-0.8876628) q[3];
sx q[3];
rz(-1.757471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5532613) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(2.8531895) q[0];
rz(-2.6092031) q[1];
sx q[1];
rz(-2.6817697) q[1];
sx q[1];
rz(2.9945701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1823746) q[0];
sx q[0];
rz(-1.9518513) q[0];
sx q[0];
rz(-0.94840886) q[0];
rz(-pi) q[1];
rz(0.23994259) q[2];
sx q[2];
rz(-1.502617) q[2];
sx q[2];
rz(2.6477674) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.30565572) q[1];
sx q[1];
rz(-2.5218997) q[1];
sx q[1];
rz(-1.6395007) q[1];
rz(-2.5649928) q[3];
sx q[3];
rz(-1.7017662) q[3];
sx q[3];
rz(3.1104345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0816575) q[2];
sx q[2];
rz(-2.7066878) q[2];
sx q[2];
rz(2.3975513) q[2];
rz(0.75731164) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(1.7448759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.025678) q[0];
sx q[0];
rz(-1.0703351) q[0];
sx q[0];
rz(-1.0967789) q[0];
rz(2.3241282) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(0.4521162) q[2];
sx q[2];
rz(-1.631626) q[2];
sx q[2];
rz(0.22125868) q[2];
rz(-1.3656473) q[3];
sx q[3];
rz(-0.57153265) q[3];
sx q[3];
rz(2.3415757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
