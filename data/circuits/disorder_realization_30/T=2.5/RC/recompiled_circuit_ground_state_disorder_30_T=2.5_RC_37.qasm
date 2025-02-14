OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71830463) q[0];
sx q[0];
rz(2.365132) q[0];
sx q[0];
rz(9.790701) q[0];
rz(0.34691063) q[1];
sx q[1];
rz(-0.28471714) q[1];
sx q[1];
rz(1.7528344) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4735218) q[0];
sx q[0];
rz(-1.5041623) q[0];
sx q[0];
rz(0.20247831) q[0];
rz(1.2728782) q[2];
sx q[2];
rz(-0.90014825) q[2];
sx q[2];
rz(0.16042319) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5386913) q[1];
sx q[1];
rz(-2.5316892) q[1];
sx q[1];
rz(-2.2944488) q[1];
rz(-2.4284989) q[3];
sx q[3];
rz(-0.72848195) q[3];
sx q[3];
rz(-0.99458867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1846788) q[2];
sx q[2];
rz(-2.2214486) q[2];
sx q[2];
rz(2.206395) q[2];
rz(-0.30185559) q[3];
sx q[3];
rz(-1.5993092) q[3];
sx q[3];
rz(-2.7479318) q[3];
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
rz(pi/2) q[0];
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
rz(2.1644156) q[0];
sx q[0];
rz(-1.875832) q[0];
sx q[0];
rz(-2.6432977) q[0];
rz(1.7138819) q[1];
sx q[1];
rz(-1.2063113) q[1];
sx q[1];
rz(1.0867585) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5887077) q[0];
sx q[0];
rz(-0.7378788) q[0];
sx q[0];
rz(0.35564519) q[0];
rz(-1.837376) q[2];
sx q[2];
rz(-1.757903) q[2];
sx q[2];
rz(3.0361036) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3698235) q[1];
sx q[1];
rz(-0.46716094) q[1];
sx q[1];
rz(-2.7792756) q[1];
rz(-pi) q[2];
rz(2.8233775) q[3];
sx q[3];
rz(-1.5822142) q[3];
sx q[3];
rz(0.59225692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19865092) q[2];
sx q[2];
rz(-2.8309839) q[2];
sx q[2];
rz(-1.3322213) q[2];
rz(1.1474991) q[3];
sx q[3];
rz(-1.5787326) q[3];
sx q[3];
rz(1.8089627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20951095) q[0];
sx q[0];
rz(-1.4028343) q[0];
sx q[0];
rz(-0.15723666) q[0];
rz(1.2278185) q[1];
sx q[1];
rz(-0.20918748) q[1];
sx q[1];
rz(0.42207178) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2651208) q[0];
sx q[0];
rz(-2.9419964) q[0];
sx q[0];
rz(-1.100698) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3968758) q[2];
sx q[2];
rz(-1.2971767) q[2];
sx q[2];
rz(1.5762941) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5383407) q[1];
sx q[1];
rz(-1.1900239) q[1];
sx q[1];
rz(-0.095007665) q[1];
rz(-pi) q[2];
rz(-2.3396569) q[3];
sx q[3];
rz(-1.3462596) q[3];
sx q[3];
rz(0.79659407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32855365) q[2];
sx q[2];
rz(-2.1812794) q[2];
sx q[2];
rz(0.73406827) q[2];
rz(-1.0775393) q[3];
sx q[3];
rz(-0.68818337) q[3];
sx q[3];
rz(-0.075210007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5587191) q[0];
sx q[0];
rz(-1.0105157) q[0];
sx q[0];
rz(0.016481312) q[0];
rz(0.074542848) q[1];
sx q[1];
rz(-0.45769474) q[1];
sx q[1];
rz(-2.8864536) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0150361) q[0];
sx q[0];
rz(-1.714596) q[0];
sx q[0];
rz(2.6565927) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33489494) q[2];
sx q[2];
rz(-1.0303632) q[2];
sx q[2];
rz(1.188736) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.29683477) q[1];
sx q[1];
rz(-2.1898309) q[1];
sx q[1];
rz(1.0356139) q[1];
rz(-1.2665073) q[3];
sx q[3];
rz(-0.36484584) q[3];
sx q[3];
rz(-1.5740652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.61465803) q[2];
sx q[2];
rz(-0.69710985) q[2];
sx q[2];
rz(0.69980168) q[2];
rz(0.091863306) q[3];
sx q[3];
rz(-1.9242761) q[3];
sx q[3];
rz(-0.45472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5192473) q[0];
sx q[0];
rz(-0.82038251) q[0];
sx q[0];
rz(1.1118332) q[0];
rz(0.19019292) q[1];
sx q[1];
rz(-0.88514248) q[1];
sx q[1];
rz(-1.9817188) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5525813) q[0];
sx q[0];
rz(-1.9144626) q[0];
sx q[0];
rz(0.036618311) q[0];
x q[1];
rz(-1.3746475) q[2];
sx q[2];
rz(-2.4939459) q[2];
sx q[2];
rz(2.7638433) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9638744) q[1];
sx q[1];
rz(-1.8184806) q[1];
sx q[1];
rz(2.5801093) q[1];
x q[2];
rz(0.044951602) q[3];
sx q[3];
rz(-2.0922497) q[3];
sx q[3];
rz(1.5779881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5707034) q[2];
sx q[2];
rz(-0.8231701) q[2];
sx q[2];
rz(-2.9111351) q[2];
rz(1.0620091) q[3];
sx q[3];
rz(-1.2758723) q[3];
sx q[3];
rz(1.5084069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.688711) q[0];
sx q[0];
rz(-2.0827561) q[0];
sx q[0];
rz(1.7857312) q[0];
rz(-2.611825) q[1];
sx q[1];
rz(-2.3593088) q[1];
sx q[1];
rz(1.1518325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760107) q[0];
sx q[0];
rz(-1.9261203) q[0];
sx q[0];
rz(1.4655345) q[0];
x q[1];
rz(-1.4131613) q[2];
sx q[2];
rz(-0.54716483) q[2];
sx q[2];
rz(0.62842759) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85688389) q[1];
sx q[1];
rz(-1.9344866) q[1];
sx q[1];
rz(-1.3631352) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6883259) q[3];
sx q[3];
rz(-1.3580048) q[3];
sx q[3];
rz(-0.26373395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6666224) q[2];
sx q[2];
rz(-0.64133659) q[2];
sx q[2];
rz(-0.46710157) q[2];
rz(-1.0111672) q[3];
sx q[3];
rz(-2.7979388) q[3];
sx q[3];
rz(1.3694192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2924627) q[0];
sx q[0];
rz(-2.9865773) q[0];
sx q[0];
rz(-2.9169061) q[0];
rz(2.977773) q[1];
sx q[1];
rz(-0.33912173) q[1];
sx q[1];
rz(0.64635578) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40898973) q[0];
sx q[0];
rz(-1.3209131) q[0];
sx q[0];
rz(2.1134209) q[0];
x q[1];
rz(1.8936526) q[2];
sx q[2];
rz(-1.7778998) q[2];
sx q[2];
rz(-2.8139909) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.028370628) q[1];
sx q[1];
rz(-2.592855) q[1];
sx q[1];
rz(-0.1446677) q[1];
rz(-1.7067753) q[3];
sx q[3];
rz(-1.9549911) q[3];
sx q[3];
rz(0.7766436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15567638) q[2];
sx q[2];
rz(-1.4480269) q[2];
sx q[2];
rz(-2.7316366) q[2];
rz(-2.3112467) q[3];
sx q[3];
rz(-0.94873077) q[3];
sx q[3];
rz(1.693694) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.603867) q[0];
sx q[0];
rz(-0.69320885) q[0];
sx q[0];
rz(2.895288) q[0];
rz(-2.314997) q[1];
sx q[1];
rz(-1.3984503) q[1];
sx q[1];
rz(-2.8776317) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6851229) q[0];
sx q[0];
rz(-1.5780266) q[0];
sx q[0];
rz(-1.475322) q[0];
x q[1];
rz(0.47672456) q[2];
sx q[2];
rz(-0.84273224) q[2];
sx q[2];
rz(-1.5155033) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3166271) q[1];
sx q[1];
rz(-2.3967603) q[1];
sx q[1];
rz(-1.1492351) q[1];
x q[2];
rz(2.5313782) q[3];
sx q[3];
rz(-2.4596678) q[3];
sx q[3];
rz(-2.1206926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0366514) q[2];
sx q[2];
rz(-1.8737326) q[2];
sx q[2];
rz(2.4264753) q[2];
rz(-0.07130833) q[3];
sx q[3];
rz(-1.5569867) q[3];
sx q[3];
rz(-0.74688545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1320553) q[0];
sx q[0];
rz(-1.6522464) q[0];
sx q[0];
rz(-0.43553964) q[0];
rz(0.14777331) q[1];
sx q[1];
rz(-1.4316033) q[1];
sx q[1];
rz(1.6730283) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72318059) q[0];
sx q[0];
rz(-1.7505129) q[0];
sx q[0];
rz(-0.43850684) q[0];
rz(-pi) q[1];
rz(-0.69738241) q[2];
sx q[2];
rz(-1.5101523) q[2];
sx q[2];
rz(-1.9549119) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.33959162) q[1];
sx q[1];
rz(-0.42515426) q[1];
sx q[1];
rz(-1.9849066) q[1];
x q[2];
rz(2.7379509) q[3];
sx q[3];
rz(-1.3714681) q[3];
sx q[3];
rz(0.9404054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.396951) q[2];
sx q[2];
rz(-1.3775185) q[2];
sx q[2];
rz(0.89293876) q[2];
rz(-1.8630155) q[3];
sx q[3];
rz(-1.1568926) q[3];
sx q[3];
rz(1.6781835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29499149) q[0];
sx q[0];
rz(-2.8104267) q[0];
sx q[0];
rz(-2.5355329) q[0];
rz(-0.010295708) q[1];
sx q[1];
rz(-0.98774397) q[1];
sx q[1];
rz(-3.0616679) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2250076) q[0];
sx q[0];
rz(-1.8226624) q[0];
sx q[0];
rz(2.1504137) q[0];
x q[1];
rz(-0.84236161) q[2];
sx q[2];
rz(-1.7879221) q[2];
sx q[2];
rz(2.6922243) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2735168) q[1];
sx q[1];
rz(-0.40765992) q[1];
sx q[1];
rz(-3.0855387) q[1];
x q[2];
rz(-2.2977423) q[3];
sx q[3];
rz(-2.8371713) q[3];
sx q[3];
rz(-2.3850887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.14110485) q[2];
sx q[2];
rz(-2.2147369) q[2];
sx q[2];
rz(2.1263988) q[2];
rz(-2.5403533) q[3];
sx q[3];
rz(-0.70178086) q[3];
sx q[3];
rz(1.0465485) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4973608) q[0];
sx q[0];
rz(-1.5174706) q[0];
sx q[0];
rz(0.68751412) q[0];
rz(0.58912206) q[1];
sx q[1];
rz(-0.67710572) q[1];
sx q[1];
rz(2.6893375) q[1];
rz(-1.4745787) q[2];
sx q[2];
rz(-2.0141891) q[2];
sx q[2];
rz(-0.63465848) q[2];
rz(2.3169869) q[3];
sx q[3];
rz(-1.6578703) q[3];
sx q[3];
rz(-0.35005611) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
