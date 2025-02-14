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
rz(0.94035971) q[0];
sx q[0];
rz(-2.0191963) q[0];
sx q[0];
rz(2.9474592) q[0];
rz(-2.0464719) q[1];
sx q[1];
rz(-1.4580589) q[1];
sx q[1];
rz(-2.3966052) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4177822) q[0];
sx q[0];
rz(-0.11427721) q[0];
sx q[0];
rz(2.6835915) q[0];
rz(-0.92496867) q[2];
sx q[2];
rz(-0.87069521) q[2];
sx q[2];
rz(-2.2440804) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6319943) q[1];
sx q[1];
rz(-0.26226362) q[1];
sx q[1];
rz(-2.7523976) q[1];
x q[2];
rz(-2.8005014) q[3];
sx q[3];
rz(-2.00019) q[3];
sx q[3];
rz(-1.856696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8038883) q[2];
sx q[2];
rz(-1.8053728) q[2];
sx q[2];
rz(1.814092) q[2];
rz(-1.368807) q[3];
sx q[3];
rz(-2.8076706) q[3];
sx q[3];
rz(-0.22148618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1526445) q[0];
sx q[0];
rz(-1.7247609) q[0];
sx q[0];
rz(3.0864571) q[0];
rz(-0.70471835) q[1];
sx q[1];
rz(-1.5635419) q[1];
sx q[1];
rz(2.096874) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25356141) q[0];
sx q[0];
rz(-0.46875254) q[0];
sx q[0];
rz(-2.5028489) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98590322) q[2];
sx q[2];
rz(-2.0587181) q[2];
sx q[2];
rz(3.121738) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.738743) q[1];
sx q[1];
rz(-0.84058769) q[1];
sx q[1];
rz(2.5806246) q[1];
x q[2];
rz(2.5651188) q[3];
sx q[3];
rz(-3.0286791) q[3];
sx q[3];
rz(2.9993331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.81950554) q[2];
sx q[2];
rz(-2.4884067) q[2];
sx q[2];
rz(1.8196003) q[2];
rz(-2.3540438) q[3];
sx q[3];
rz(-1.4044263) q[3];
sx q[3];
rz(1.0549841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9287441) q[0];
sx q[0];
rz(-0.66900122) q[0];
sx q[0];
rz(-0.71841946) q[0];
rz(0.33572117) q[1];
sx q[1];
rz(-1.6751143) q[1];
sx q[1];
rz(-2.1885923) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7756726) q[0];
sx q[0];
rz(-1.6073415) q[0];
sx q[0];
rz(-2.1120872) q[0];
rz(-pi) q[1];
rz(-1.0510873) q[2];
sx q[2];
rz(-1.698709) q[2];
sx q[2];
rz(-3.0232883) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53297199) q[1];
sx q[1];
rz(-2.2257881) q[1];
sx q[1];
rz(2.6650732) q[1];
rz(-pi) q[2];
rz(2.5312214) q[3];
sx q[3];
rz(-1.1412176) q[3];
sx q[3];
rz(0.80738941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0042051729) q[2];
sx q[2];
rz(-1.1342099) q[2];
sx q[2];
rz(2.4887264) q[2];
rz(2.3434434) q[3];
sx q[3];
rz(-1.8316725) q[3];
sx q[3];
rz(0.64364141) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2586486) q[0];
sx q[0];
rz(-0.84027165) q[0];
sx q[0];
rz(2.3249481) q[0];
rz(0.69149292) q[1];
sx q[1];
rz(-1.34812) q[1];
sx q[1];
rz(-0.74849558) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2484528) q[0];
sx q[0];
rz(-1.7692385) q[0];
sx q[0];
rz(-0.31427419) q[0];
rz(-pi) q[1];
rz(1.9709936) q[2];
sx q[2];
rz(-2.9098401) q[2];
sx q[2];
rz(2.0830926) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.15259057) q[1];
sx q[1];
rz(-2.3034597) q[1];
sx q[1];
rz(0.92392606) q[1];
rz(-pi) q[2];
rz(-2.9660712) q[3];
sx q[3];
rz(-1.3427882) q[3];
sx q[3];
rz(-1.0358178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3889435) q[2];
sx q[2];
rz(-0.77072531) q[2];
sx q[2];
rz(-0.79461092) q[2];
rz(-1.4258344) q[3];
sx q[3];
rz(-2.1013997) q[3];
sx q[3];
rz(0.37253255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36885095) q[0];
sx q[0];
rz(-1.0798825) q[0];
sx q[0];
rz(0.4976196) q[0];
rz(2.3720062) q[1];
sx q[1];
rz(-1.9895357) q[1];
sx q[1];
rz(2.41113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4896547) q[0];
sx q[0];
rz(-1.7853072) q[0];
sx q[0];
rz(0.69846054) q[0];
rz(1.2271757) q[2];
sx q[2];
rz(-2.4848652) q[2];
sx q[2];
rz(-1.3948139) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2368612) q[1];
sx q[1];
rz(-1.7805459) q[1];
sx q[1];
rz(0.81717234) q[1];
rz(-pi) q[2];
rz(-0.83010556) q[3];
sx q[3];
rz(-2.5200994) q[3];
sx q[3];
rz(1.1661539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2490425) q[2];
sx q[2];
rz(-1.5183134) q[2];
sx q[2];
rz(2.6882233) q[2];
rz(-1.310965) q[3];
sx q[3];
rz(-2.9679208) q[3];
sx q[3];
rz(3.0176676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740358) q[0];
sx q[0];
rz(-2.0581764) q[0];
sx q[0];
rz(1.3483082) q[0];
rz(1.096161) q[1];
sx q[1];
rz(-1.4827012) q[1];
sx q[1];
rz(-2.5199264) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8529315) q[0];
sx q[0];
rz(-0.71446361) q[0];
sx q[0];
rz(-3.1174401) q[0];
rz(-pi) q[1];
rz(-2.9078324) q[2];
sx q[2];
rz(-0.40567652) q[2];
sx q[2];
rz(-2.7004227) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.409103) q[1];
sx q[1];
rz(-1.81899) q[1];
sx q[1];
rz(2.3677738) q[1];
x q[2];
rz(2.4824941) q[3];
sx q[3];
rz(-2.3602398) q[3];
sx q[3];
rz(2.3574061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2094476) q[2];
sx q[2];
rz(-0.79938447) q[2];
sx q[2];
rz(1.4383593) q[2];
rz(-0.98073331) q[3];
sx q[3];
rz(-1.8756198) q[3];
sx q[3];
rz(-1.2119306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.88406968) q[0];
sx q[0];
rz(-1.7634089) q[0];
sx q[0];
rz(-2.8164465) q[0];
rz(-2.6349321) q[1];
sx q[1];
rz(-1.0136565) q[1];
sx q[1];
rz(1.8258757) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0804101) q[0];
sx q[0];
rz(-1.4482486) q[0];
sx q[0];
rz(2.2558801) q[0];
rz(-pi) q[1];
rz(-2.1998134) q[2];
sx q[2];
rz(-2.468839) q[2];
sx q[2];
rz(0.37601177) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3218477) q[1];
sx q[1];
rz(-0.79230601) q[1];
sx q[1];
rz(2.5755432) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0990155) q[3];
sx q[3];
rz(-2.2025962) q[3];
sx q[3];
rz(-2.7160983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.331984) q[2];
sx q[2];
rz(-1.0386244) q[2];
sx q[2];
rz(1.5254376) q[2];
rz(-1.7841313) q[3];
sx q[3];
rz(-1.6922035) q[3];
sx q[3];
rz(-1.6239032) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9759393) q[0];
sx q[0];
rz(-1.140056) q[0];
sx q[0];
rz(-0.64919382) q[0];
rz(-0.80750418) q[1];
sx q[1];
rz(-2.0048001) q[1];
sx q[1];
rz(-1.3884707) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20304414) q[0];
sx q[0];
rz(-1.6106092) q[0];
sx q[0];
rz(-2.1757102) q[0];
rz(-pi) q[1];
rz(1.2818579) q[2];
sx q[2];
rz(-1.7035489) q[2];
sx q[2];
rz(-1.3049187) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.35508868) q[1];
sx q[1];
rz(-1.5032217) q[1];
sx q[1];
rz(1.3994292) q[1];
rz(-pi) q[2];
rz(0.47775538) q[3];
sx q[3];
rz(-0.43583187) q[3];
sx q[3];
rz(1.2660668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2041152) q[2];
sx q[2];
rz(-0.14675879) q[2];
sx q[2];
rz(0.80782962) q[2];
rz(0.67148036) q[3];
sx q[3];
rz(-1.6356133) q[3];
sx q[3];
rz(1.5988662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2272187) q[0];
sx q[0];
rz(-0.53633538) q[0];
sx q[0];
rz(-2.6887023) q[0];
rz(0.29531404) q[1];
sx q[1];
rz(-1.2295877) q[1];
sx q[1];
rz(-1.7426851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5495396) q[0];
sx q[0];
rz(-0.73707923) q[0];
sx q[0];
rz(-1.4455185) q[0];
x q[1];
rz(-1.8434502) q[2];
sx q[2];
rz(-0.52992199) q[2];
sx q[2];
rz(-0.31604813) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8321243) q[1];
sx q[1];
rz(-1.1946284) q[1];
sx q[1];
rz(2.6539567) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4612551) q[3];
sx q[3];
rz(-0.19831144) q[3];
sx q[3];
rz(1.2480145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7704775) q[2];
sx q[2];
rz(-2.2687843) q[2];
sx q[2];
rz(0.57658833) q[2];
rz(2.9234486) q[3];
sx q[3];
rz(-2.348867) q[3];
sx q[3];
rz(-1.920759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.6169154) q[0];
sx q[0];
rz(-0.24743947) q[0];
sx q[0];
rz(-3.0531378) q[0];
rz(0.5312008) q[1];
sx q[1];
rz(-0.4041268) q[1];
sx q[1];
rz(-1.521135) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5514497) q[0];
sx q[0];
rz(-1.9634754) q[0];
sx q[0];
rz(-3.0513502) q[0];
rz(-pi) q[1];
rz(2.0481061) q[2];
sx q[2];
rz(-1.6132406) q[2];
sx q[2];
rz(2.1316656) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6437592) q[1];
sx q[1];
rz(-1.6439142) q[1];
sx q[1];
rz(-0.22511367) q[1];
rz(-1.1559256) q[3];
sx q[3];
rz(-1.0153099) q[3];
sx q[3];
rz(-1.6270571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8497808) q[2];
sx q[2];
rz(-1.3498638) q[2];
sx q[2];
rz(1.527171) q[2];
rz(1.269086) q[3];
sx q[3];
rz(-2.1554155) q[3];
sx q[3];
rz(2.0932978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.108718) q[0];
sx q[0];
rz(-1.1058818) q[0];
sx q[0];
rz(-0.82431128) q[0];
rz(2.6678008) q[1];
sx q[1];
rz(-1.5722678) q[1];
sx q[1];
rz(-1.5617465) q[1];
rz(-0.59745477) q[2];
sx q[2];
rz(-1.93004) q[2];
sx q[2];
rz(-1.5782028) q[2];
rz(2.577435) q[3];
sx q[3];
rz(-2.0750792) q[3];
sx q[3];
rz(0.32479494) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
