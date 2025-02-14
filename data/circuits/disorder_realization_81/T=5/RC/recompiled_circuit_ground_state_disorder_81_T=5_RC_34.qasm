OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.12299744) q[0];
sx q[0];
rz(-1.7641492) q[0];
sx q[0];
rz(-2.7121845) q[0];
rz(1.6027066) q[1];
sx q[1];
rz(-1.4338926) q[1];
sx q[1];
rz(-1.5996999) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60889757) q[0];
sx q[0];
rz(-2.7809394) q[0];
sx q[0];
rz(0.92443539) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6060271) q[2];
sx q[2];
rz(-1.4270947) q[2];
sx q[2];
rz(-1.8409539) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1356537) q[1];
sx q[1];
rz(-2.7126813) q[1];
sx q[1];
rz(0.072838293) q[1];
rz(0.67486169) q[3];
sx q[3];
rz(-0.45025846) q[3];
sx q[3];
rz(-0.78801149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8371007) q[2];
sx q[2];
rz(-1.203275) q[2];
sx q[2];
rz(-0.692918) q[2];
rz(-0.20017008) q[3];
sx q[3];
rz(-0.19014159) q[3];
sx q[3];
rz(-2.0076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3226586) q[0];
sx q[0];
rz(-0.19911961) q[0];
sx q[0];
rz(-1.06485) q[0];
rz(-0.62243593) q[1];
sx q[1];
rz(-1.348) q[1];
sx q[1];
rz(-2.650824) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3929185) q[0];
sx q[0];
rz(-3.1160917) q[0];
sx q[0];
rz(-0.37125094) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42224531) q[2];
sx q[2];
rz(-1.7755055) q[2];
sx q[2];
rz(-1.4377571) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.789114) q[1];
sx q[1];
rz(-1.4755469) q[1];
sx q[1];
rz(0.75525166) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4014259) q[3];
sx q[3];
rz(-1.5233831) q[3];
sx q[3];
rz(-2.5518038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6530767) q[2];
sx q[2];
rz(-1.642903) q[2];
sx q[2];
rz(0.24031362) q[2];
rz(-0.49992418) q[3];
sx q[3];
rz(-2.5559055) q[3];
sx q[3];
rz(1.2691931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8847454) q[0];
sx q[0];
rz(-0.95504967) q[0];
sx q[0];
rz(0.98168674) q[0];
rz(0.49199545) q[1];
sx q[1];
rz(-1.0849846) q[1];
sx q[1];
rz(-1.6384151) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5444191) q[0];
sx q[0];
rz(-2.0453296) q[0];
sx q[0];
rz(-2.3719792) q[0];
rz(3.0122224) q[2];
sx q[2];
rz(-0.40310848) q[2];
sx q[2];
rz(1.0327686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99276421) q[1];
sx q[1];
rz(-1.1849313) q[1];
sx q[1];
rz(-2.7817315) q[1];
x q[2];
rz(-1.7812626) q[3];
sx q[3];
rz(-0.099848824) q[3];
sx q[3];
rz(2.9304402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65695277) q[2];
sx q[2];
rz(-0.52565614) q[2];
sx q[2];
rz(3.0204115) q[2];
rz(2.1335404) q[3];
sx q[3];
rz(-1.1153778) q[3];
sx q[3];
rz(-2.0602267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0693531) q[0];
sx q[0];
rz(-0.00088748137) q[0];
sx q[0];
rz(-2.6278611) q[0];
rz(0.95798245) q[1];
sx q[1];
rz(-1.999141) q[1];
sx q[1];
rz(1.6987919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.643647) q[0];
sx q[0];
rz(-1.4063324) q[0];
sx q[0];
rz(2.7269159) q[0];
rz(2.7755205) q[2];
sx q[2];
rz(-2.5806081) q[2];
sx q[2];
rz(1.6268886) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2580766) q[1];
sx q[1];
rz(-0.95702584) q[1];
sx q[1];
rz(0.10054904) q[1];
rz(-pi) q[2];
rz(0.6616627) q[3];
sx q[3];
rz(-1.6473624) q[3];
sx q[3];
rz(0.14684248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.55502597) q[2];
sx q[2];
rz(-2.1777966) q[2];
sx q[2];
rz(1.2653992) q[2];
rz(-1.966656) q[3];
sx q[3];
rz(-0.73052162) q[3];
sx q[3];
rz(-1.9780212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88605276) q[0];
sx q[0];
rz(-0.20296725) q[0];
sx q[0];
rz(1.8170005) q[0];
rz(-0.96877226) q[1];
sx q[1];
rz(-1.4854393) q[1];
sx q[1];
rz(1.0383777) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2823855) q[0];
sx q[0];
rz(-2.7250054) q[0];
sx q[0];
rz(-2.7340915) q[0];
x q[1];
rz(-2.3164301) q[2];
sx q[2];
rz(-1.0010011) q[2];
sx q[2];
rz(-1.93376) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9818791) q[1];
sx q[1];
rz(-0.37283373) q[1];
sx q[1];
rz(-3.0283958) q[1];
x q[2];
rz(-2.5143751) q[3];
sx q[3];
rz(-2.3786847) q[3];
sx q[3];
rz(-1.4351532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6394627) q[2];
sx q[2];
rz(-2.837193) q[2];
sx q[2];
rz(-0.89140618) q[2];
rz(-1.8799479) q[3];
sx q[3];
rz(-1.5104048) q[3];
sx q[3];
rz(2.4752899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0618184) q[0];
sx q[0];
rz(-0.49538716) q[0];
sx q[0];
rz(0.99639446) q[0];
rz(2.2415316) q[1];
sx q[1];
rz(-0.78949094) q[1];
sx q[1];
rz(2.9642504) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3729438) q[0];
sx q[0];
rz(-1.5710314) q[0];
sx q[0];
rz(1.7223486) q[0];
rz(1.5992237) q[2];
sx q[2];
rz(-2.5027983) q[2];
sx q[2];
rz(2.6306689) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4882433) q[1];
sx q[1];
rz(-1.0945787) q[1];
sx q[1];
rz(-2.819092) q[1];
rz(1.7436036) q[3];
sx q[3];
rz(-2.6263642) q[3];
sx q[3];
rz(0.15095161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8511054) q[2];
sx q[2];
rz(-1.1815716) q[2];
sx q[2];
rz(2.8372852) q[2];
rz(2.815222) q[3];
sx q[3];
rz(-2.5984952) q[3];
sx q[3];
rz(0.91199005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07311634) q[0];
sx q[0];
rz(-0.40962064) q[0];
sx q[0];
rz(-2.2090744) q[0];
rz(0.069123507) q[1];
sx q[1];
rz(-0.2176452) q[1];
sx q[1];
rz(-2.1014012) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0737586) q[0];
sx q[0];
rz(-0.76601765) q[0];
sx q[0];
rz(-0.037184663) q[0];
rz(-1.8273152) q[2];
sx q[2];
rz(-1.7396915) q[2];
sx q[2];
rz(-2.0467364) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.8833163) q[1];
sx q[1];
rz(-0.81352106) q[1];
sx q[1];
rz(-1.3431647) q[1];
rz(-pi) q[2];
rz(-2.4428818) q[3];
sx q[3];
rz(-1.1160399) q[3];
sx q[3];
rz(0.48331279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13504623) q[2];
sx q[2];
rz(-2.4032205) q[2];
sx q[2];
rz(-2.4086003) q[2];
rz(-2.0586355) q[3];
sx q[3];
rz(-2.0854009) q[3];
sx q[3];
rz(-2.1347031) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55423823) q[0];
sx q[0];
rz(-0.08389689) q[0];
sx q[0];
rz(-2.6485637) q[0];
rz(2.9250277) q[1];
sx q[1];
rz(-2.000587) q[1];
sx q[1];
rz(2.1176178) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1676583) q[0];
sx q[0];
rz(-1.5588781) q[0];
sx q[0];
rz(2.5456136) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31957133) q[2];
sx q[2];
rz(-2.7929578) q[2];
sx q[2];
rz(2.2120668) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1962016) q[1];
sx q[1];
rz(-1.7111254) q[1];
sx q[1];
rz(-0.36044557) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29924691) q[3];
sx q[3];
rz(-0.95029921) q[3];
sx q[3];
rz(0.79122358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0563125) q[2];
sx q[2];
rz(-1.6322501) q[2];
sx q[2];
rz(2.4658266) q[2];
rz(2.7285649) q[3];
sx q[3];
rz(-1.690381) q[3];
sx q[3];
rz(0.54615027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6290879) q[0];
sx q[0];
rz(-0.65852037) q[0];
sx q[0];
rz(3.107048) q[0];
rz(-0.054556219) q[1];
sx q[1];
rz(-1.9787534) q[1];
sx q[1];
rz(-2.3202855) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9050582) q[0];
sx q[0];
rz(-0.58249677) q[0];
sx q[0];
rz(-2.9597917) q[0];
rz(-pi) q[1];
rz(-0.088772687) q[2];
sx q[2];
rz(-1.5800416) q[2];
sx q[2];
rz(-0.68999664) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.37335884) q[1];
sx q[1];
rz(-1.7204509) q[1];
sx q[1];
rz(-1.6878769) q[1];
rz(-0.48252941) q[3];
sx q[3];
rz(-2.165613) q[3];
sx q[3];
rz(1.1728668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.36046946) q[2];
sx q[2];
rz(-2.1253773) q[2];
sx q[2];
rz(0.13295573) q[2];
rz(-0.41108701) q[3];
sx q[3];
rz(-1.6229595) q[3];
sx q[3];
rz(-1.0857922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0947615) q[0];
sx q[0];
rz(-2.5346041) q[0];
sx q[0];
rz(1.2588311) q[0];
rz(-1.6350485) q[1];
sx q[1];
rz(-2.4848487) q[1];
sx q[1];
rz(-2.3283995) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98018089) q[0];
sx q[0];
rz(-1.4927571) q[0];
sx q[0];
rz(-1.6059931) q[0];
rz(-pi) q[1];
rz(-1.4237464) q[2];
sx q[2];
rz(-0.67710256) q[2];
sx q[2];
rz(-2.8092217) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.836536) q[1];
sx q[1];
rz(-0.99545762) q[1];
sx q[1];
rz(1.025455) q[1];
rz(-0.45222262) q[3];
sx q[3];
rz(-1.4503696) q[3];
sx q[3];
rz(-2.6129006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64858156) q[2];
sx q[2];
rz(-2.0040671) q[2];
sx q[2];
rz(2.0900334) q[2];
rz(-1.018853) q[3];
sx q[3];
rz(-0.84210432) q[3];
sx q[3];
rz(2.4837608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1558341) q[0];
sx q[0];
rz(-2.4875165) q[0];
sx q[0];
rz(-2.2610337) q[0];
rz(-1.8219933) q[1];
sx q[1];
rz(-1.1450014) q[1];
sx q[1];
rz(-3.0897279) q[1];
rz(-2.9319977) q[2];
sx q[2];
rz(-1.4996698) q[2];
sx q[2];
rz(-1.3347129) q[2];
rz(2.9444957) q[3];
sx q[3];
rz(-1.5119678) q[3];
sx q[3];
rz(1.5869535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
