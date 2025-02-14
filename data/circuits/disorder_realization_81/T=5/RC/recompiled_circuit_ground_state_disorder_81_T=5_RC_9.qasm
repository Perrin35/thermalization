OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0185952) q[0];
sx q[0];
rz(-1.3774435) q[0];
sx q[0];
rz(-0.42940816) q[0];
rz(1.6027066) q[1];
sx q[1];
rz(-1.4338926) q[1];
sx q[1];
rz(1.5418928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2875007) q[0];
sx q[0];
rz(-1.2852291) q[0];
sx q[0];
rz(0.2233611) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23878204) q[2];
sx q[2];
rz(-2.9936643) q[2];
sx q[2];
rz(-1.5419568) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0555901) q[1];
sx q[1];
rz(-1.9984954) q[1];
sx q[1];
rz(-1.5375288) q[1];
x q[2];
rz(-0.67486169) q[3];
sx q[3];
rz(-0.45025846) q[3];
sx q[3];
rz(-2.3535812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.30449197) q[2];
sx q[2];
rz(-1.203275) q[2];
sx q[2];
rz(0.692918) q[2];
rz(-2.9414226) q[3];
sx q[3];
rz(-0.19014159) q[3];
sx q[3];
rz(-1.1339124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-1.8189341) q[0];
sx q[0];
rz(-2.942473) q[0];
sx q[0];
rz(-2.0767427) q[0];
rz(2.5191567) q[1];
sx q[1];
rz(-1.7935926) q[1];
sx q[1];
rz(-0.49076864) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1932632) q[0];
sx q[0];
rz(-1.5800467) q[0];
sx q[0];
rz(3.1178283) q[0];
rz(-pi) q[1];
rz(-0.42224531) q[2];
sx q[2];
rz(-1.7755055) q[2];
sx q[2];
rz(1.7038356) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.789114) q[1];
sx q[1];
rz(-1.4755469) q[1];
sx q[1];
rz(-2.386341) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.048100483) q[3];
sx q[3];
rz(-1.7399746) q[3];
sx q[3];
rz(-2.1524803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48851594) q[2];
sx q[2];
rz(-1.4986897) q[2];
sx q[2];
rz(-2.901279) q[2];
rz(-2.6416685) q[3];
sx q[3];
rz(-0.58568716) q[3];
sx q[3];
rz(1.2691931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.2568473) q[0];
sx q[0];
rz(-0.95504967) q[0];
sx q[0];
rz(-2.1599059) q[0];
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
rz(-0.41427754) q[0];
sx q[0];
rz(-2.2637667) q[0];
sx q[0];
rz(-2.5056865) q[0];
rz(-pi) q[1];
x q[1];
rz(1.515834) q[2];
sx q[2];
rz(-1.171249) q[2];
sx q[2];
rz(1.9683226) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7780077) q[1];
sx q[1];
rz(-2.6202046) q[1];
sx q[1];
rz(-0.85659493) q[1];
rz(-1.4731367) q[3];
sx q[3];
rz(-1.5499695) q[3];
sx q[3];
rz(-1.5690924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4846399) q[2];
sx q[2];
rz(-2.6159365) q[2];
sx q[2];
rz(3.0204115) q[2];
rz(-2.1335404) q[3];
sx q[3];
rz(-1.1153778) q[3];
sx q[3];
rz(2.0602267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0722395) q[0];
sx q[0];
rz(-0.00088748137) q[0];
sx q[0];
rz(2.6278611) q[0];
rz(-0.95798245) q[1];
sx q[1];
rz(-1.1424516) q[1];
sx q[1];
rz(1.6987919) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0009036) q[0];
sx q[0];
rz(-1.1620518) q[0];
sx q[0];
rz(1.391414) q[0];
rz(-2.7755205) q[2];
sx q[2];
rz(-2.5806081) q[2];
sx q[2];
rz(1.5147041) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.084667) q[1];
sx q[1];
rz(-0.62090331) q[1];
sx q[1];
rz(-1.4292745) q[1];
x q[2];
rz(0.12421457) q[3];
sx q[3];
rz(-2.476177) q[3];
sx q[3];
rz(1.6197325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5865667) q[2];
sx q[2];
rz(-2.1777966) q[2];
sx q[2];
rz(-1.8761934) q[2];
rz(1.1749367) q[3];
sx q[3];
rz(-2.411071) q[3];
sx q[3];
rz(-1.1635715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2555399) q[0];
sx q[0];
rz(-2.9386254) q[0];
sx q[0];
rz(-1.3245921) q[0];
rz(2.1728204) q[1];
sx q[1];
rz(-1.6561534) q[1];
sx q[1];
rz(-1.0383777) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85920716) q[0];
sx q[0];
rz(-2.7250054) q[0];
sx q[0];
rz(-0.4075012) q[0];
rz(-pi) q[1];
rz(0.81401396) q[2];
sx q[2];
rz(-2.2377439) q[2];
sx q[2];
rz(-2.2499354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6250302) q[1];
sx q[1];
rz(-1.5296401) q[1];
sx q[1];
rz(-0.37066083) q[1];
x q[2];
rz(2.5143751) q[3];
sx q[3];
rz(-0.76290799) q[3];
sx q[3];
rz(-1.4351532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6394627) q[2];
sx q[2];
rz(-2.837193) q[2];
sx q[2];
rz(0.89140618) q[2];
rz(-1.8799479) q[3];
sx q[3];
rz(-1.5104048) q[3];
sx q[3];
rz(2.4752899) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0618184) q[0];
sx q[0];
rz(-0.49538716) q[0];
sx q[0];
rz(2.1451982) q[0];
rz(-2.2415316) q[1];
sx q[1];
rz(-2.3521017) q[1];
sx q[1];
rz(-0.17734227) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19781659) q[0];
sx q[0];
rz(-1.7223486) q[0];
sx q[0];
rz(-3.1413548) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5992237) q[2];
sx q[2];
rz(-2.5027983) q[2];
sx q[2];
rz(0.5109238) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23456043) q[1];
sx q[1];
rz(-1.8563358) q[1];
sx q[1];
rz(2.0689059) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0796156) q[3];
sx q[3];
rz(-1.6556228) q[3];
sx q[3];
rz(-1.2691085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29048723) q[2];
sx q[2];
rz(-1.9600211) q[2];
sx q[2];
rz(0.30430749) q[2];
rz(-0.32637063) q[3];
sx q[3];
rz(-2.5984952) q[3];
sx q[3];
rz(-2.2296026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0684763) q[0];
sx q[0];
rz(-2.731972) q[0];
sx q[0];
rz(0.9325183) q[0];
rz(0.069123507) q[1];
sx q[1];
rz(-0.2176452) q[1];
sx q[1];
rz(-2.1014012) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6118263) q[0];
sx q[0];
rz(-1.5450204) q[0];
sx q[0];
rz(2.3759205) q[0];
rz(0.17449504) q[2];
sx q[2];
rz(-1.3180079) q[2];
sx q[2];
rz(-2.6215907) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2963075) q[1];
sx q[1];
rz(-1.4060548) q[1];
sx q[1];
rz(0.77033019) q[1];
rz(-pi) q[2];
rz(-2.4916441) q[3];
sx q[3];
rz(-0.81226617) q[3];
sx q[3];
rz(2.5358729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13504623) q[2];
sx q[2];
rz(-0.73837215) q[2];
sx q[2];
rz(0.73299232) q[2];
rz(-2.0586355) q[3];
sx q[3];
rz(-2.0854009) q[3];
sx q[3];
rz(1.0068896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55423823) q[0];
sx q[0];
rz(-3.0576958) q[0];
sx q[0];
rz(2.6485637) q[0];
rz(0.21656491) q[1];
sx q[1];
rz(-2.000587) q[1];
sx q[1];
rz(-2.1176178) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7465377) q[0];
sx q[0];
rz(-0.97486541) q[0];
sx q[0];
rz(-1.5563957) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31957133) q[2];
sx q[2];
rz(-0.34863483) q[2];
sx q[2];
rz(-0.92952585) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4272642) q[1];
sx q[1];
rz(-1.9275394) q[1];
sx q[1];
rz(1.4209695) q[1];
rz(-pi) q[2];
rz(2.2129503) q[3];
sx q[3];
rz(-1.3285884) q[3];
sx q[3];
rz(0.95706576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0852802) q[2];
sx q[2];
rz(-1.6322501) q[2];
sx q[2];
rz(2.4658266) q[2];
rz(-2.7285649) q[3];
sx q[3];
rz(-1.4512117) q[3];
sx q[3];
rz(-2.5954424) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6290879) q[0];
sx q[0];
rz(-2.4830723) q[0];
sx q[0];
rz(3.107048) q[0];
rz(-0.054556219) q[1];
sx q[1];
rz(-1.9787534) q[1];
sx q[1];
rz(0.82130718) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23653447) q[0];
sx q[0];
rz(-0.58249677) q[0];
sx q[0];
rz(0.18180099) q[0];
x q[1];
rz(1.5615145) q[2];
sx q[2];
rz(-1.6595652) q[2];
sx q[2];
rz(-0.87997681) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9266204) q[1];
sx q[1];
rz(-1.6865623) q[1];
sx q[1];
rz(-2.9909219) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91854878) q[3];
sx q[3];
rz(-1.1763185) q[3];
sx q[3];
rz(-0.6834417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7811232) q[2];
sx q[2];
rz(-1.0162153) q[2];
sx q[2];
rz(-0.13295573) q[2];
rz(2.7305056) q[3];
sx q[3];
rz(-1.6229595) q[3];
sx q[3];
rz(2.0558004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046831176) q[0];
sx q[0];
rz(-0.60698858) q[0];
sx q[0];
rz(-1.2588311) q[0];
rz(1.5065441) q[1];
sx q[1];
rz(-0.656744) q[1];
sx q[1];
rz(-0.81319317) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5933605) q[0];
sx q[0];
rz(-1.5357067) q[0];
sx q[0];
rz(0.078087383) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4237464) q[2];
sx q[2];
rz(-2.4644901) q[2];
sx q[2];
rz(2.8092217) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5845994) q[1];
sx q[1];
rz(-2.0210588) q[1];
sx q[1];
rz(-2.492639) q[1];
rz(2.8714259) q[3];
sx q[3];
rz(-2.674683) q[3];
sx q[3];
rz(-1.8570569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64858156) q[2];
sx q[2];
rz(-1.1375256) q[2];
sx q[2];
rz(1.0515593) q[2];
rz(-2.1227396) q[3];
sx q[3];
rz(-2.2994883) q[3];
sx q[3];
rz(2.4837608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1558341) q[0];
sx q[0];
rz(-2.4875165) q[0];
sx q[0];
rz(-2.2610337) q[0];
rz(-1.3195994) q[1];
sx q[1];
rz(-1.9965912) q[1];
sx q[1];
rz(0.051864787) q[1];
rz(-2.8116799) q[2];
sx q[2];
rz(-0.22116667) q[2];
sx q[2];
rz(0.55851182) q[2];
rz(1.5108091) q[3];
sx q[3];
rz(-1.7675478) q[3];
sx q[3];
rz(0.02789733) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
