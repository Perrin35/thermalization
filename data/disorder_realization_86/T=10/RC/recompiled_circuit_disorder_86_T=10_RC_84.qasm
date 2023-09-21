OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(-2.2280333) q[0];
sx q[0];
rz(1.7295184) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(5.6875416) q[1];
sx q[1];
rz(7.7653801) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55573758) q[0];
sx q[0];
rz(-1.1475539) q[0];
sx q[0];
rz(-1.7413571) q[0];
x q[1];
rz(1.4380768) q[2];
sx q[2];
rz(-1.8258397) q[2];
sx q[2];
rz(2.1357352) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9285779) q[1];
sx q[1];
rz(-1.9933812) q[1];
sx q[1];
rz(-1.0282474) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0406038) q[3];
sx q[3];
rz(-2.1312993) q[3];
sx q[3];
rz(1.3274173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1564864) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(0.86581725) q[2];
rz(2.1872897) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(-1.8538063) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99825478) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(-3.1153733) q[0];
rz(-1.6014618) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(0.96347934) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.11461) q[0];
sx q[0];
rz(-2.5275505) q[0];
sx q[0];
rz(-1.5747889) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0785525) q[2];
sx q[2];
rz(-2.5139367) q[2];
sx q[2];
rz(2.9693659) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.023167921) q[1];
sx q[1];
rz(-1.2512565) q[1];
sx q[1];
rz(2.097514) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2198592) q[3];
sx q[3];
rz(-1.8027455) q[3];
sx q[3];
rz(-1.2082781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5144689) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(-0.13452402) q[2];
rz(0.7450122) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(-2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2117675) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(-0.79743687) q[0];
rz(-2.0939317) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(-0.55999666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6277498) q[0];
sx q[0];
rz(-1.7729513) q[0];
sx q[0];
rz(-1.9613128) q[0];
x q[1];
rz(-2.838344) q[2];
sx q[2];
rz(-1.5911284) q[2];
sx q[2];
rz(-0.081239935) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3126038) q[1];
sx q[1];
rz(-1.7692411) q[1];
sx q[1];
rz(2.7752084) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79389823) q[3];
sx q[3];
rz(-0.63459914) q[3];
sx q[3];
rz(-1.3018228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.75227633) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(2.9690572) q[2];
rz(2.1595188) q[3];
sx q[3];
rz(-1.7445824) q[3];
sx q[3];
rz(1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.003222) q[0];
sx q[0];
rz(-2.0439742) q[0];
sx q[0];
rz(0.28451434) q[0];
rz(0.31670397) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(-1.8428615) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239149) q[0];
sx q[0];
rz(-2.0664584) q[0];
sx q[0];
rz(-0.25650521) q[0];
rz(0.0012871731) q[2];
sx q[2];
rz(-0.80438559) q[2];
sx q[2];
rz(3.121701) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5609834) q[1];
sx q[1];
rz(-2.1316075) q[1];
sx q[1];
rz(-0.22052712) q[1];
x q[2];
rz(-0.11233791) q[3];
sx q[3];
rz(-1.2445407) q[3];
sx q[3];
rz(2.5374967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(0.38468012) q[2];
rz(2.3875333) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(-1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5383179) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(1.7549365) q[0];
rz(2.9105913) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(-2.8447661) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1332069) q[0];
sx q[0];
rz(-1.0002245) q[0];
sx q[0];
rz(1.617336) q[0];
rz(-pi) q[1];
rz(-1.417824) q[2];
sx q[2];
rz(-0.83380552) q[2];
sx q[2];
rz(-1.0066777) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0592812) q[1];
sx q[1];
rz(-1.1059522) q[1];
sx q[1];
rz(0.44537284) q[1];
x q[2];
rz(-1.6597219) q[3];
sx q[3];
rz(-2.2825135) q[3];
sx q[3];
rz(0.54642788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(0.39247593) q[2];
rz(1.1522419) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(-2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0157938) q[0];
sx q[0];
rz(-1.5690465) q[0];
sx q[0];
rz(2.3902067) q[0];
rz(1.3279351) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(-2.5352535) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5099112) q[0];
sx q[0];
rz(-1.5874377) q[0];
sx q[0];
rz(-3.095754) q[0];
rz(-pi) q[1];
rz(-2.2491127) q[2];
sx q[2];
rz(-1.2391029) q[2];
sx q[2];
rz(-1.734037) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50748435) q[1];
sx q[1];
rz(-2.1180696) q[1];
sx q[1];
rz(0.21168153) q[1];
rz(-pi) q[2];
rz(1.8984406) q[3];
sx q[3];
rz(-2.9122105) q[3];
sx q[3];
rz(-2.4250507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(2.0022557) q[2];
rz(-1.4849439) q[3];
sx q[3];
rz(-1.1805725) q[3];
sx q[3];
rz(0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5320324) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(-0.50810057) q[0];
rz(1.5628901) q[1];
sx q[1];
rz(-1.0888313) q[1];
sx q[1];
rz(-0.79024822) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9433141) q[0];
sx q[0];
rz(-2.4662848) q[0];
sx q[0];
rz(0.4839464) q[0];
rz(-pi) q[1];
rz(-3.087567) q[2];
sx q[2];
rz(-2.9276491) q[2];
sx q[2];
rz(-0.33188785) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3527457) q[1];
sx q[1];
rz(-0.48301304) q[1];
sx q[1];
rz(-2.5301945) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1931476) q[3];
sx q[3];
rz(-0.95818633) q[3];
sx q[3];
rz(2.2539504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0885075) q[2];
sx q[2];
rz(-2.6997456) q[2];
sx q[2];
rz(1.7283758) q[2];
rz(-1.6648071) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(-3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168032) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(1.0472263) q[0];
rz(-0.60910243) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0213288) q[0];
sx q[0];
rz(-1.6006032) q[0];
sx q[0];
rz(2.1304312) q[0];
rz(-1.1461166) q[2];
sx q[2];
rz(-2.0600024) q[2];
sx q[2];
rz(-0.81791544) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.70506239) q[1];
sx q[1];
rz(-0.092518004) q[1];
sx q[1];
rz(-3.0180879) q[1];
x q[2];
rz(-1.9854529) q[3];
sx q[3];
rz(-2.4932043) q[3];
sx q[3];
rz(-0.63601953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9528815) q[2];
sx q[2];
rz(-2.7313576) q[2];
sx q[2];
rz(-2.2593373) q[2];
rz(1.7404209) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.9410979) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27154487) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(1.7154988) q[0];
rz(3.0601314) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(-0.55823278) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0575858) q[0];
sx q[0];
rz(-1.1450197) q[0];
sx q[0];
rz(2.7140679) q[0];
rz(-pi) q[1];
rz(-1.105537) q[2];
sx q[2];
rz(-0.12013398) q[2];
sx q[2];
rz(0.58819729) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3840752) q[1];
sx q[1];
rz(-1.6346524) q[1];
sx q[1];
rz(-0.82660316) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6551412) q[3];
sx q[3];
rz(-1.7351741) q[3];
sx q[3];
rz(2.1742976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2925064) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(2.9294087) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(1.2020483) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6367209) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(-1.5378392) q[0];
rz(-2.3161855) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(-0.5232946) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6494203) q[0];
sx q[0];
rz(-1.5086552) q[0];
sx q[0];
rz(0.6092351) q[0];
x q[1];
rz(-2.142749) q[2];
sx q[2];
rz(-1.5339282) q[2];
sx q[2];
rz(1.8429304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2460675) q[1];
sx q[1];
rz(-1.9290036) q[1];
sx q[1];
rz(1.2342865) q[1];
rz(1.80199) q[3];
sx q[3];
rz(-1.3320859) q[3];
sx q[3];
rz(-1.064144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.837073) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(0.26930299) q[2];
rz(0.4942016) q[3];
sx q[3];
rz(-0.84635693) q[3];
sx q[3];
rz(-1.0860898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3257278) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(-1.6745463) q[1];
sx q[1];
rz(-2.84927) q[1];
sx q[1];
rz(-1.897859) q[1];
rz(0.011209839) q[2];
sx q[2];
rz(-1.8335473) q[2];
sx q[2];
rz(2.0469472) q[2];
rz(0.23444093) q[3];
sx q[3];
rz(-0.81080484) q[3];
sx q[3];
rz(-0.57639359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];