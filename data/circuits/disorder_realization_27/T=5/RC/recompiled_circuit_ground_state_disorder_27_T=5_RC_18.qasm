OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34484997) q[0];
sx q[0];
rz(-0.27422187) q[0];
sx q[0];
rz(-2.5728777) q[0];
rz(1.2110127) q[1];
sx q[1];
rz(-2.14415) q[1];
sx q[1];
rz(2.8740191) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8519605) q[0];
sx q[0];
rz(-1.3697764) q[0];
sx q[0];
rz(0.55390771) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2368343) q[2];
sx q[2];
rz(-1.4638454) q[2];
sx q[2];
rz(0.013205139) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1700344) q[1];
sx q[1];
rz(-1.5791248) q[1];
sx q[1];
rz(-0.002408601) q[1];
rz(-pi) q[2];
x q[2];
rz(0.062248793) q[3];
sx q[3];
rz(-1.4302974) q[3];
sx q[3];
rz(1.922091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25207511) q[2];
sx q[2];
rz(-1.0675665) q[2];
sx q[2];
rz(-1.1313103) q[2];
rz(-2.6913397) q[3];
sx q[3];
rz(-2.4501652) q[3];
sx q[3];
rz(0.94436193) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4985519) q[0];
sx q[0];
rz(-2.2096071) q[0];
sx q[0];
rz(-1.407628) q[0];
rz(-0.44250008) q[1];
sx q[1];
rz(-1.718037) q[1];
sx q[1];
rz(0.59534591) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.567826) q[0];
sx q[0];
rz(-1.8086595) q[0];
sx q[0];
rz(0.25517558) q[0];
x q[1];
rz(-0.96244241) q[2];
sx q[2];
rz(-2.3909937) q[2];
sx q[2];
rz(-0.19718328) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9756979) q[1];
sx q[1];
rz(-1.5656316) q[1];
sx q[1];
rz(2.7126606) q[1];
x q[2];
rz(-0.54259681) q[3];
sx q[3];
rz(-0.68477453) q[3];
sx q[3];
rz(2.7921576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.88196102) q[2];
sx q[2];
rz(-0.61820784) q[2];
sx q[2];
rz(2.372443) q[2];
rz(1.7800356) q[3];
sx q[3];
rz(-2.1741314) q[3];
sx q[3];
rz(-2.5462525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24460569) q[0];
sx q[0];
rz(-0.91490442) q[0];
sx q[0];
rz(-2.8047674) q[0];
rz(1.7103851) q[1];
sx q[1];
rz(-2.2957048) q[1];
sx q[1];
rz(3.0444042) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6568156) q[0];
sx q[0];
rz(-2.6069399) q[0];
sx q[0];
rz(-1.9770245) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4732185) q[2];
sx q[2];
rz(-1.1910025) q[2];
sx q[2];
rz(-1.9623836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0984368) q[1];
sx q[1];
rz(-1.7088026) q[1];
sx q[1];
rz(-0.34855493) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8840553) q[3];
sx q[3];
rz(-1.9935605) q[3];
sx q[3];
rz(0.54856578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.97239697) q[2];
sx q[2];
rz(-1.376386) q[2];
sx q[2];
rz(-0.81673679) q[2];
rz(2.2190602) q[3];
sx q[3];
rz(-0.43729344) q[3];
sx q[3];
rz(1.9492662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794491) q[0];
sx q[0];
rz(-1.3379931) q[0];
sx q[0];
rz(2.2747967) q[0];
rz(-1.9056994) q[1];
sx q[1];
rz(-1.1150603) q[1];
sx q[1];
rz(2.8996276) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2529255) q[0];
sx q[0];
rz(-0.87742108) q[0];
sx q[0];
rz(2.7916551) q[0];
rz(-pi) q[1];
rz(1.1461805) q[2];
sx q[2];
rz(-0.17689366) q[2];
sx q[2];
rz(-1.5695656) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3212525) q[1];
sx q[1];
rz(-0.42916052) q[1];
sx q[1];
rz(2.8080775) q[1];
x q[2];
rz(0.49698835) q[3];
sx q[3];
rz(-1.6145633) q[3];
sx q[3];
rz(2.5183293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8786826) q[2];
sx q[2];
rz(-1.6308558) q[2];
sx q[2];
rz(-0.53544694) q[2];
rz(2.6483436) q[3];
sx q[3];
rz(-2.2507164) q[3];
sx q[3];
rz(2.6935327) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066147476) q[0];
sx q[0];
rz(-2.9904521) q[0];
sx q[0];
rz(-0.84392631) q[0];
rz(-1.4752202) q[1];
sx q[1];
rz(-1.4862783) q[1];
sx q[1];
rz(-0.39316887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61680865) q[0];
sx q[0];
rz(-2.3032585) q[0];
sx q[0];
rz(-0.34466593) q[0];
x q[1];
rz(-1.3530144) q[2];
sx q[2];
rz(-0.94557191) q[2];
sx q[2];
rz(-0.54131258) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.95994334) q[1];
sx q[1];
rz(-1.8487329) q[1];
sx q[1];
rz(-2.6647869) q[1];
rz(-pi) q[2];
rz(0.78759463) q[3];
sx q[3];
rz(-1.9301842) q[3];
sx q[3];
rz(-1.5415292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4142485) q[2];
sx q[2];
rz(-0.20432893) q[2];
sx q[2];
rz(2.7434529) q[2];
rz(-2.5967755) q[3];
sx q[3];
rz(-0.78854338) q[3];
sx q[3];
rz(-1.8383693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9141465) q[0];
sx q[0];
rz(-2.1307724) q[0];
sx q[0];
rz(-2.2858802) q[0];
rz(0.69264597) q[1];
sx q[1];
rz(-2.1470862) q[1];
sx q[1];
rz(2.8725502) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5554898) q[0];
sx q[0];
rz(-1.5464096) q[0];
sx q[0];
rz(1.3792319) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0641818) q[2];
sx q[2];
rz(-1.6010188) q[2];
sx q[2];
rz(-0.53080785) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66722357) q[1];
sx q[1];
rz(-1.3880236) q[1];
sx q[1];
rz(2.3007042) q[1];
x q[2];
rz(1.7625436) q[3];
sx q[3];
rz(-1.7320398) q[3];
sx q[3];
rz(-1.2210326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1143703) q[2];
sx q[2];
rz(-1.8222858) q[2];
sx q[2];
rz(-0.11030062) q[2];
rz(2.2731764) q[3];
sx q[3];
rz(-1.1766368) q[3];
sx q[3];
rz(-1.4364012) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7480943) q[0];
sx q[0];
rz(-0.49999923) q[0];
sx q[0];
rz(2.2802343) q[0];
rz(1.6294847) q[1];
sx q[1];
rz(-1.6810828) q[1];
sx q[1];
rz(-0.82383627) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.829946) q[0];
sx q[0];
rz(-2.5426546) q[0];
sx q[0];
rz(1.3715368) q[0];
rz(-2.3917213) q[2];
sx q[2];
rz(-1.2629384) q[2];
sx q[2];
rz(1.6689545) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7303402) q[1];
sx q[1];
rz(-1.3578051) q[1];
sx q[1];
rz(0.61720444) q[1];
rz(-pi) q[2];
rz(-2.7384375) q[3];
sx q[3];
rz(-2.4488827) q[3];
sx q[3];
rz(1.5628536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6350101) q[2];
sx q[2];
rz(-0.51598769) q[2];
sx q[2];
rz(-0.79279509) q[2];
rz(-0.005216287) q[3];
sx q[3];
rz(-2.3514533) q[3];
sx q[3];
rz(-1.4504356) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.97013) q[0];
sx q[0];
rz(-1.9487533) q[0];
sx q[0];
rz(0.18950732) q[0];
rz(2.3686523) q[1];
sx q[1];
rz(-2.645292) q[1];
sx q[1];
rz(0.596284) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8533289) q[0];
sx q[0];
rz(-1.402308) q[0];
sx q[0];
rz(-0.69126076) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79769602) q[2];
sx q[2];
rz(-2.0828649) q[2];
sx q[2];
rz(-1.2810117) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.50508037) q[1];
sx q[1];
rz(-0.13992913) q[1];
sx q[1];
rz(-0.72461463) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1607481) q[3];
sx q[3];
rz(-0.80435565) q[3];
sx q[3];
rz(0.70589069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72558713) q[2];
sx q[2];
rz(-3.1276939) q[2];
sx q[2];
rz(1.928398) q[2];
rz(0.98617918) q[3];
sx q[3];
rz(-1.4158019) q[3];
sx q[3];
rz(-1.2590316) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1213433) q[0];
sx q[0];
rz(-1.5859402) q[0];
sx q[0];
rz(1.1100618) q[0];
rz(-2.8129261) q[1];
sx q[1];
rz(-1.5549436) q[1];
sx q[1];
rz(1.2967671) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.219562) q[0];
sx q[0];
rz(-2.2264105) q[0];
sx q[0];
rz(-3.0232885) q[0];
rz(-3.0795394) q[2];
sx q[2];
rz(-1.9665355) q[2];
sx q[2];
rz(1.526265) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4420524) q[1];
sx q[1];
rz(-1.6332383) q[1];
sx q[1];
rz(0.49380912) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27196692) q[3];
sx q[3];
rz(-0.68373954) q[3];
sx q[3];
rz(-2.7955987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0729596) q[2];
sx q[2];
rz(-2.1559842) q[2];
sx q[2];
rz(0.10406058) q[2];
rz(2.0067298) q[3];
sx q[3];
rz(-1.3645423) q[3];
sx q[3];
rz(-2.9141736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65570152) q[0];
sx q[0];
rz(-2.1567397) q[0];
sx q[0];
rz(0.73053288) q[0];
rz(2.5841374) q[1];
sx q[1];
rz(-1.1703706) q[1];
sx q[1];
rz(2.7117859) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.902122) q[0];
sx q[0];
rz(-0.80901635) q[0];
sx q[0];
rz(1.5714963) q[0];
rz(-pi) q[1];
rz(-2.5556106) q[2];
sx q[2];
rz(-0.44011099) q[2];
sx q[2];
rz(-2.7173079) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.7260202) q[1];
sx q[1];
rz(-1.4910798) q[1];
sx q[1];
rz(1.1989051) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9409995) q[3];
sx q[3];
rz(-1.4552081) q[3];
sx q[3];
rz(-0.59578958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3320015) q[2];
sx q[2];
rz(-2.1198699) q[2];
sx q[2];
rz(-0.29279718) q[2];
rz(3.0012567) q[3];
sx q[3];
rz(-2.1122746) q[3];
sx q[3];
rz(2.6405507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3705227) q[0];
sx q[0];
rz(-1.6765544) q[0];
sx q[0];
rz(0.21677207) q[0];
rz(0.71939214) q[1];
sx q[1];
rz(-1.2750625) q[1];
sx q[1];
rz(0.19663179) q[1];
rz(-1.8875296) q[2];
sx q[2];
rz(-0.52957305) q[2];
sx q[2];
rz(-1.1323462) q[2];
rz(-1.1872798) q[3];
sx q[3];
rz(-2.7957932) q[3];
sx q[3];
rz(2.3898706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
