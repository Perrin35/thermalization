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
rz(-1.93058) q[1];
sx q[1];
rz(-0.99744263) q[1];
sx q[1];
rz(-2.8740191) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28963213) q[0];
sx q[0];
rz(-1.3697764) q[0];
sx q[0];
rz(-0.55390771) q[0];
x q[1];
rz(1.2368343) q[2];
sx q[2];
rz(-1.6777473) q[2];
sx q[2];
rz(-0.013205139) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5423746) q[1];
sx q[1];
rz(-1.5732048) q[1];
sx q[1];
rz(1.5791248) q[1];
rz(-pi) q[2];
rz(3.0793439) q[3];
sx q[3];
rz(-1.7112952) q[3];
sx q[3];
rz(1.922091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8895175) q[2];
sx q[2];
rz(-1.0675665) q[2];
sx q[2];
rz(2.0102823) q[2];
rz(-2.6913397) q[3];
sx q[3];
rz(-2.4501652) q[3];
sx q[3];
rz(0.94436193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4985519) q[0];
sx q[0];
rz(-2.2096071) q[0];
sx q[0];
rz(-1.407628) q[0];
rz(-2.6990926) q[1];
sx q[1];
rz(-1.4235556) q[1];
sx q[1];
rz(-2.5462467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.731643) q[0];
sx q[0];
rz(-0.34706693) q[0];
sx q[0];
rz(-0.76526977) q[0];
rz(0.91752618) q[2];
sx q[2];
rz(-1.9712312) q[2];
sx q[2];
rz(1.296907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1658948) q[1];
sx q[1];
rz(-1.575961) q[1];
sx q[1];
rz(-2.7126606) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61025783) q[3];
sx q[3];
rz(-1.2380945) q[3];
sx q[3];
rz(-1.4833029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.88196102) q[2];
sx q[2];
rz(-0.61820784) q[2];
sx q[2];
rz(-2.372443) q[2];
rz(1.7800356) q[3];
sx q[3];
rz(-0.96746126) q[3];
sx q[3];
rz(-0.59534016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.896987) q[0];
sx q[0];
rz(-0.91490442) q[0];
sx q[0];
rz(-0.33682522) q[0];
rz(1.7103851) q[1];
sx q[1];
rz(-0.84588784) q[1];
sx q[1];
rz(-3.0444042) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1203994) q[0];
sx q[0];
rz(-1.0836856) q[0];
sx q[0];
rz(2.9117286) q[0];
rz(0.38143845) q[2];
sx q[2];
rz(-1.4801916) q[2];
sx q[2];
rz(-0.355313) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0984368) q[1];
sx q[1];
rz(-1.43279) q[1];
sx q[1];
rz(0.34855493) q[1];
rz(-pi) q[2];
rz(0.52690701) q[3];
sx q[3];
rz(-2.1873173) q[3];
sx q[3];
rz(1.7948727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1691957) q[2];
sx q[2];
rz(-1.7652067) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794491) q[0];
sx q[0];
rz(-1.3379931) q[0];
sx q[0];
rz(-0.86679593) q[0];
rz(1.9056994) q[1];
sx q[1];
rz(-1.1150603) q[1];
sx q[1];
rz(0.24196504) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2529255) q[0];
sx q[0];
rz(-0.87742108) q[0];
sx q[0];
rz(0.34993751) q[0];
x q[1];
rz(-0.073512065) q[2];
sx q[2];
rz(-1.7318372) q[2];
sx q[2];
rz(2.0025776) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82034012) q[1];
sx q[1];
rz(-2.7124321) q[1];
sx q[1];
rz(-2.8080775) q[1];
rz(-pi) q[2];
rz(2.6446043) q[3];
sx q[3];
rz(-1.5270294) q[3];
sx q[3];
rz(2.5183293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8786826) q[2];
sx q[2];
rz(-1.5107369) q[2];
sx q[2];
rz(2.6061457) q[2];
rz(-0.49324909) q[3];
sx q[3];
rz(-2.2507164) q[3];
sx q[3];
rz(-0.44805995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0754452) q[0];
sx q[0];
rz(-0.15114052) q[0];
sx q[0];
rz(-2.2976663) q[0];
rz(1.4752202) q[1];
sx q[1];
rz(-1.4862783) q[1];
sx q[1];
rz(0.39316887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0174781) q[0];
sx q[0];
rz(-0.79567608) q[0];
sx q[0];
rz(1.2114197) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3530144) q[2];
sx q[2];
rz(-0.94557191) q[2];
sx q[2];
rz(-0.54131258) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1816493) q[1];
sx q[1];
rz(-1.2928597) q[1];
sx q[1];
rz(0.47680579) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48748555) q[3];
sx q[3];
rz(-2.2922487) q[3];
sx q[3];
rz(2.8340428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4142485) q[2];
sx q[2];
rz(-0.20432893) q[2];
sx q[2];
rz(0.3981398) q[2];
rz(-2.5967755) q[3];
sx q[3];
rz(-2.3530493) q[3];
sx q[3];
rz(1.8383693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9141465) q[0];
sx q[0];
rz(-2.1307724) q[0];
sx q[0];
rz(0.8557125) q[0];
rz(0.69264597) q[1];
sx q[1];
rz(-0.99450642) q[1];
sx q[1];
rz(-2.8725502) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1097668) q[0];
sx q[0];
rz(-2.9485011) q[0];
sx q[0];
rz(1.4433799) q[0];
x q[1];
rz(-2.0774109) q[2];
sx q[2];
rz(-1.5405739) q[2];
sx q[2];
rz(-0.53080785) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0768041) q[1];
sx q[1];
rz(-0.85569438) q[1];
sx q[1];
rz(-2.8984757) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86427291) q[3];
sx q[3];
rz(-0.24989299) q[3];
sx q[3];
rz(-2.8003729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1143703) q[2];
sx q[2];
rz(-1.3193069) q[2];
sx q[2];
rz(-0.11030062) q[2];
rz(-2.2731764) q[3];
sx q[3];
rz(-1.9649558) q[3];
sx q[3];
rz(1.7051914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7480943) q[0];
sx q[0];
rz(-2.6415934) q[0];
sx q[0];
rz(-0.86135832) q[0];
rz(-1.6294847) q[1];
sx q[1];
rz(-1.4605099) q[1];
sx q[1];
rz(2.3177564) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55144009) q[0];
sx q[0];
rz(-2.1562898) q[0];
sx q[0];
rz(0.13429883) q[0];
x q[1];
rz(-0.4365224) q[2];
sx q[2];
rz(-2.3425205) q[2];
sx q[2];
rz(-2.7288849) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8331404) q[1];
sx q[1];
rz(-2.1720533) q[1];
sx q[1];
rz(-1.8300301) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8855612) q[3];
sx q[3];
rz(-2.1986695) q[3];
sx q[3];
rz(1.0726269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5065826) q[2];
sx q[2];
rz(-0.51598769) q[2];
sx q[2];
rz(-0.79279509) q[2];
rz(0.005216287) q[3];
sx q[3];
rz(-0.79013932) q[3];
sx q[3];
rz(-1.4504356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.1714627) q[0];
sx q[0];
rz(-1.9487533) q[0];
sx q[0];
rz(-0.18950732) q[0];
rz(-0.7729404) q[1];
sx q[1];
rz(-2.645292) q[1];
sx q[1];
rz(0.596284) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7211821) q[0];
sx q[0];
rz(-0.89119688) q[0];
sx q[0];
rz(-1.7880938) q[0];
rz(-pi) q[1];
rz(0.66571301) q[2];
sx q[2];
rz(-2.2253198) q[2];
sx q[2];
rz(-0.15617035) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50508037) q[1];
sx q[1];
rz(-0.13992913) q[1];
sx q[1];
rz(-2.416978) q[1];
x q[2];
rz(0.98084456) q[3];
sx q[3];
rz(-0.80435565) q[3];
sx q[3];
rz(2.435702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72558713) q[2];
sx q[2];
rz(-3.1276939) q[2];
sx q[2];
rz(-1.928398) q[2];
rz(2.1554135) q[3];
sx q[3];
rz(-1.7257907) q[3];
sx q[3];
rz(-1.2590316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1213433) q[0];
sx q[0];
rz(-1.5859402) q[0];
sx q[0];
rz(-2.0315309) q[0];
rz(-2.8129261) q[1];
sx q[1];
rz(-1.5549436) q[1];
sx q[1];
rz(1.2967671) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7294819) q[0];
sx q[0];
rz(-0.66464948) q[0];
sx q[0];
rz(1.7230711) q[0];
rz(-pi) q[1];
rz(-1.4234366) q[2];
sx q[2];
rz(-2.7412716) q[2];
sx q[2];
rz(1.4555228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3039141) q[1];
sx q[1];
rz(-2.0635567) q[1];
sx q[1];
rz(-1.4999092) q[1];
rz(-2.4761183) q[3];
sx q[3];
rz(-1.4002808) q[3];
sx q[3];
rz(-2.1297034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0729596) q[2];
sx q[2];
rz(-0.98560846) q[2];
sx q[2];
rz(-3.0375321) q[2];
rz(1.1348628) q[3];
sx q[3];
rz(-1.3645423) q[3];
sx q[3];
rz(-0.22741905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65570152) q[0];
sx q[0];
rz(-0.98485297) q[0];
sx q[0];
rz(2.4110598) q[0];
rz(-0.55745521) q[1];
sx q[1];
rz(-1.1703706) q[1];
sx q[1];
rz(2.7117859) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8107501) q[0];
sx q[0];
rz(-1.5713028) q[0];
sx q[0];
rz(2.3798126) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7676959) q[2];
sx q[2];
rz(-1.8086401) q[2];
sx q[2];
rz(1.4542945) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0954206) q[1];
sx q[1];
rz(-2.7616427) q[1];
sx q[1];
rz(-1.7871961) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12390512) q[3];
sx q[3];
rz(-1.2031816) q[3];
sx q[3];
rz(2.1218561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.80959117) q[2];
sx q[2];
rz(-2.1198699) q[2];
sx q[2];
rz(-2.8487955) q[2];
rz(0.14033595) q[3];
sx q[3];
rz(-2.1122746) q[3];
sx q[3];
rz(0.50104195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.77107) q[0];
sx q[0];
rz(-1.4650383) q[0];
sx q[0];
rz(-2.9248206) q[0];
rz(-0.71939214) q[1];
sx q[1];
rz(-1.8665301) q[1];
sx q[1];
rz(-2.9449609) q[1];
rz(1.2540631) q[2];
sx q[2];
rz(-0.52957305) q[2];
sx q[2];
rz(-1.1323462) q[2];
rz(3.0075913) q[3];
sx q[3];
rz(-1.8905427) q[3];
sx q[3];
rz(-0.34656634) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
