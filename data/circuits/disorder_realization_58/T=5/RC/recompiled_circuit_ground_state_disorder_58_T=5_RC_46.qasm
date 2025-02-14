OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4492884) q[0];
sx q[0];
rz(3.6078499) q[0];
sx q[0];
rz(10.719263) q[0];
rz(2.8757088) q[1];
sx q[1];
rz(-0.20800132) q[1];
sx q[1];
rz(1.8880358) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82845518) q[0];
sx q[0];
rz(-1.0424409) q[0];
sx q[0];
rz(1.7067451) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0090583) q[2];
sx q[2];
rz(-1.1992559) q[2];
sx q[2];
rz(-1.3309935) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80501295) q[1];
sx q[1];
rz(-2.4396371) q[1];
sx q[1];
rz(1.9483267) q[1];
rz(-pi) q[2];
rz(3.1224566) q[3];
sx q[3];
rz(-1.8950671) q[3];
sx q[3];
rz(-1.1391885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1050538) q[2];
sx q[2];
rz(-1.4194856) q[2];
sx q[2];
rz(-1.6671906) q[2];
rz(-0.27753943) q[3];
sx q[3];
rz(-1.2801291) q[3];
sx q[3];
rz(2.2138219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14811806) q[0];
sx q[0];
rz(-0.19522218) q[0];
sx q[0];
rz(-1.8147722) q[0];
rz(0.38816342) q[1];
sx q[1];
rz(-1.5048985) q[1];
sx q[1];
rz(2.6249552) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1735794) q[0];
sx q[0];
rz(-1.359419) q[0];
sx q[0];
rz(-1.4175936) q[0];
x q[1];
rz(-1.1519127) q[2];
sx q[2];
rz(-1.7563213) q[2];
sx q[2];
rz(-2.1806661) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0320466) q[1];
sx q[1];
rz(-1.418772) q[1];
sx q[1];
rz(-0.8117453) q[1];
rz(-1.850892) q[3];
sx q[3];
rz(-1.8098272) q[3];
sx q[3];
rz(1.7861507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8119729) q[2];
sx q[2];
rz(-0.5054349) q[2];
sx q[2];
rz(2.2774515) q[2];
rz(-0.69178528) q[3];
sx q[3];
rz(-2.0103879) q[3];
sx q[3];
rz(-2.2085021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.316204) q[0];
sx q[0];
rz(-0.80556691) q[0];
sx q[0];
rz(-1.4027931) q[0];
rz(-0.56651506) q[1];
sx q[1];
rz(-1.5048051) q[1];
sx q[1];
rz(-1.8623955) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7544781) q[0];
sx q[0];
rz(-2.8818948) q[0];
sx q[0];
rz(1.5952871) q[0];
x q[1];
rz(-2.6635936) q[2];
sx q[2];
rz(-2.1076116) q[2];
sx q[2];
rz(1.2144685) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0356939) q[1];
sx q[1];
rz(-2.0146684) q[1];
sx q[1];
rz(-3.1015293) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6353929) q[3];
sx q[3];
rz(-0.69033128) q[3];
sx q[3];
rz(-1.3439646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0684356) q[2];
sx q[2];
rz(-0.34978875) q[2];
sx q[2];
rz(1.5514099) q[2];
rz(-2.6326211) q[3];
sx q[3];
rz(-0.83566982) q[3];
sx q[3];
rz(1.6905748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1125672) q[0];
sx q[0];
rz(-1.644716) q[0];
sx q[0];
rz(0.56370869) q[0];
rz(1.6911223) q[1];
sx q[1];
rz(-1.9861954) q[1];
sx q[1];
rz(-0.82121003) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7477363) q[0];
sx q[0];
rz(-0.85490037) q[0];
sx q[0];
rz(-1.5564043) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8217373) q[2];
sx q[2];
rz(-1.8531688) q[2];
sx q[2];
rz(-0.47114633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5933696) q[1];
sx q[1];
rz(-1.0564305) q[1];
sx q[1];
rz(0.84131188) q[1];
rz(1.0276919) q[3];
sx q[3];
rz(-1.6031577) q[3];
sx q[3];
rz(-1.1278271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44068367) q[2];
sx q[2];
rz(-2.4606885) q[2];
sx q[2];
rz(1.7636501) q[2];
rz(2.9772229) q[3];
sx q[3];
rz(-1.5458958) q[3];
sx q[3];
rz(0.36851287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.565777) q[0];
sx q[0];
rz(-1.6933279) q[0];
sx q[0];
rz(-2.6337295) q[0];
rz(2.2841618) q[1];
sx q[1];
rz(-1.6897759) q[1];
sx q[1];
rz(-2.6618777) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3025359) q[0];
sx q[0];
rz(-0.73213644) q[0];
sx q[0];
rz(-0.49473543) q[0];
rz(-2.7577997) q[2];
sx q[2];
rz(-1.6365882) q[2];
sx q[2];
rz(-3.095997) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.060412571) q[1];
sx q[1];
rz(-2.4890889) q[1];
sx q[1];
rz(-1.1954225) q[1];
rz(-1.671319) q[3];
sx q[3];
rz(-0.28083235) q[3];
sx q[3];
rz(-0.36666825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5053284) q[2];
sx q[2];
rz(-1.1050861) q[2];
sx q[2];
rz(-0.21978933) q[2];
rz(2.9124741) q[3];
sx q[3];
rz(-0.90355211) q[3];
sx q[3];
rz(0.98178274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55766469) q[0];
sx q[0];
rz(-2.5983577) q[0];
sx q[0];
rz(-1.1274717) q[0];
rz(-0.30578956) q[1];
sx q[1];
rz(-1.1916279) q[1];
sx q[1];
rz(-2.9514899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0677493) q[0];
sx q[0];
rz(-1.1007358) q[0];
sx q[0];
rz(0.36202927) q[0];
rz(-pi) q[1];
rz(0.37895149) q[2];
sx q[2];
rz(-2.3109461) q[2];
sx q[2];
rz(-1.3791858) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.44675436) q[1];
sx q[1];
rz(-0.86320048) q[1];
sx q[1];
rz(-2.7622499) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6473057) q[3];
sx q[3];
rz(-1.8222295) q[3];
sx q[3];
rz(0.86800324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.019913435) q[2];
sx q[2];
rz(-1.6807669) q[2];
sx q[2];
rz(-1.100568) q[2];
rz(-1.3044283) q[3];
sx q[3];
rz(-1.5893987) q[3];
sx q[3];
rz(1.7884375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8562427) q[0];
sx q[0];
rz(-0.42735639) q[0];
sx q[0];
rz(-2.5653978) q[0];
rz(1.0528437) q[1];
sx q[1];
rz(-2.5710227) q[1];
sx q[1];
rz(-1.0594692) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5813025) q[0];
sx q[0];
rz(-1.6972739) q[0];
sx q[0];
rz(2.4471388) q[0];
x q[1];
rz(-1.1722819) q[2];
sx q[2];
rz(-1.544612) q[2];
sx q[2];
rz(0.46253935) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1394315) q[1];
sx q[1];
rz(-2.0116848) q[1];
sx q[1];
rz(1.3590165) q[1];
x q[2];
rz(-0.7110157) q[3];
sx q[3];
rz(-1.064015) q[3];
sx q[3];
rz(-1.8975951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8656371) q[2];
sx q[2];
rz(-0.58791462) q[2];
sx q[2];
rz(2.0666583) q[2];
rz(2.2771207) q[3];
sx q[3];
rz(-1.8755707) q[3];
sx q[3];
rz(1.5629432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7445755) q[0];
sx q[0];
rz(-0.96857849) q[0];
sx q[0];
rz(2.6343935) q[0];
rz(1.1811258) q[1];
sx q[1];
rz(-1.487251) q[1];
sx q[1];
rz(-2.2307253) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1483123) q[0];
sx q[0];
rz(-2.1851843) q[0];
sx q[0];
rz(-2.5117158) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1520432) q[2];
sx q[2];
rz(-1.74161) q[2];
sx q[2];
rz(0.51994158) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4994608) q[1];
sx q[1];
rz(-1.5822463) q[1];
sx q[1];
rz(1.7246555) q[1];
x q[2];
rz(2.0924545) q[3];
sx q[3];
rz(-0.46346617) q[3];
sx q[3];
rz(1.3090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7562423) q[2];
sx q[2];
rz(-1.664398) q[2];
sx q[2];
rz(0.64794668) q[2];
rz(-0.09662763) q[3];
sx q[3];
rz(-2.1164618) q[3];
sx q[3];
rz(-0.9534165) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18428093) q[0];
sx q[0];
rz(-0.98783699) q[0];
sx q[0];
rz(1.8010358) q[0];
rz(0.52349177) q[1];
sx q[1];
rz(-1.5308056) q[1];
sx q[1];
rz(1.2045822) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0055375) q[0];
sx q[0];
rz(-0.82804543) q[0];
sx q[0];
rz(2.8715747) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9926489) q[2];
sx q[2];
rz(-1.3729549) q[2];
sx q[2];
rz(2.1959675) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.5960418) q[1];
sx q[1];
rz(-1.4117472) q[1];
sx q[1];
rz(-0.56175128) q[1];
rz(-pi) q[2];
rz(2.0818578) q[3];
sx q[3];
rz(-1.9283224) q[3];
sx q[3];
rz(-0.75393576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.039310731) q[2];
sx q[2];
rz(-1.0264341) q[2];
sx q[2];
rz(0.2956051) q[2];
rz(-3.0292656) q[3];
sx q[3];
rz(-1.5528409) q[3];
sx q[3];
rz(-2.2789392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2304147) q[0];
sx q[0];
rz(-2.6417612) q[0];
sx q[0];
rz(-2.3590132) q[0];
rz(0.29620194) q[1];
sx q[1];
rz(-1.240851) q[1];
sx q[1];
rz(2.9414419) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81129942) q[0];
sx q[0];
rz(-1.4084554) q[0];
sx q[0];
rz(2.3672124) q[0];
rz(-pi) q[1];
rz(-1.0050943) q[2];
sx q[2];
rz(-1.8100693) q[2];
sx q[2];
rz(-2.4636961) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0072277) q[1];
sx q[1];
rz(-1.0505465) q[1];
sx q[1];
rz(-0.84905973) q[1];
x q[2];
rz(-2.5298821) q[3];
sx q[3];
rz(-2.1107622) q[3];
sx q[3];
rz(-1.5276791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.062181648) q[2];
sx q[2];
rz(-0.14180413) q[2];
sx q[2];
rz(1.0395435) q[2];
rz(3.1145575) q[3];
sx q[3];
rz(-0.98374933) q[3];
sx q[3];
rz(-0.99193096) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1886002) q[0];
sx q[0];
rz(-0.66467265) q[0];
sx q[0];
rz(-1.959214) q[0];
rz(1.7091119) q[1];
sx q[1];
rz(-1.8791589) q[1];
sx q[1];
rz(-1.5502677) q[1];
rz(-0.15122945) q[2];
sx q[2];
rz(-1.4561903) q[2];
sx q[2];
rz(-1.0753808) q[2];
rz(0.19743528) q[3];
sx q[3];
rz(-2.6876269) q[3];
sx q[3];
rz(-0.85314565) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
