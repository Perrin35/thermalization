OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52842927) q[0];
sx q[0];
rz(2.0818721) q[0];
sx q[0];
rz(11.835397) q[0];
rz(-1.5001186) q[1];
sx q[1];
rz(-2.1067696) q[1];
sx q[1];
rz(0.94354454) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7109414) q[0];
sx q[0];
rz(-1.5358155) q[0];
sx q[0];
rz(-3.0745688) q[0];
rz(0.65096345) q[2];
sx q[2];
rz(-1.5032094) q[2];
sx q[2];
rz(-2.0219321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4409677) q[1];
sx q[1];
rz(-1.1807627) q[1];
sx q[1];
rz(-2.7729211) q[1];
rz(-pi) q[2];
rz(2.5991873) q[3];
sx q[3];
rz(-2.7895658) q[3];
sx q[3];
rz(-2.5792518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.25508183) q[2];
sx q[2];
rz(-1.7604897) q[2];
sx q[2];
rz(-1.250766) q[2];
rz(-1.7154153) q[3];
sx q[3];
rz(-2.2255247) q[3];
sx q[3];
rz(-0.9799408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.132906) q[0];
sx q[0];
rz(-2.0693021) q[0];
sx q[0];
rz(2.5426478) q[0];
rz(1.3409746) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(-0.96639955) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4616529) q[0];
sx q[0];
rz(-2.0113809) q[0];
sx q[0];
rz(-0.64042129) q[0];
x q[1];
rz(-2.5194089) q[2];
sx q[2];
rz(-1.6112279) q[2];
sx q[2];
rz(-1.7379023) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.84345531) q[1];
sx q[1];
rz(-1.4569439) q[1];
sx q[1];
rz(1.650412) q[1];
rz(2.3923799) q[3];
sx q[3];
rz(-1.8945165) q[3];
sx q[3];
rz(-0.92326984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0559343) q[2];
sx q[2];
rz(-2.3264383) q[2];
sx q[2];
rz(2.8404964) q[2];
rz(1.1931233) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(-1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0369204) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(1.954129) q[0];
rz(-1.2359515) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(1.3175861) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8225107) q[0];
sx q[0];
rz(-1.7829478) q[0];
sx q[0];
rz(0.69885079) q[0];
x q[1];
rz(1.7682398) q[2];
sx q[2];
rz(-1.3609481) q[2];
sx q[2];
rz(0.91453493) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.31955645) q[1];
sx q[1];
rz(-1.469538) q[1];
sx q[1];
rz(2.7201369) q[1];
x q[2];
rz(2.0882323) q[3];
sx q[3];
rz(-1.468717) q[3];
sx q[3];
rz(-1.7774338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1304156) q[2];
sx q[2];
rz(-1.3935564) q[2];
sx q[2];
rz(0.2066361) q[2];
rz(-2.4335499) q[3];
sx q[3];
rz(-0.20773023) q[3];
sx q[3];
rz(-2.2487683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77465039) q[0];
sx q[0];
rz(-2.9270524) q[0];
sx q[0];
rz(0.88622093) q[0];
rz(-2.1318502) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(-1.2264235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0013989) q[0];
sx q[0];
rz(-1.5045907) q[0];
sx q[0];
rz(-2.2982236) q[0];
rz(-pi) q[1];
rz(-2.6817276) q[2];
sx q[2];
rz(-0.93732873) q[2];
sx q[2];
rz(-1.8749274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7734079) q[1];
sx q[1];
rz(-0.28973026) q[1];
sx q[1];
rz(0.75399953) q[1];
x q[2];
rz(2.5049582) q[3];
sx q[3];
rz(-2.5599179) q[3];
sx q[3];
rz(-0.10891529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4975138) q[2];
sx q[2];
rz(-1.3138694) q[2];
sx q[2];
rz(2.148596) q[2];
rz(1.8289061) q[3];
sx q[3];
rz(-2.1195181) q[3];
sx q[3];
rz(-0.48373568) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-2.826773) q[0];
sx q[0];
rz(1.1859878) q[0];
rz(1.4233308) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(2.5591154) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156508) q[0];
sx q[0];
rz(-1.6424718) q[0];
sx q[0];
rz(-2.1732251) q[0];
rz(2.992606) q[2];
sx q[2];
rz(-0.89846957) q[2];
sx q[2];
rz(-0.42299262) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40618784) q[1];
sx q[1];
rz(-1.5899854) q[1];
sx q[1];
rz(2.4270942) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1547847) q[3];
sx q[3];
rz(-0.41066658) q[3];
sx q[3];
rz(-1.9656903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4867268) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(-0.081710903) q[2];
rz(-0.47406667) q[3];
sx q[3];
rz(-1.877955) q[3];
sx q[3];
rz(1.6430395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.896647) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(-0.0078049302) q[0];
rz(1.4004978) q[1];
sx q[1];
rz(-0.85406071) q[1];
sx q[1];
rz(-1.1046462) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7846851) q[0];
sx q[0];
rz(-2.4570358) q[0];
sx q[0];
rz(0.25951578) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8137663) q[2];
sx q[2];
rz(-2.1967078) q[2];
sx q[2];
rz(1.6210131) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0749952) q[1];
sx q[1];
rz(-2.1146333) q[1];
sx q[1];
rz(0.72504136) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8940582) q[3];
sx q[3];
rz(-0.48719104) q[3];
sx q[3];
rz(-2.4647453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2509987) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(2.5863623) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.3135066) q[3];
sx q[3];
rz(1.593332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5531439) q[0];
sx q[0];
rz(-2.6597436) q[0];
sx q[0];
rz(-3.0798262) q[0];
rz(0.24208367) q[1];
sx q[1];
rz(-0.37529072) q[1];
sx q[1];
rz(1.1118836) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0430849) q[0];
sx q[0];
rz(-0.83139172) q[0];
sx q[0];
rz(0.99494536) q[0];
x q[1];
rz(0.013826088) q[2];
sx q[2];
rz(-1.6484043) q[2];
sx q[2];
rz(-1.8344049) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.249914) q[1];
sx q[1];
rz(-2.2644682) q[1];
sx q[1];
rz(-1.9338884) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8133409) q[3];
sx q[3];
rz(-1.8040856) q[3];
sx q[3];
rz(1.0321898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8093439) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(-2.3279482) q[2];
rz(1.404473) q[3];
sx q[3];
rz(-2.9128894) q[3];
sx q[3];
rz(2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-0.62548816) q[0];
sx q[0];
rz(-1.3502716) q[0];
sx q[0];
rz(-1.425449) q[0];
rz(1.6199934) q[1];
sx q[1];
rz(-0.75192538) q[1];
sx q[1];
rz(-2.5040748) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.92256) q[0];
sx q[0];
rz(-1.0244644) q[0];
sx q[0];
rz(-0.68740293) q[0];
rz(-pi) q[1];
rz(-1.2136739) q[2];
sx q[2];
rz(-2.2631858) q[2];
sx q[2];
rz(1.1536319) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.058640826) q[1];
sx q[1];
rz(-2.4070027) q[1];
sx q[1];
rz(-2.525108) q[1];
rz(-pi) q[2];
rz(0.57070891) q[3];
sx q[3];
rz(-1.686704) q[3];
sx q[3];
rz(-1.4012208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1311538) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(2.0054224) q[2];
rz(-1.4853959) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(-1.2095399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5159601) q[0];
sx q[0];
rz(-0.78173286) q[0];
sx q[0];
rz(0.67614722) q[0];
rz(-2.3162084) q[1];
sx q[1];
rz(-0.4239347) q[1];
sx q[1];
rz(-2.1264145) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7207986) q[0];
sx q[0];
rz(-1.4402698) q[0];
sx q[0];
rz(0.41850787) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4887772) q[2];
sx q[2];
rz(-1.4482933) q[2];
sx q[2];
rz(2.0147689) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3489192) q[1];
sx q[1];
rz(-1.9209849) q[1];
sx q[1];
rz(2.7677571) q[1];
rz(2.3750651) q[3];
sx q[3];
rz(-1.818728) q[3];
sx q[3];
rz(-0.31162308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4201346) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(2.1742415) q[2];
rz(1.5445276) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(-2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6185146) q[0];
sx q[0];
rz(-1.5627562) q[0];
sx q[0];
rz(0.05649795) q[0];
rz(2.0691195) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(1.4046232) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4182189) q[0];
sx q[0];
rz(-1.9310111) q[0];
sx q[0];
rz(0.13803137) q[0];
rz(-pi) q[1];
rz(-1.4490453) q[2];
sx q[2];
rz(-1.8728313) q[2];
sx q[2];
rz(-1.9027325) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9060668) q[1];
sx q[1];
rz(-2.0350254) q[1];
sx q[1];
rz(0.81495754) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0839755) q[3];
sx q[3];
rz(-0.15078292) q[3];
sx q[3];
rz(-2.5654716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.848032) q[2];
sx q[2];
rz(-2.3876987) q[2];
sx q[2];
rz(2.8354697) q[2];
rz(-0.39811578) q[3];
sx q[3];
rz(-1.476215) q[3];
sx q[3];
rz(-2.117363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33481471) q[0];
sx q[0];
rz(-2.0756742) q[0];
sx q[0];
rz(-2.6609127) q[0];
rz(-1.9819992) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(2.2904916) q[2];
sx q[2];
rz(-2.7534178) q[2];
sx q[2];
rz(-0.22321246) q[2];
rz(1.0767827) q[3];
sx q[3];
rz(-1.8643338) q[3];
sx q[3];
rz(-2.1716933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];