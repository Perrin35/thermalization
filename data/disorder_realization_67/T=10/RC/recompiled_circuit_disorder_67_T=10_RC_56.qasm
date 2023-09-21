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
rz(-1.0597205) q[0];
sx q[0];
rz(-2.4106195) q[0];
rz(1.641474) q[1];
sx q[1];
rz(-1.0348231) q[1];
sx q[1];
rz(2.1980481) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65981728) q[0];
sx q[0];
rz(-0.075591139) q[0];
sx q[0];
rz(-2.6600921) q[0];
rz(-pi) q[1];
rz(1.655683) q[2];
sx q[2];
rz(-2.2200218) q[2];
sx q[2];
rz(-0.50253403) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4409677) q[1];
sx q[1];
rz(-1.96083) q[1];
sx q[1];
rz(0.36867152) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7581851) q[3];
sx q[3];
rz(-1.2710147) q[3];
sx q[3];
rz(-0.0084458394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8865108) q[2];
sx q[2];
rz(-1.7604897) q[2];
sx q[2];
rz(-1.8908267) q[2];
rz(1.4261774) q[3];
sx q[3];
rz(-2.2255247) q[3];
sx q[3];
rz(2.1616518) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0086867) q[0];
sx q[0];
rz(-1.0722906) q[0];
sx q[0];
rz(-2.5426478) q[0];
rz(-1.8006181) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(2.1751931) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9430267) q[0];
sx q[0];
rz(-0.99980027) q[0];
sx q[0];
rz(2.1023554) q[0];
x q[1];
rz(0.069300058) q[2];
sx q[2];
rz(-0.6233223) q[2];
sx q[2];
rz(0.22340439) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2981373) q[1];
sx q[1];
rz(-1.6846488) q[1];
sx q[1];
rz(1.650412) q[1];
rz(-pi) q[2];
rz(0.74921272) q[3];
sx q[3];
rz(-1.8945165) q[3];
sx q[3];
rz(-2.2183228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0559343) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(-2.8404964) q[2];
rz(1.9484693) q[3];
sx q[3];
rz(-1.5501225) q[3];
sx q[3];
rz(-1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0369204) q[0];
sx q[0];
rz(-1.5043229) q[0];
sx q[0];
rz(-1.954129) q[0];
rz(-1.2359515) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(1.3175861) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353969) q[0];
sx q[0];
rz(-2.4164696) q[0];
sx q[0];
rz(-0.32307415) q[0];
rz(-2.927711) q[2];
sx q[2];
rz(-1.3777395) q[2];
sx q[2];
rz(-2.4436827) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2059523) q[1];
sx q[1];
rz(-1.9899568) q[1];
sx q[1];
rz(-1.6816891) q[1];
rz(3.0242689) q[3];
sx q[3];
rz(-2.085272) q[3];
sx q[3];
rz(-2.8770212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1304156) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(0.2066361) q[2];
rz(2.4335499) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(0.89282435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(0.84336908) q[0];
x q[1];
rz(-1.0272155) q[2];
sx q[2];
rz(-2.377844) q[2];
sx q[2];
rz(0.56994146) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.36818477) q[1];
sx q[1];
rz(-2.8518624) q[1];
sx q[1];
rz(0.75399953) q[1];
rz(1.9434483) q[3];
sx q[3];
rz(-1.1132006) q[3];
sx q[3];
rz(-2.3082993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6440789) q[2];
sx q[2];
rz(-1.8277233) q[2];
sx q[2];
rz(-2.148596) q[2];
rz(-1.8289061) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-0.48373568) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(-1.1859878) q[0];
rz(1.4233308) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(-0.58247724) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3956446) q[0];
sx q[0];
rz(-0.97013226) q[0];
sx q[0];
rz(-3.0546741) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3864473) q[2];
sx q[2];
rz(-0.68612387) q[2];
sx q[2];
rz(2.9550936) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7354048) q[1];
sx q[1];
rz(-1.5899854) q[1];
sx q[1];
rz(-2.4270942) q[1];
rz(-pi) q[2];
rz(-1.5037687) q[3];
sx q[3];
rz(-1.1653295) q[3];
sx q[3];
rz(2.1342579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4867268) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(3.0598818) q[2];
rz(-2.667526) q[3];
sx q[3];
rz(-1.877955) q[3];
sx q[3];
rz(-1.6430395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24494568) q[0];
sx q[0];
rz(-2.0136254) q[0];
sx q[0];
rz(-3.1337877) q[0];
rz(-1.4004978) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(-1.1046462) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68708006) q[0];
sx q[0];
rz(-2.2283163) q[0];
sx q[0];
rz(-1.3643273) q[0];
rz(1.3278264) q[2];
sx q[2];
rz(-2.1967078) q[2];
sx q[2];
rz(1.6210131) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9340583) q[1];
sx q[1];
rz(-0.96734069) q[1];
sx q[1];
rz(-2.250309) q[1];
x q[2];
rz(1.2475345) q[3];
sx q[3];
rz(-2.6544016) q[3];
sx q[3];
rz(2.4647453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2509987) q[2];
sx q[2];
rz(-2.4119174) q[2];
sx q[2];
rz(2.5863623) q[2];
rz(-2.9688719) q[3];
sx q[3];
rz(-1.3135066) q[3];
sx q[3];
rz(1.5482607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5884488) q[0];
sx q[0];
rz(-2.6597436) q[0];
sx q[0];
rz(0.061766457) q[0];
rz(0.24208367) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(-1.1118836) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0985078) q[0];
sx q[0];
rz(-2.3102009) q[0];
sx q[0];
rz(0.99494536) q[0];
x q[1];
rz(-3.1277666) q[2];
sx q[2];
rz(-1.4931884) q[2];
sx q[2];
rz(-1.3071878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.249914) q[1];
sx q[1];
rz(-0.87712446) q[1];
sx q[1];
rz(-1.2077043) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63515969) q[3];
sx q[3];
rz(-0.40024647) q[3];
sx q[3];
rz(2.0066341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8093439) q[2];
sx q[2];
rz(-0.76433864) q[2];
sx q[2];
rz(0.81364441) q[2];
rz(-1.7371197) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(-2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5161045) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(1.7161436) q[0];
rz(1.5215993) q[1];
sx q[1];
rz(-0.75192538) q[1];
sx q[1];
rz(-0.63751784) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.92256) q[0];
sx q[0];
rz(-2.1171283) q[0];
sx q[0];
rz(2.4541897) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2136739) q[2];
sx q[2];
rz(-2.2631858) q[2];
sx q[2];
rz(-1.1536319) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70367614) q[1];
sx q[1];
rz(-2.1494467) q[1];
sx q[1];
rz(2.0520567) q[1];
rz(2.9293204) q[3];
sx q[3];
rz(-2.5605154) q[3];
sx q[3];
rz(0.0086590954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1311538) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(2.0054224) q[2];
rz(-1.6561967) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(1.2095399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6256325) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(-2.4654454) q[0];
rz(-2.3162084) q[1];
sx q[1];
rz(-0.4239347) q[1];
sx q[1];
rz(1.0151781) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4207941) q[0];
sx q[0];
rz(-1.7013229) q[0];
sx q[0];
rz(-0.41850787) q[0];
x q[1];
rz(-0.65281547) q[2];
sx q[2];
rz(-1.4482933) q[2];
sx q[2];
rz(-1.1268238) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2296914) q[1];
sx q[1];
rz(-1.9209314) q[1];
sx q[1];
rz(-1.196911) q[1];
x q[2];
rz(1.9087402) q[3];
sx q[3];
rz(-2.3082808) q[3];
sx q[3];
rz(-2.1144707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4201346) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(-2.1742415) q[2];
rz(1.5445276) q[3];
sx q[3];
rz(-1.8287851) q[3];
sx q[3];
rz(2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6185146) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(0.05649795) q[0];
rz(-2.0691195) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(1.7369695) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0987941) q[0];
sx q[0];
rz(-2.7569175) q[0];
sx q[0];
rz(-1.2205475) q[0];
rz(-pi) q[1];
rz(2.769906) q[2];
sx q[2];
rz(-2.8166397) q[2];
sx q[2];
rz(-1.5124958) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.23552588) q[1];
sx q[1];
rz(-1.1065673) q[1];
sx q[1];
rz(0.81495754) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5620473) q[3];
sx q[3];
rz(-1.7213271) q[3];
sx q[3];
rz(0.51784408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.29356062) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(-2.8354697) q[2];
rz(2.7434769) q[3];
sx q[3];
rz(-1.476215) q[3];
sx q[3];
rz(-2.117363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33481471) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(-1.1595935) q[1];
sx q[1];
rz(-1.0212785) q[1];
sx q[1];
rz(-2.8550128) q[1];
rz(1.8691312) q[2];
sx q[2];
rz(-1.8229501) q[2];
sx q[2];
rz(2.0291871) q[2];
rz(2.8108834) q[3];
sx q[3];
rz(-1.0996795) q[3];
sx q[3];
rz(-0.75547937) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
