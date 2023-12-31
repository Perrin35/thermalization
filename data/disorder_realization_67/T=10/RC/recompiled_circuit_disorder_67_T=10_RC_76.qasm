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
rz(0.73097316) q[0];
rz(-1.5001186) q[1];
sx q[1];
rz(-2.1067696) q[1];
sx q[1];
rz(0.94354454) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9991) q[0];
sx q[0];
rz(-1.5038135) q[0];
sx q[0];
rz(-1.5357369) q[0];
rz(3.0303454) q[2];
sx q[2];
rz(-2.4876378) q[2];
sx q[2];
rz(-0.36270579) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.01602068) q[1];
sx q[1];
rz(-1.2309832) q[1];
sx q[1];
rz(-1.9858951) q[1];
x q[2];
rz(-2.8367917) q[3];
sx q[3];
rz(-1.3918575) q[3];
sx q[3];
rz(-1.5233056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8865108) q[2];
sx q[2];
rz(-1.3811029) q[2];
sx q[2];
rz(1.8908267) q[2];
rz(1.4261774) q[3];
sx q[3];
rz(-2.2255247) q[3];
sx q[3];
rz(-0.9799408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
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
rz(0.59894484) q[0];
rz(-1.3409746) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(0.96639955) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37106284) q[0];
sx q[0];
rz(-0.75936985) q[0];
sx q[0];
rz(-2.4735527) q[0];
x q[1];
rz(-0.069300058) q[2];
sx q[2];
rz(-2.5182704) q[2];
sx q[2];
rz(0.22340439) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2981373) q[1];
sx q[1];
rz(-1.6846488) q[1];
sx q[1];
rz(-1.650412) q[1];
rz(-pi) q[2];
rz(-0.74921272) q[3];
sx q[3];
rz(-1.8945165) q[3];
sx q[3];
rz(-0.92326984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0559343) q[2];
sx q[2];
rz(-2.3264383) q[2];
sx q[2];
rz(-2.8404964) q[2];
rz(-1.9484693) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(1.1863856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10467228) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(-1.1874636) q[0];
rz(-1.2359515) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(1.3175861) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31908195) q[0];
sx q[0];
rz(-1.7829478) q[0];
sx q[0];
rz(0.69885079) q[0];
rz(2.927711) q[2];
sx q[2];
rz(-1.7638532) q[2];
sx q[2];
rz(-2.4436827) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9356404) q[1];
sx q[1];
rz(-1.9899568) q[1];
sx q[1];
rz(-1.6816891) q[1];
rz(1.3665974) q[3];
sx q[3];
rz(-0.52650982) q[3];
sx q[3];
rz(-3.112117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0111771) q[2];
sx q[2];
rz(-1.3935564) q[2];
sx q[2];
rz(-0.2066361) q[2];
rz(-0.7080428) q[3];
sx q[3];
rz(-0.20773023) q[3];
sx q[3];
rz(2.2487683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3669423) q[0];
sx q[0];
rz(-0.21454021) q[0];
sx q[0];
rz(-0.88622093) q[0];
rz(-2.1318502) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(1.9151691) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.489434) q[0];
sx q[0];
rz(-0.84531784) q[0];
sx q[0];
rz(3.0530531) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2573651) q[2];
sx q[2];
rz(-1.9366169) q[2];
sx q[2];
rz(-2.5522752) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1434506) q[1];
sx q[1];
rz(-1.7805903) q[1];
sx q[1];
rz(1.7721121) q[1];
rz(1.9434483) q[3];
sx q[3];
rz(-2.028392) q[3];
sx q[3];
rz(2.3082993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6440789) q[2];
sx q[2];
rz(-1.3138694) q[2];
sx q[2];
rz(2.148596) q[2];
rz(1.8289061) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-2.657857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(1.9556048) q[0];
rz(1.7182619) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(-2.5591154) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3956446) q[0];
sx q[0];
rz(-0.97013226) q[0];
sx q[0];
rz(-0.086918513) q[0];
x q[1];
rz(2.2485579) q[2];
sx q[2];
rz(-1.4544011) q[2];
sx q[2];
rz(-1.900577) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40618784) q[1];
sx q[1];
rz(-1.5899854) q[1];
sx q[1];
rz(-0.71449844) q[1];
x q[2];
rz(0.4062823) q[3];
sx q[3];
rz(-1.5092106) q[3];
sx q[3];
rz(-0.53698925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4867268) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(-0.081710903) q[2];
rz(-0.47406667) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(1.4985532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.896647) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(0.0078049302) q[0];
rz(-1.7410949) q[1];
sx q[1];
rz(-0.85406071) q[1];
sx q[1];
rz(-1.1046462) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4545126) q[0];
sx q[0];
rz(-0.91327635) q[0];
sx q[0];
rz(-1.7772654) q[0];
rz(-pi) q[1];
rz(-2.8203037) q[2];
sx q[2];
rz(-2.4761204) q[2];
sx q[2];
rz(1.2207536) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0749952) q[1];
sx q[1];
rz(-2.1146333) q[1];
sx q[1];
rz(2.4165513) q[1];
x q[2];
rz(1.8940582) q[3];
sx q[3];
rz(-2.6544016) q[3];
sx q[3];
rz(0.67684735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8905939) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(2.5863623) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(-1.593332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
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
rz(-1.5531439) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(-3.0798262) q[0];
rz(-0.24208367) q[1];
sx q[1];
rz(-0.37529072) q[1];
sx q[1];
rz(-1.1118836) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0430849) q[0];
sx q[0];
rz(-0.83139172) q[0];
sx q[0];
rz(-0.99494536) q[0];
rz(-pi) q[1];
rz(-1.493181) q[2];
sx q[2];
rz(-1.5570119) q[2];
sx q[2];
rz(2.8790561) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55918499) q[1];
sx q[1];
rz(-1.294194) q[1];
sx q[1];
rz(-0.72699593) q[1];
x q[2];
rz(2.8133409) q[3];
sx q[3];
rz(-1.337507) q[3];
sx q[3];
rz(-2.1094028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8093439) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(-2.3279482) q[2];
rz(-1.7371197) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(0.61541921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5161045) q[0];
sx q[0];
rz(-1.3502716) q[0];
sx q[0];
rz(-1.425449) q[0];
rz(1.5215993) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(0.63751784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.92256) q[0];
sx q[0];
rz(-1.0244644) q[0];
sx q[0];
rz(-0.68740293) q[0];
rz(-pi) q[1];
rz(-1.9279187) q[2];
sx q[2];
rz(-0.87840688) q[2];
sx q[2];
rz(-1.9879607) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.058640826) q[1];
sx q[1];
rz(-2.4070027) q[1];
sx q[1];
rz(0.61648468) q[1];
rz(-pi) q[2];
rz(0.57070891) q[3];
sx q[3];
rz(-1.4548886) q[3];
sx q[3];
rz(1.4012208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0104388) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(1.1361702) q[2];
rz(1.4853959) q[3];
sx q[3];
rz(-2.6172726) q[3];
sx q[3];
rz(-1.2095399) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6256325) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(0.67614722) q[0];
rz(-2.3162084) q[1];
sx q[1];
rz(-0.4239347) q[1];
sx q[1];
rz(-2.1264145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0494173) q[0];
sx q[0];
rz(-1.1560688) q[0];
sx q[0];
rz(-1.7134922) q[0];
rz(2.9416111) q[2];
sx q[2];
rz(-2.4790384) q[2];
sx q[2];
rz(2.8560864) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79267348) q[1];
sx q[1];
rz(-1.2206077) q[1];
sx q[1];
rz(-2.7677571) q[1];
rz(-2.3750651) q[3];
sx q[3];
rz(-1.3228647) q[3];
sx q[3];
rz(-0.31162308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72145808) q[2];
sx q[2];
rz(-2.9253503) q[2];
sx q[2];
rz(2.1742415) q[2];
rz(1.5445276) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(0.35287228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0427986) q[0];
sx q[0];
rz(-0.3846752) q[0];
sx q[0];
rz(1.2205475) q[0];
rz(-pi) q[1];
rz(1.6925473) q[2];
sx q[2];
rz(-1.8728313) q[2];
sx q[2];
rz(1.2388602) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9060668) q[1];
sx q[1];
rz(-2.0350254) q[1];
sx q[1];
rz(2.3266351) q[1];
rz(-pi) q[2];
x q[2];
rz(0.057617188) q[3];
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
rz(-0.30612293) q[2];
rz(-0.39811578) q[3];
sx q[3];
rz(-1.476215) q[3];
sx q[3];
rz(-2.117363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
