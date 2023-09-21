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
rz(1.641474) q[1];
sx q[1];
rz(5.2483622) q[1];
sx q[1];
rz(11.622826) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65981728) q[0];
sx q[0];
rz(-3.0660015) q[0];
sx q[0];
rz(2.6600921) q[0];
rz(3.0303454) q[2];
sx q[2];
rz(-2.4876378) q[2];
sx q[2];
rz(-0.36270579) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.125572) q[1];
sx q[1];
rz(-1.2309832) q[1];
sx q[1];
rz(-1.1556975) q[1];
x q[2];
rz(-1.7581851) q[3];
sx q[3];
rz(-1.2710147) q[3];
sx q[3];
rz(-0.0084458394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8865108) q[2];
sx q[2];
rz(-1.3811029) q[2];
sx q[2];
rz(-1.250766) q[2];
rz(1.4261774) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(0.9799408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.132906) q[0];
sx q[0];
rz(-1.0722906) q[0];
sx q[0];
rz(-2.5426478) q[0];
rz(-1.3409746) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(-2.1751931) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37106284) q[0];
sx q[0];
rz(-2.3822228) q[0];
sx q[0];
rz(-0.66803996) q[0];
x q[1];
rz(3.0722926) q[2];
sx q[2];
rz(-0.6233223) q[2];
sx q[2];
rz(-0.22340439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.71827723) q[1];
sx q[1];
rz(-1.4916972) q[1];
sx q[1];
rz(0.11421108) q[1];
rz(-2.0004683) q[3];
sx q[3];
rz(-0.86887348) q[3];
sx q[3];
rz(-2.7817291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.085658375) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(0.30109626) q[2];
rz(-1.9484693) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(-1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10467228) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(1.954129) q[0];
rz(1.9056412) q[1];
sx q[1];
rz(-1.0373479) q[1];
sx q[1];
rz(-1.3175861) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353969) q[0];
sx q[0];
rz(-2.4164696) q[0];
sx q[0];
rz(2.8185185) q[0];
rz(2.927711) q[2];
sx q[2];
rz(-1.7638532) q[2];
sx q[2];
rz(-2.4436827) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2059523) q[1];
sx q[1];
rz(-1.1516358) q[1];
sx q[1];
rz(1.4599035) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1173238) q[3];
sx q[3];
rz(-2.085272) q[3];
sx q[3];
rz(2.8770212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0111771) q[2];
sx q[2];
rz(-1.3935564) q[2];
sx q[2];
rz(2.9349566) q[2];
rz(-2.4335499) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(2.2487683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77465039) q[0];
sx q[0];
rz(-2.9270524) q[0];
sx q[0];
rz(-0.88622093) q[0];
rz(-2.1318502) q[1];
sx q[1];
rz(-0.90615288) q[1];
sx q[1];
rz(-1.9151691) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6521586) q[0];
sx q[0];
rz(-2.2962748) q[0];
sx q[0];
rz(3.0530531) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45986508) q[2];
sx q[2];
rz(-2.2042639) q[2];
sx q[2];
rz(-1.2666653) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.671771) q[1];
sx q[1];
rz(-1.3739532) q[1];
sx q[1];
rz(2.927604) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6552116) q[3];
sx q[3];
rz(-1.9035305) q[3];
sx q[3];
rz(-2.2330724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4975138) q[2];
sx q[2];
rz(-1.8277233) q[2];
sx q[2];
rz(-0.99299661) q[2];
rz(-1.3126866) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(0.48373568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(0.58247724) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156508) q[0];
sx q[0];
rz(-1.4991209) q[0];
sx q[0];
rz(-0.96836758) q[0];
rz(-pi) q[1];
rz(0.8930348) q[2];
sx q[2];
rz(-1.4544011) q[2];
sx q[2];
rz(-1.2410156) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9936258) q[1];
sx q[1];
rz(-2.2851351) q[1];
sx q[1];
rz(-1.5453969) q[1];
rz(2.7353103) q[3];
sx q[3];
rz(-1.5092106) q[3];
sx q[3];
rz(-2.6046034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4867268) q[2];
sx q[2];
rz(-1.606769) q[2];
sx q[2];
rz(-0.081710903) q[2];
rz(0.47406667) q[3];
sx q[3];
rz(-1.877955) q[3];
sx q[3];
rz(-1.6430395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24494568) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(3.1337877) q[0];
rz(-1.4004978) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(2.0369464) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4545126) q[0];
sx q[0];
rz(-2.2283163) q[0];
sx q[0];
rz(1.3643273) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32128895) q[2];
sx q[2];
rz(-2.4761204) q[2];
sx q[2];
rz(-1.9208391) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9340583) q[1];
sx q[1];
rz(-2.174252) q[1];
sx q[1];
rz(0.8912837) q[1];
x q[2];
rz(2.9748627) q[3];
sx q[3];
rz(-1.1108526) q[3];
sx q[3];
rz(-1.0392287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2509987) q[2];
sx q[2];
rz(-2.4119174) q[2];
sx q[2];
rz(2.5863623) q[2];
rz(-2.9688719) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(-1.5482607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5884488) q[0];
sx q[0];
rz(-2.6597436) q[0];
sx q[0];
rz(0.061766457) q[0];
rz(2.899509) q[1];
sx q[1];
rz(-0.37529072) q[1];
sx q[1];
rz(-1.1118836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2568946) q[0];
sx q[0];
rz(-1.1567133) q[0];
sx q[0];
rz(2.3143682) q[0];
rz(-pi) q[1];
rz(-3.1277666) q[2];
sx q[2];
rz(-1.6484043) q[2];
sx q[2];
rz(1.3071878) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5824077) q[1];
sx q[1];
rz(-1.8473986) q[1];
sx q[1];
rz(0.72699593) q[1];
rz(1.3248596) q[3];
sx q[3];
rz(-1.8898367) q[3];
sx q[3];
rz(-0.46003534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3322488) q[2];
sx q[2];
rz(-0.76433864) q[2];
sx q[2];
rz(-2.3279482) q[2];
rz(-1.7371197) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(-2.5261734) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5161045) q[0];
sx q[0];
rz(-1.3502716) q[0];
sx q[0];
rz(1.7161436) q[0];
rz(1.5215993) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(-2.5040748) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051415074) q[0];
sx q[0];
rz(-2.1438103) q[0];
sx q[0];
rz(0.90419241) q[0];
rz(-pi) q[1];
rz(0.72458467) q[2];
sx q[2];
rz(-1.2983592) q[2];
sx q[2];
rz(2.9582634) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70367614) q[1];
sx q[1];
rz(-0.99214593) q[1];
sx q[1];
rz(2.0520567) q[1];
rz(1.7082801) q[3];
sx q[3];
rz(-1.004389) q[3];
sx q[3];
rz(-2.8979104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1311538) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(-1.1361702) q[2];
rz(1.6561967) q[3];
sx q[3];
rz(-2.6172726) q[3];
sx q[3];
rz(1.2095399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.6256325) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(2.4654454) q[0];
rz(0.82538429) q[1];
sx q[1];
rz(-0.4239347) q[1];
sx q[1];
rz(1.0151781) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43468201) q[0];
sx q[0];
rz(-2.7043531) q[0];
sx q[0];
rz(-2.8291563) q[0];
x q[1];
rz(0.19998156) q[2];
sx q[2];
rz(-0.66255424) q[2];
sx q[2];
rz(-0.28550622) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.059541313) q[1];
sx q[1];
rz(-2.6350628) q[1];
sx q[1];
rz(-2.3561213) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3750651) q[3];
sx q[3];
rz(-1.3228647) q[3];
sx q[3];
rz(0.31162308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4201346) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(-2.1742415) q[2];
rz(-1.5445276) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(-0.35287228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6185146) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(-0.05649795) q[0];
rz(2.0691195) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(1.4046232) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79849762) q[0];
sx q[0];
rz(-1.4416749) q[0];
sx q[0];
rz(-1.9341747) q[0];
rz(-pi) q[1];
rz(2.769906) q[2];
sx q[2];
rz(-0.32495299) q[2];
sx q[2];
rz(1.5124958) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7787331) q[1];
sx q[1];
rz(-2.2792788) q[1];
sx q[1];
rz(-0.94019903) q[1];
rz(-3.0839755) q[3];
sx q[3];
rz(-2.9908097) q[3];
sx q[3];
rz(-0.57612102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.29356062) q[2];
sx q[2];
rz(-2.3876987) q[2];
sx q[2];
rz(-2.8354697) q[2];
rz(-2.7434769) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067779) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(1.9819992) q[1];
sx q[1];
rz(-1.0212785) q[1];
sx q[1];
rz(-2.8550128) q[1];
rz(-2.878306) q[2];
sx q[2];
rz(-1.2821715) q[2];
sx q[2];
rz(0.53496219) q[2];
rz(-0.33070926) q[3];
sx q[3];
rz(-1.0996795) q[3];
sx q[3];
rz(-0.75547937) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];