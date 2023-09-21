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
rz(-2.1980481) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4817754) q[0];
sx q[0];
rz(-3.0660015) q[0];
sx q[0];
rz(0.48150058) q[0];
rz(-pi) q[1];
x q[1];
rz(1.655683) q[2];
sx q[2];
rz(-0.92157084) q[2];
sx q[2];
rz(0.50253403) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.700625) q[1];
sx q[1];
rz(-1.1807627) q[1];
sx q[1];
rz(0.36867152) q[1];
rz(0.30480095) q[3];
sx q[3];
rz(-1.3918575) q[3];
sx q[3];
rz(-1.5233056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8865108) q[2];
sx q[2];
rz(-1.7604897) q[2];
sx q[2];
rz(-1.8908267) q[2];
rz(1.4261774) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(-2.1616518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.132906) q[0];
sx q[0];
rz(-1.0722906) q[0];
sx q[0];
rz(-0.59894484) q[0];
rz(-1.3409746) q[1];
sx q[1];
rz(-2.1913765) q[1];
sx q[1];
rz(2.1751931) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9430267) q[0];
sx q[0];
rz(-2.1417924) q[0];
sx q[0];
rz(2.1023554) q[0];
rz(-pi) q[1];
rz(-1.6205377) q[2];
sx q[2];
rz(-0.94919862) q[2];
sx q[2];
rz(0.13812401) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2981373) q[1];
sx q[1];
rz(-1.4569439) q[1];
sx q[1];
rz(-1.4911806) q[1];
x q[2];
rz(0.74921272) q[3];
sx q[3];
rz(-1.8945165) q[3];
sx q[3];
rz(0.92326984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.085658375) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(-0.30109626) q[2];
rz(1.9484693) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10467228) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(1.1874636) q[0];
rz(1.9056412) q[1];
sx q[1];
rz(-1.0373479) q[1];
sx q[1];
rz(1.8240066) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353969) q[0];
sx q[0];
rz(-0.72512308) q[0];
sx q[0];
rz(-0.32307415) q[0];
rz(-pi) q[1];
rz(-1.3733528) q[2];
sx q[2];
rz(-1.3609481) q[2];
sx q[2];
rz(-2.2270577) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31955645) q[1];
sx q[1];
rz(-1.469538) q[1];
sx q[1];
rz(2.7201369) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1173238) q[3];
sx q[3];
rz(-1.0563207) q[3];
sx q[3];
rz(-0.26457149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1304156) q[2];
sx q[2];
rz(-1.3935564) q[2];
sx q[2];
rz(-2.9349566) q[2];
rz(-2.4335499) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(2.2487683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
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
rz(-2.2553717) q[0];
rz(2.1318502) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(1.2264235) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521586) q[0];
sx q[0];
rz(-0.84531784) q[0];
sx q[0];
rz(3.0530531) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0272155) q[2];
sx q[2];
rz(-2.377844) q[2];
sx q[2];
rz(-0.56994146) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.46982161) q[1];
sx q[1];
rz(-1.7676395) q[1];
sx q[1];
rz(-2.927604) q[1];
x q[2];
rz(1.9434483) q[3];
sx q[3];
rz(-1.1132006) q[3];
sx q[3];
rz(-2.3082993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4975138) q[2];
sx q[2];
rz(-1.3138694) q[2];
sx q[2];
rz(-0.99299661) q[2];
rz(-1.8289061) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-0.48373568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.118367) q[0];
sx q[0];
rz(-2.826773) q[0];
sx q[0];
rz(-1.9556048) q[0];
rz(1.7182619) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(-2.5591154) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3956446) q[0];
sx q[0];
rz(-2.1714604) q[0];
sx q[0];
rz(3.0546741) q[0];
x q[1];
rz(-1.7551454) q[2];
sx q[2];
rz(-2.4554688) q[2];
sx q[2];
rz(-0.18649907) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1867265) q[1];
sx q[1];
rz(-2.426882) q[1];
sx q[1];
rz(-0.029280854) q[1];
x q[2];
rz(1.6378239) q[3];
sx q[3];
rz(-1.9762632) q[3];
sx q[3];
rz(-2.1342579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4867268) q[2];
sx q[2];
rz(-1.606769) q[2];
sx q[2];
rz(0.081710903) q[2];
rz(-0.47406667) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(-1.6430395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.896647) q[0];
sx q[0];
rz(-2.0136254) q[0];
sx q[0];
rz(-0.0078049302) q[0];
rz(-1.4004978) q[1];
sx q[1];
rz(-2.2875319) q[1];
sx q[1];
rz(2.0369464) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0110328) q[0];
sx q[0];
rz(-1.4078119) q[0];
sx q[0];
rz(2.4736604) q[0];
rz(1.8137663) q[2];
sx q[2];
rz(-0.94488482) q[2];
sx q[2];
rz(1.6210131) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9756411) q[1];
sx q[1];
rz(-0.8756606) q[1];
sx q[1];
rz(0.73928164) q[1];
rz(2.036318) q[3];
sx q[3];
rz(-1.7200617) q[3];
sx q[3];
rz(-0.6061337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8905939) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(-2.5863623) q[2];
rz(2.9688719) q[3];
sx q[3];
rz(-1.3135066) q[3];
sx q[3];
rz(-1.5482607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5531439) q[0];
sx q[0];
rz(-2.6597436) q[0];
sx q[0];
rz(0.061766457) q[0];
rz(-2.899509) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(-1.1118836) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8098973) q[0];
sx q[0];
rz(-2.2391717) q[0];
sx q[0];
rz(-2.6033127) q[0];
x q[1];
rz(1.7467473) q[2];
sx q[2];
rz(-3.0627652) q[2];
sx q[2];
rz(-1.4836756) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4278533) q[1];
sx q[1];
rz(-0.76875988) q[1];
sx q[1];
rz(-0.40366918) q[1];
x q[2];
rz(-0.63515969) q[3];
sx q[3];
rz(-2.7413462) q[3];
sx q[3];
rz(1.1349585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3322488) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(-2.3279482) q[2];
rz(-1.404473) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5161045) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(-1.7161436) q[0];
rz(1.5215993) q[1];
sx q[1];
rz(-2.3896673) q[1];
sx q[1];
rz(-2.5040748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2190327) q[0];
sx q[0];
rz(-2.1171283) q[0];
sx q[0];
rz(2.4541897) q[0];
rz(-pi) q[1];
rz(-1.9279187) q[2];
sx q[2];
rz(-0.87840688) q[2];
sx q[2];
rz(1.1536319) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9962822) q[1];
sx q[1];
rz(-1.9687555) q[1];
sx q[1];
rz(0.63509649) q[1];
x q[2];
rz(-0.57070891) q[3];
sx q[3];
rz(-1.4548886) q[3];
sx q[3];
rz(1.7403719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0104388) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(-2.0054224) q[2];
rz(-1.4853959) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(-1.2095399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5159601) q[0];
sx q[0];
rz(-2.3598598) q[0];
sx q[0];
rz(-2.4654454) q[0];
rz(2.3162084) q[1];
sx q[1];
rz(-0.4239347) q[1];
sx q[1];
rz(-1.0151781) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092175352) q[0];
sx q[0];
rz(-1.9855238) q[0];
sx q[0];
rz(-1.7134922) q[0];
x q[1];
rz(-2.4887772) q[2];
sx q[2];
rz(-1.4482933) q[2];
sx q[2];
rz(-2.0147689) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0820513) q[1];
sx q[1];
rz(-2.6350628) q[1];
sx q[1];
rz(2.3561213) q[1];
x q[2];
rz(0.76652758) q[3];
sx q[3];
rz(-1.3228647) q[3];
sx q[3];
rz(2.8299696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4201346) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(-0.96735111) q[2];
rz(1.5445276) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(-2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5230781) q[0];
sx q[0];
rz(-1.5627562) q[0];
sx q[0];
rz(-3.0850947) q[0];
rz(-1.0724732) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(-1.7369695) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4182189) q[0];
sx q[0];
rz(-1.2105816) q[0];
sx q[0];
rz(-0.13803137) q[0];
rz(-pi) q[1];
x q[1];
rz(2.769906) q[2];
sx q[2];
rz(-0.32495299) q[2];
sx q[2];
rz(-1.6290968) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.9359293) q[1];
sx q[1];
rz(-2.231039) q[1];
sx q[1];
rz(-2.5388989) q[1];
x q[2];
rz(-0.057617188) q[3];
sx q[3];
rz(-2.9908097) q[3];
sx q[3];
rz(-2.5654716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.848032) q[2];
sx q[2];
rz(-2.3876987) q[2];
sx q[2];
rz(0.30612293) q[2];
rz(2.7434769) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(2.117363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8067779) q[0];
sx q[0];
rz(-2.0756742) q[0];
sx q[0];
rz(-2.6609127) q[0];
rz(-1.9819992) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(0.85110101) q[2];
sx q[2];
rz(-0.38817482) q[2];
sx q[2];
rz(2.9183802) q[2];
rz(-2.8108834) q[3];
sx q[3];
rz(-2.0419131) q[3];
sx q[3];
rz(2.3861133) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];