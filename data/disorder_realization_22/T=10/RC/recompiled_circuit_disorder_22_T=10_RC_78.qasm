OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7093631) q[0];
sx q[0];
rz(-2.1837283) q[0];
sx q[0];
rz(-0.14444484) q[0];
rz(-2.5748409) q[1];
sx q[1];
rz(-2.6161939) q[1];
sx q[1];
rz(2.1638343) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5241961) q[0];
sx q[0];
rz(-2.8688736) q[0];
sx q[0];
rz(0.32193907) q[0];
x q[1];
rz(-2.0852282) q[2];
sx q[2];
rz(-1.8407514) q[2];
sx q[2];
rz(1.1411238) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.3111167) q[1];
sx q[1];
rz(-1.4314572) q[1];
sx q[1];
rz(1.4275101) q[1];
rz(-2.9631859) q[3];
sx q[3];
rz(-2.3462786) q[3];
sx q[3];
rz(2.2291396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4686761) q[2];
sx q[2];
rz(-1.9925995) q[2];
sx q[2];
rz(0.93227512) q[2];
rz(0.19876924) q[3];
sx q[3];
rz(-2.0310183) q[3];
sx q[3];
rz(-2.1762302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4345877) q[0];
sx q[0];
rz(-0.90536896) q[0];
sx q[0];
rz(0.36112753) q[0];
rz(-1.4350285) q[1];
sx q[1];
rz(-1.7838493) q[1];
sx q[1];
rz(2.3235869) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5404856) q[0];
sx q[0];
rz(-1.2955106) q[0];
sx q[0];
rz(2.7406373) q[0];
rz(-pi) q[1];
rz(1.9423219) q[2];
sx q[2];
rz(-1.6852334) q[2];
sx q[2];
rz(0.83100806) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4600196) q[1];
sx q[1];
rz(-0.96833723) q[1];
sx q[1];
rz(2.6006992) q[1];
rz(-pi) q[2];
rz(1.9345476) q[3];
sx q[3];
rz(-1.517429) q[3];
sx q[3];
rz(1.3747017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8890185) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(-1.7209631) q[2];
rz(-1.8255) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(2.7533598) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7398359) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(2.2706568) q[0];
rz(-2.8254106) q[1];
sx q[1];
rz(-2.8600287) q[1];
sx q[1];
rz(2.8443764) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7446211) q[0];
sx q[0];
rz(-0.40980761) q[0];
sx q[0];
rz(2.9817392) q[0];
rz(-pi) q[1];
rz(-0.88163968) q[2];
sx q[2];
rz(-0.92703968) q[2];
sx q[2];
rz(1.8635441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19993648) q[1];
sx q[1];
rz(-1.3841277) q[1];
sx q[1];
rz(0.048528683) q[1];
rz(-2.3171595) q[3];
sx q[3];
rz(-2.094305) q[3];
sx q[3];
rz(1.5200966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7685984) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(-2.7139943) q[2];
rz(1.1887431) q[3];
sx q[3];
rz(-2.5116428) q[3];
sx q[3];
rz(0.52743131) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7052085) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(0.77392459) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(2.4598222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6503158) q[0];
sx q[0];
rz(-2.804545) q[0];
sx q[0];
rz(-2.2469673) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8800456) q[2];
sx q[2];
rz(-0.62710688) q[2];
sx q[2];
rz(2.4583465) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7766429) q[1];
sx q[1];
rz(-0.72307359) q[1];
sx q[1];
rz(0.88100453) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.79822) q[3];
sx q[3];
rz(-1.9150754) q[3];
sx q[3];
rz(1.6791277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.008808) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(2.183389) q[2];
rz(1.0664252) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(-2.7619894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5565857) q[0];
sx q[0];
rz(-1.7946694) q[0];
sx q[0];
rz(-0.56030309) q[0];
rz(2.141748) q[1];
sx q[1];
rz(-0.20345774) q[1];
sx q[1];
rz(-1.6220185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9902089) q[0];
sx q[0];
rz(-2.4575893) q[0];
sx q[0];
rz(0.74599501) q[0];
rz(-pi) q[1];
rz(-2.5298169) q[2];
sx q[2];
rz(-1.099274) q[2];
sx q[2];
rz(-0.39786354) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2677742) q[1];
sx q[1];
rz(-2.3582637) q[1];
sx q[1];
rz(0.28247139) q[1];
rz(1.7468466) q[3];
sx q[3];
rz(-0.62024833) q[3];
sx q[3];
rz(-0.32905096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4524298) q[2];
sx q[2];
rz(-2.5897557) q[2];
sx q[2];
rz(2.9186644) q[2];
rz(3.1068504) q[3];
sx q[3];
rz(-1.3870753) q[3];
sx q[3];
rz(-0.071578659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8054304) q[0];
sx q[0];
rz(-2.8283089) q[0];
sx q[0];
rz(1.0700595) q[0];
rz(1.3609715) q[1];
sx q[1];
rz(-2.762251) q[1];
sx q[1];
rz(-1.8575352) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7362471) q[0];
sx q[0];
rz(-2.1655472) q[0];
sx q[0];
rz(0.88881641) q[0];
rz(2.1806296) q[2];
sx q[2];
rz(-2.6764538) q[2];
sx q[2];
rz(-2.7679408) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9575427) q[1];
sx q[1];
rz(-1.8584538) q[1];
sx q[1];
rz(-2.6581453) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0809903) q[3];
sx q[3];
rz(-2.5488857) q[3];
sx q[3];
rz(0.61471516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47508919) q[2];
sx q[2];
rz(-1.6978426) q[2];
sx q[2];
rz(-2.5040023) q[2];
rz(2.3049138) q[3];
sx q[3];
rz(-1.1058608) q[3];
sx q[3];
rz(1.7975413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5752983) q[0];
sx q[0];
rz(-0.42995444) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(0.42770806) q[1];
sx q[1];
rz(-1.6141012) q[1];
sx q[1];
rz(0.93820757) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.339401) q[0];
sx q[0];
rz(-1.8759449) q[0];
sx q[0];
rz(-1.1423654) q[0];
rz(-1.0023408) q[2];
sx q[2];
rz(-2.217514) q[2];
sx q[2];
rz(-1.2127753) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.92749121) q[1];
sx q[1];
rz(-2.2117105) q[1];
sx q[1];
rz(-0.44968857) q[1];
rz(-pi) q[2];
rz(-0.35258099) q[3];
sx q[3];
rz(-0.73871021) q[3];
sx q[3];
rz(1.1383566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2618835) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(-1.3860469) q[2];
rz(1.322768) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(-3.1125606) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26738527) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(-1.4338795) q[0];
rz(-1.8677615) q[1];
sx q[1];
rz(-2.0055983) q[1];
sx q[1];
rz(2.3103255) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48152637) q[0];
sx q[0];
rz(-2.7114081) q[0];
sx q[0];
rz(-0.35920401) q[0];
rz(-pi) q[1];
rz(-0.93584658) q[2];
sx q[2];
rz(-1.6747464) q[2];
sx q[2];
rz(2.7452591) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23849328) q[1];
sx q[1];
rz(-1.0817173) q[1];
sx q[1];
rz(2.3563983) q[1];
rz(-pi) q[2];
rz(-0.76922272) q[3];
sx q[3];
rz(-2.8260494) q[3];
sx q[3];
rz(-3.0358918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5780118) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(2.1772299) q[2];
rz(1.1635121) q[3];
sx q[3];
rz(-0.61146277) q[3];
sx q[3];
rz(0.65892974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37373856) q[0];
sx q[0];
rz(-0.73497325) q[0];
sx q[0];
rz(-2.1642165) q[0];
rz(1.3865698) q[1];
sx q[1];
rz(-1.3061378) q[1];
sx q[1];
rz(1.1057373) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7439197) q[0];
sx q[0];
rz(-2.0154675) q[0];
sx q[0];
rz(1.4414653) q[0];
rz(-1.703891) q[2];
sx q[2];
rz(-0.67099748) q[2];
sx q[2];
rz(0.69588307) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.22509174) q[1];
sx q[1];
rz(-1.8019925) q[1];
sx q[1];
rz(0.1901615) q[1];
rz(-pi) q[2];
rz(1.2826142) q[3];
sx q[3];
rz(-1.2328706) q[3];
sx q[3];
rz(-1.5087138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5332807) q[2];
sx q[2];
rz(-0.38485843) q[2];
sx q[2];
rz(-2.9917955) q[2];
rz(-1.3730565) q[3];
sx q[3];
rz(-1.405973) q[3];
sx q[3];
rz(1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49366429) q[0];
sx q[0];
rz(-0.51769185) q[0];
sx q[0];
rz(0.1272442) q[0];
rz(1.4808902) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(-0.1677992) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0906592) q[0];
sx q[0];
rz(-1.5899842) q[0];
sx q[0];
rz(-0.030254342) q[0];
rz(-1.7255515) q[2];
sx q[2];
rz(-0.96826474) q[2];
sx q[2];
rz(-0.11520152) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92229453) q[1];
sx q[1];
rz(-2.1592327) q[1];
sx q[1];
rz(1.5625619) q[1];
rz(-pi) q[2];
rz(-0.76370244) q[3];
sx q[3];
rz(-2.9716431) q[3];
sx q[3];
rz(-0.22596879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1148791) q[2];
sx q[2];
rz(-2.202704) q[2];
sx q[2];
rz(-2.3804469) q[2];
rz(0.090027697) q[3];
sx q[3];
rz(-2.138425) q[3];
sx q[3];
rz(-2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54820838) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(0.38800115) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(2.3315196) q[2];
sx q[2];
rz(-2.0740866) q[2];
sx q[2];
rz(-1.496051) q[2];
rz(-1.2896982) q[3];
sx q[3];
rz(-2.5847808) q[3];
sx q[3];
rz(1.9980711) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];