OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8347972) q[0];
sx q[0];
rz(-2.2936294) q[0];
sx q[0];
rz(0.037394878) q[0];
rz(2.1283863) q[1];
sx q[1];
rz(-2.6599045) q[1];
sx q[1];
rz(-1.4451292) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28222154) q[0];
sx q[0];
rz(-1.1423151) q[0];
sx q[0];
rz(1.3708049) q[0];
rz(2.1712473) q[2];
sx q[2];
rz(-3.0171347) q[2];
sx q[2];
rz(-0.70792922) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6283704) q[1];
sx q[1];
rz(-2.0525825) q[1];
sx q[1];
rz(-0.95487736) q[1];
rz(-pi) q[2];
rz(-1.3919034) q[3];
sx q[3];
rz(-0.57235347) q[3];
sx q[3];
rz(1.5327041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7279613) q[2];
sx q[2];
rz(-3.0484338) q[2];
sx q[2];
rz(-1.6348582) q[2];
rz(2.9480751) q[3];
sx q[3];
rz(-2.3573124) q[3];
sx q[3];
rz(-0.94582742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0993318) q[0];
sx q[0];
rz(-1.7758545) q[0];
sx q[0];
rz(-2.9625764) q[0];
rz(1.5366813) q[1];
sx q[1];
rz(-2.8004526) q[1];
sx q[1];
rz(1.5515597) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59373271) q[0];
sx q[0];
rz(-1.6681801) q[0];
sx q[0];
rz(-1.2800526) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29921542) q[2];
sx q[2];
rz(-2.4096074) q[2];
sx q[2];
rz(0.74886403) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4228004) q[1];
sx q[1];
rz(-1.3072291) q[1];
sx q[1];
rz(1.9843319) q[1];
rz(-2.670774) q[3];
sx q[3];
rz(-2.0662466) q[3];
sx q[3];
rz(1.9268056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4701074) q[2];
sx q[2];
rz(-1.476373) q[2];
sx q[2];
rz(0.022484953) q[2];
rz(0.89618987) q[3];
sx q[3];
rz(-0.78138566) q[3];
sx q[3];
rz(-1.5145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80980587) q[0];
sx q[0];
rz(-0.2066732) q[0];
sx q[0];
rz(-1.5326387) q[0];
rz(1.1398075) q[1];
sx q[1];
rz(-2.8228357) q[1];
sx q[1];
rz(0.87361139) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43481091) q[0];
sx q[0];
rz(-1.0651089) q[0];
sx q[0];
rz(-0.75904738) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0192835) q[2];
sx q[2];
rz(-2.1423369) q[2];
sx q[2];
rz(0.71046605) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1908299) q[1];
sx q[1];
rz(-1.0894945) q[1];
sx q[1];
rz(-1.3973305) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6461883) q[3];
sx q[3];
rz(-1.8429379) q[3];
sx q[3];
rz(2.4484602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5736299) q[2];
sx q[2];
rz(-2.6831388) q[2];
sx q[2];
rz(0.47631329) q[2];
rz(1.025398) q[3];
sx q[3];
rz(-1.7849331) q[3];
sx q[3];
rz(-0.29388139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.77182257) q[0];
sx q[0];
rz(-0.85404587) q[0];
sx q[0];
rz(2.5452132) q[0];
rz(2.3884933) q[1];
sx q[1];
rz(-1.3558847) q[1];
sx q[1];
rz(-2.7900043) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10710005) q[0];
sx q[0];
rz(-1.7020683) q[0];
sx q[0];
rz(-2.6126325) q[0];
rz(-pi) q[1];
rz(-2.142316) q[2];
sx q[2];
rz(-2.3722931) q[2];
sx q[2];
rz(0.77924773) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.233577) q[1];
sx q[1];
rz(-1.963841) q[1];
sx q[1];
rz(0.24891757) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80644108) q[3];
sx q[3];
rz(-2.3262352) q[3];
sx q[3];
rz(-1.7585332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.72254649) q[2];
sx q[2];
rz(-1.498035) q[2];
sx q[2];
rz(0.92918116) q[2];
rz(-1.2292713) q[3];
sx q[3];
rz(-3.1049187) q[3];
sx q[3];
rz(-1.8721972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.1790328) q[0];
sx q[0];
rz(-2.5293009) q[0];
sx q[0];
rz(2.6137733) q[0];
rz(2.5714696) q[1];
sx q[1];
rz(-0.56990439) q[1];
sx q[1];
rz(2.747587) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45371374) q[0];
sx q[0];
rz(-1.7648932) q[0];
sx q[0];
rz(2.6537623) q[0];
rz(2.5486331) q[2];
sx q[2];
rz(-1.2622132) q[2];
sx q[2];
rz(-1.9401996) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7841508) q[1];
sx q[1];
rz(-1.3489745) q[1];
sx q[1];
rz(-3.0206969) q[1];
x q[2];
rz(-2.4179055) q[3];
sx q[3];
rz(-0.95382788) q[3];
sx q[3];
rz(-2.8599515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3923308) q[2];
sx q[2];
rz(-1.8529842) q[2];
sx q[2];
rz(-1.1345081) q[2];
rz(-3.115861) q[3];
sx q[3];
rz(-1.0545571) q[3];
sx q[3];
rz(2.499685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57539097) q[0];
sx q[0];
rz(-1.3302777) q[0];
sx q[0];
rz(2.6824685) q[0];
rz(0.8210012) q[1];
sx q[1];
rz(-0.88290015) q[1];
sx q[1];
rz(0.51148907) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59508649) q[0];
sx q[0];
rz(-2.0430142) q[0];
sx q[0];
rz(-2.8123463) q[0];
x q[1];
rz(-1.8051926) q[2];
sx q[2];
rz(-2.2210741) q[2];
sx q[2];
rz(-0.63287193) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3977511) q[1];
sx q[1];
rz(-1.4201179) q[1];
sx q[1];
rz(0.57445261) q[1];
rz(-pi) q[2];
rz(1.4019764) q[3];
sx q[3];
rz(-0.69069117) q[3];
sx q[3];
rz(2.9258779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56200999) q[2];
sx q[2];
rz(-2.0333813) q[2];
sx q[2];
rz(-0.75135922) q[2];
rz(0.56685081) q[3];
sx q[3];
rz(-1.5375429) q[3];
sx q[3];
rz(1.3273299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.421627) q[0];
sx q[0];
rz(-0.1816853) q[0];
sx q[0];
rz(0.0075465329) q[0];
rz(2.201572) q[1];
sx q[1];
rz(-2.4766141) q[1];
sx q[1];
rz(-3.0090581) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7791393) q[0];
sx q[0];
rz(-2.4041688) q[0];
sx q[0];
rz(0.93517424) q[0];
rz(-pi) q[1];
rz(-1.7039677) q[2];
sx q[2];
rz(-1.2107163) q[2];
sx q[2];
rz(-3.0293426) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1906084) q[1];
sx q[1];
rz(-1.3845456) q[1];
sx q[1];
rz(-0.3979072) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9136732) q[3];
sx q[3];
rz(-1.4063121) q[3];
sx q[3];
rz(3.0066688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.995945) q[2];
sx q[2];
rz(-1.6463248) q[2];
sx q[2];
rz(3.0403467) q[2];
rz(-2.4920987) q[3];
sx q[3];
rz(-0.5972623) q[3];
sx q[3];
rz(-0.35972843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41530171) q[0];
sx q[0];
rz(-1.8446209) q[0];
sx q[0];
rz(-1.8743961) q[0];
rz(-1.2665117) q[1];
sx q[1];
rz(-2.9837065) q[1];
sx q[1];
rz(2.3983541) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.04674656) q[0];
sx q[0];
rz(-2.136793) q[0];
sx q[0];
rz(-1.3309684) q[0];
x q[1];
rz(0.79106462) q[2];
sx q[2];
rz(-1.3692489) q[2];
sx q[2];
rz(2.9746051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0580045) q[1];
sx q[1];
rz(-2.8561855) q[1];
sx q[1];
rz(-2.1607375) q[1];
rz(-pi) q[2];
rz(-2.2970639) q[3];
sx q[3];
rz(-2.2457079) q[3];
sx q[3];
rz(0.19899398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1242975) q[2];
sx q[2];
rz(-2.9399019) q[2];
sx q[2];
rz(-1.7151493) q[2];
rz(1.0157478) q[3];
sx q[3];
rz(-1.4371212) q[3];
sx q[3];
rz(-0.82971853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5468686) q[0];
sx q[0];
rz(-0.65383738) q[0];
sx q[0];
rz(3.1256909) q[0];
rz(1.6390027) q[1];
sx q[1];
rz(-2.465261) q[1];
sx q[1];
rz(-0.99448386) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8581945) q[0];
sx q[0];
rz(-2.1680729) q[0];
sx q[0];
rz(-2.6676763) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5076958) q[2];
sx q[2];
rz(-1.1931538) q[2];
sx q[2];
rz(-1.7971731) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82518903) q[1];
sx q[1];
rz(-2.0540407) q[1];
sx q[1];
rz(-0.8512599) q[1];
rz(2.9904305) q[3];
sx q[3];
rz(-2.0796607) q[3];
sx q[3];
rz(-0.14948949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1556306) q[2];
sx q[2];
rz(-1.7113643) q[2];
sx q[2];
rz(-2.963781) q[2];
rz(-1.4583679) q[3];
sx q[3];
rz(-2.7198313) q[3];
sx q[3];
rz(1.2678857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6630702) q[0];
sx q[0];
rz(-1.8926184) q[0];
sx q[0];
rz(-1.857969) q[0];
rz(2.3578857) q[1];
sx q[1];
rz(-0.77373928) q[1];
sx q[1];
rz(-1.3073889) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4625774) q[0];
sx q[0];
rz(-2.0082012) q[0];
sx q[0];
rz(-1.4694197) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4965335) q[2];
sx q[2];
rz(-0.38887244) q[2];
sx q[2];
rz(0.14479862) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2210113) q[1];
sx q[1];
rz(-2.9579109) q[1];
sx q[1];
rz(1.8964975) q[1];
rz(1.0211759) q[3];
sx q[3];
rz(-1.2116287) q[3];
sx q[3];
rz(0.041426126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6002097) q[2];
sx q[2];
rz(-1.4098488) q[2];
sx q[2];
rz(0.48995885) q[2];
rz(-1.1146891) q[3];
sx q[3];
rz(-1.9652818) q[3];
sx q[3];
rz(-2.8118242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020323044) q[0];
sx q[0];
rz(-1.9135495) q[0];
sx q[0];
rz(1.3187153) q[0];
rz(1.2976788) q[1];
sx q[1];
rz(-2.4083125) q[1];
sx q[1];
rz(-0.42793035) q[1];
rz(0.25925706) q[2];
sx q[2];
rz(-1.8198063) q[2];
sx q[2];
rz(1.168269) q[2];
rz(2.9419327) q[3];
sx q[3];
rz(-0.87423751) q[3];
sx q[3];
rz(-2.8578491) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
