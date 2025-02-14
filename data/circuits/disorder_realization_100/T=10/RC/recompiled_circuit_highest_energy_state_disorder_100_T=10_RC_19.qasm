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
rz(1.1178782) q[0];
sx q[0];
rz(-2.0562545) q[0];
sx q[0];
rz(0.35882741) q[0];
rz(1.2598414) q[1];
sx q[1];
rz(-1.0955732) q[1];
sx q[1];
rz(2.1020558) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0721999) q[0];
sx q[0];
rz(-2.2041) q[0];
sx q[0];
rz(1.1266438) q[0];
rz(-pi) q[1];
rz(0.99604692) q[2];
sx q[2];
rz(-0.050499126) q[2];
sx q[2];
rz(2.8157774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.032385739) q[1];
sx q[1];
rz(-0.47697526) q[1];
sx q[1];
rz(2.2677648) q[1];
x q[2];
rz(-0.67010309) q[3];
sx q[3];
rz(-1.6466318) q[3];
sx q[3];
rz(-2.1431654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.871668) q[2];
sx q[2];
rz(-2.2871127) q[2];
sx q[2];
rz(-0.16647767) q[2];
rz(3.0050333) q[3];
sx q[3];
rz(-0.86524335) q[3];
sx q[3];
rz(-1.4599919) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1305552) q[0];
sx q[0];
rz(-2.8189973) q[0];
sx q[0];
rz(2.7815797) q[0];
rz(-0.75045466) q[1];
sx q[1];
rz(-2.2513335) q[1];
sx q[1];
rz(0.043936122) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30653954) q[0];
sx q[0];
rz(-2.9322698) q[0];
sx q[0];
rz(-2.3002982) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92655757) q[2];
sx q[2];
rz(-1.6100002) q[2];
sx q[2];
rz(-2.7317177) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.94823882) q[1];
sx q[1];
rz(-1.3557503) q[1];
sx q[1];
rz(0.7521473) q[1];
rz(-1.3055849) q[3];
sx q[3];
rz(-1.3797292) q[3];
sx q[3];
rz(1.1950584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3540078) q[2];
sx q[2];
rz(-0.44469357) q[2];
sx q[2];
rz(2.1237874) q[2];
rz(-2.9338037) q[3];
sx q[3];
rz(-0.60616797) q[3];
sx q[3];
rz(-0.35473216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.6924234) q[0];
sx q[0];
rz(-2.848931) q[0];
sx q[0];
rz(2.0617275) q[0];
rz(-1.0049459) q[1];
sx q[1];
rz(-0.69098538) q[1];
sx q[1];
rz(-1.8963922) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8522135) q[0];
sx q[0];
rz(-0.32158467) q[0];
sx q[0];
rz(-2.8711161) q[0];
rz(1.7342308) q[2];
sx q[2];
rz(-1.15112) q[2];
sx q[2];
rz(0.070022665) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1670523) q[1];
sx q[1];
rz(-2.5220832) q[1];
sx q[1];
rz(-2.0761247) q[1];
rz(2.9726282) q[3];
sx q[3];
rz(-2.1848483) q[3];
sx q[3];
rz(1.3425296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.24091844) q[2];
sx q[2];
rz(-1.8697048) q[2];
sx q[2];
rz(-0.34350485) q[2];
rz(-0.014723012) q[3];
sx q[3];
rz(-2.4462409) q[3];
sx q[3];
rz(-1.9007614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0111276) q[0];
sx q[0];
rz(-1.0437597) q[0];
sx q[0];
rz(-0.68487942) q[0];
rz(0.20826805) q[1];
sx q[1];
rz(-1.8332278) q[1];
sx q[1];
rz(-2.3484255) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8858636) q[0];
sx q[0];
rz(-2.1128214) q[0];
sx q[0];
rz(-0.11491187) q[0];
rz(2.9689413) q[2];
sx q[2];
rz(-1.358629) q[2];
sx q[2];
rz(1.4868975) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85911688) q[1];
sx q[1];
rz(-2.1819677) q[1];
sx q[1];
rz(-0.019655991) q[1];
x q[2];
rz(-2.0871619) q[3];
sx q[3];
rz(-0.41958365) q[3];
sx q[3];
rz(-2.6100067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76306474) q[2];
sx q[2];
rz(-1.6918162) q[2];
sx q[2];
rz(-2.6587963) q[2];
rz(1.637508) q[3];
sx q[3];
rz(-2.5129694) q[3];
sx q[3];
rz(-0.32046902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.095161) q[0];
sx q[0];
rz(-0.23155364) q[0];
sx q[0];
rz(-3.0401373) q[0];
rz(-2.9265112) q[1];
sx q[1];
rz(-1.9317893) q[1];
sx q[1];
rz(1.869841) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95071424) q[0];
sx q[0];
rz(-2.029065) q[0];
sx q[0];
rz(-0.31744297) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59753363) q[2];
sx q[2];
rz(-1.1206822) q[2];
sx q[2];
rz(-0.91757739) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3569361) q[1];
sx q[1];
rz(-1.1788834) q[1];
sx q[1];
rz(2.980113) q[1];
rz(-pi) q[2];
rz(0.21274626) q[3];
sx q[3];
rz(-1.233056) q[3];
sx q[3];
rz(-2.4500873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7868906) q[2];
sx q[2];
rz(-2.6933647) q[2];
sx q[2];
rz(0.1498214) q[2];
rz(-0.24770728) q[3];
sx q[3];
rz(-2.7713573) q[3];
sx q[3];
rz(-2.3410334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.596452) q[0];
sx q[0];
rz(-0.62262028) q[0];
sx q[0];
rz(-1.8555634) q[0];
rz(2.3226358) q[1];
sx q[1];
rz(-0.30636925) q[1];
sx q[1];
rz(2.9584296) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4990627) q[0];
sx q[0];
rz(-2.1234491) q[0];
sx q[0];
rz(1.3324609) q[0];
rz(-pi) q[1];
rz(2.8122453) q[2];
sx q[2];
rz(-2.9463049) q[2];
sx q[2];
rz(-1.600071) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.10658856) q[1];
sx q[1];
rz(-0.87751203) q[1];
sx q[1];
rz(2.6412852) q[1];
x q[2];
rz(2.6425903) q[3];
sx q[3];
rz(-1.8516282) q[3];
sx q[3];
rz(-2.4256191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5554819) q[2];
sx q[2];
rz(-1.7038944) q[2];
sx q[2];
rz(2.6044031) q[2];
rz(1.4026583) q[3];
sx q[3];
rz(-1.2638998) q[3];
sx q[3];
rz(2.5130443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9888159) q[0];
sx q[0];
rz(-0.39867908) q[0];
sx q[0];
rz(0.11845778) q[0];
rz(-1.1064103) q[1];
sx q[1];
rz(-1.2245919) q[1];
sx q[1];
rz(2.4647663) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40672228) q[0];
sx q[0];
rz(-3.0483735) q[0];
sx q[0];
rz(-0.63424503) q[0];
rz(-pi) q[1];
rz(-2.9666128) q[2];
sx q[2];
rz(-0.6616486) q[2];
sx q[2];
rz(1.6296737) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21444328) q[1];
sx q[1];
rz(-0.24733686) q[1];
sx q[1];
rz(-0.74199583) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6366129) q[3];
sx q[3];
rz(-0.01465791) q[3];
sx q[3];
rz(2.3848118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38830882) q[2];
sx q[2];
rz(-2.1948185) q[2];
sx q[2];
rz(-3.0906265) q[2];
rz(-2.986159) q[3];
sx q[3];
rz(-1.6499465) q[3];
sx q[3];
rz(-1.0203993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8841356) q[0];
sx q[0];
rz(-1.1549042) q[0];
sx q[0];
rz(0.98404032) q[0];
rz(-2.6718196) q[1];
sx q[1];
rz(-1.4642986) q[1];
sx q[1];
rz(2.9205172) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.675908) q[0];
sx q[0];
rz(-3.0455112) q[0];
sx q[0];
rz(0.19739993) q[0];
rz(-pi) q[1];
x q[1];
rz(2.661411) q[2];
sx q[2];
rz(-0.70887414) q[2];
sx q[2];
rz(-0.95979662) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1103777) q[1];
sx q[1];
rz(-1.3450012) q[1];
sx q[1];
rz(1.2658735) q[1];
rz(-0.99346907) q[3];
sx q[3];
rz(-1.1599891) q[3];
sx q[3];
rz(2.0604756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85349637) q[2];
sx q[2];
rz(-0.6820389) q[2];
sx q[2];
rz(-1.3300396) q[2];
rz(0.62411493) q[3];
sx q[3];
rz(-2.1835486) q[3];
sx q[3];
rz(-0.38930711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57988155) q[0];
sx q[0];
rz(-1.6498673) q[0];
sx q[0];
rz(-1.0505744) q[0];
rz(2.1613278) q[1];
sx q[1];
rz(-1.0824208) q[1];
sx q[1];
rz(1.6260653) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3526392) q[0];
sx q[0];
rz(-1.727956) q[0];
sx q[0];
rz(2.609014) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47591803) q[2];
sx q[2];
rz(-1.7810155) q[2];
sx q[2];
rz(1.1016358) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5298898) q[1];
sx q[1];
rz(-1.262959) q[1];
sx q[1];
rz(-0.3605538) q[1];
rz(-1.100176) q[3];
sx q[3];
rz(-1.017316) q[3];
sx q[3];
rz(0.054892232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3773697) q[2];
sx q[2];
rz(-1.2219967) q[2];
sx q[2];
rz(-2.4019305) q[2];
rz(2.9861279) q[3];
sx q[3];
rz(-2.7630617) q[3];
sx q[3];
rz(0.19708656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4231606) q[0];
sx q[0];
rz(-2.6078556) q[0];
sx q[0];
rz(2.0045795) q[0];
rz(2.5088189) q[1];
sx q[1];
rz(-2.4686765) q[1];
sx q[1];
rz(2.4957962) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5006913) q[0];
sx q[0];
rz(-0.37129042) q[0];
sx q[0];
rz(-2.6616606) q[0];
rz(1.8027202) q[2];
sx q[2];
rz(-2.337079) q[2];
sx q[2];
rz(0.58959145) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91612047) q[1];
sx q[1];
rz(-1.5016306) q[1];
sx q[1];
rz(-0.81838496) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0460351) q[3];
sx q[3];
rz(-1.2933795) q[3];
sx q[3];
rz(-0.21464561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.656903) q[2];
sx q[2];
rz(-0.36269665) q[2];
sx q[2];
rz(-3.1032069) q[2];
rz(0.11135993) q[3];
sx q[3];
rz(-0.42698082) q[3];
sx q[3];
rz(0.67489433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4511694) q[0];
sx q[0];
rz(-1.4372062) q[0];
sx q[0];
rz(1.3896598) q[0];
rz(-1.8545064) q[1];
sx q[1];
rz(-0.99936395) q[1];
sx q[1];
rz(-2.0157464) q[1];
rz(-2.1846175) q[2];
sx q[2];
rz(-1.0697094) q[2];
sx q[2];
rz(-0.30164837) q[2];
rz(0.75931924) q[3];
sx q[3];
rz(-1.6990468) q[3];
sx q[3];
rz(2.852462) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
