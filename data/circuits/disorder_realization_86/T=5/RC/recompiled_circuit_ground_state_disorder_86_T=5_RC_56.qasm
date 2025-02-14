OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.779939) q[0];
sx q[0];
rz(-0.10804636) q[0];
sx q[0];
rz(-1.2471696) q[0];
rz(-0.78020686) q[1];
sx q[1];
rz(-1.127004) q[1];
sx q[1];
rz(2.3957774) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5818258) q[0];
sx q[0];
rz(-2.7493461) q[0];
sx q[0];
rz(-2.0779209) q[0];
rz(3.132944) q[2];
sx q[2];
rz(-1.4100572) q[2];
sx q[2];
rz(0.5102821) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.31065658) q[1];
sx q[1];
rz(-1.6864221) q[1];
sx q[1];
rz(-1.7831037) q[1];
rz(1.7379825) q[3];
sx q[3];
rz(-2.0847528) q[3];
sx q[3];
rz(1.592976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9790393) q[2];
sx q[2];
rz(-2.5472842) q[2];
sx q[2];
rz(-1.3794559) q[2];
rz(0.72921324) q[3];
sx q[3];
rz(-0.76494923) q[3];
sx q[3];
rz(2.2123857) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58296975) q[0];
sx q[0];
rz(-1.0635149) q[0];
sx q[0];
rz(0.24120086) q[0];
rz(-2.8855715) q[1];
sx q[1];
rz(-1.0216917) q[1];
sx q[1];
rz(-0.78114885) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9140919) q[0];
sx q[0];
rz(-1.0100528) q[0];
sx q[0];
rz(0.95277159) q[0];
x q[1];
rz(-3.0623661) q[2];
sx q[2];
rz(-2.5753655) q[2];
sx q[2];
rz(1.5528284) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9731739) q[1];
sx q[1];
rz(-2.2464464) q[1];
sx q[1];
rz(-2.3526117) q[1];
rz(-pi) q[2];
rz(-1.433302) q[3];
sx q[3];
rz(-1.6148657) q[3];
sx q[3];
rz(1.5397244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8476734) q[2];
sx q[2];
rz(-1.5519698) q[2];
sx q[2];
rz(-1.5576564) q[2];
rz(3.1368351) q[3];
sx q[3];
rz(-0.59499756) q[3];
sx q[3];
rz(-2.7958272) q[3];
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
rz(-pi/2) q[3];
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
rz(1.3874338) q[0];
sx q[0];
rz(-1.6922981) q[0];
sx q[0];
rz(0.19101983) q[0];
rz(-2.7905131) q[1];
sx q[1];
rz(-0.43126884) q[1];
sx q[1];
rz(-1.4494337) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3397749) q[0];
sx q[0];
rz(-1.2317941) q[0];
sx q[0];
rz(1.6247686) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6963553) q[2];
sx q[2];
rz(-1.1876593) q[2];
sx q[2];
rz(-0.68510011) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3354794) q[1];
sx q[1];
rz(-1.8361143) q[1];
sx q[1];
rz(-0.47228864) q[1];
rz(-pi) q[2];
rz(1.9123069) q[3];
sx q[3];
rz(-2.1355029) q[3];
sx q[3];
rz(-0.17819861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8423975) q[2];
sx q[2];
rz(-1.7123875) q[2];
sx q[2];
rz(-3.0136285) q[2];
rz(-1.1658824) q[3];
sx q[3];
rz(-0.88763014) q[3];
sx q[3];
rz(-0.33795801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8562451) q[0];
sx q[0];
rz(-1.0839533) q[0];
sx q[0];
rz(1.4803084) q[0];
rz(3.0604494) q[1];
sx q[1];
rz(-2.0442043) q[1];
sx q[1];
rz(-2.2611484) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9362088) q[0];
sx q[0];
rz(-0.025553731) q[0];
sx q[0];
rz(-0.25746246) q[0];
rz(3.092406) q[2];
sx q[2];
rz(-0.95251673) q[2];
sx q[2];
rz(0.82727369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26759155) q[1];
sx q[1];
rz(-1.0387712) q[1];
sx q[1];
rz(1.4417955) q[1];
rz(-pi) q[2];
rz(1.5945928) q[3];
sx q[3];
rz(-2.1283983) q[3];
sx q[3];
rz(1.3187283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.60260281) q[2];
sx q[2];
rz(-2.4675641) q[2];
sx q[2];
rz(1.7146141) q[2];
rz(1.1566628) q[3];
sx q[3];
rz(-1.4256698) q[3];
sx q[3];
rz(0.15779933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91931525) q[0];
sx q[0];
rz(-2.3265525) q[0];
sx q[0];
rz(2.2015233) q[0];
rz(-0.4982416) q[1];
sx q[1];
rz(-2.2254641) q[1];
sx q[1];
rz(1.1246276) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33512188) q[0];
sx q[0];
rz(-0.86804077) q[0];
sx q[0];
rz(3.0897365) q[0];
rz(-pi) q[1];
rz(1.1246936) q[2];
sx q[2];
rz(-0.44844018) q[2];
sx q[2];
rz(2.5288127) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4044721) q[1];
sx q[1];
rz(-1.1386765) q[1];
sx q[1];
rz(0.67373709) q[1];
rz(-pi) q[2];
rz(1.5494969) q[3];
sx q[3];
rz(-2.3562585) q[3];
sx q[3];
rz(1.2011647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.70790946) q[2];
sx q[2];
rz(-2.3072672) q[2];
sx q[2];
rz(1.4617807) q[2];
rz(-0.24623571) q[3];
sx q[3];
rz(-2.2610531) q[3];
sx q[3];
rz(1.0730526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(1.6980625) q[0];
sx q[0];
rz(-1.9140697) q[0];
sx q[0];
rz(-0.45502934) q[0];
rz(2.9235234) q[1];
sx q[1];
rz(-0.90556216) q[1];
sx q[1];
rz(-2.6630482) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1878596) q[0];
sx q[0];
rz(-1.6687013) q[0];
sx q[0];
rz(-0.66524532) q[0];
rz(-1.9265106) q[2];
sx q[2];
rz(-0.96181574) q[2];
sx q[2];
rz(0.019542309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8560972) q[1];
sx q[1];
rz(-1.8988653) q[1];
sx q[1];
rz(-1.0873919) q[1];
rz(-pi) q[2];
rz(-2.4535937) q[3];
sx q[3];
rz(-2.0562226) q[3];
sx q[3];
rz(-1.7958876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1856508) q[2];
sx q[2];
rz(-1.4943244) q[2];
sx q[2];
rz(0.68598023) q[2];
rz(-1.3395122) q[3];
sx q[3];
rz(-1.9728262) q[3];
sx q[3];
rz(1.7495988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5550391) q[0];
sx q[0];
rz(-1.3864484) q[0];
sx q[0];
rz(2.6626124) q[0];
rz(2.2827177) q[1];
sx q[1];
rz(-2.5996467) q[1];
sx q[1];
rz(3.0368793) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0506347) q[0];
sx q[0];
rz(-2.7085811) q[0];
sx q[0];
rz(2.6754624) q[0];
rz(-pi) q[1];
rz(-2.6557561) q[2];
sx q[2];
rz(-1.8981032) q[2];
sx q[2];
rz(-0.32020928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4762209) q[1];
sx q[1];
rz(-1.7867861) q[1];
sx q[1];
rz(2.1588227) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0927622) q[3];
sx q[3];
rz(-1.5198738) q[3];
sx q[3];
rz(1.6855406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0571478) q[2];
sx q[2];
rz(-1.8859325) q[2];
sx q[2];
rz(-2.2596333) q[2];
rz(0.070934892) q[3];
sx q[3];
rz(-1.8788012) q[3];
sx q[3];
rz(-0.96496636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1607745) q[0];
sx q[0];
rz(-2.6698298) q[0];
sx q[0];
rz(1.2534575) q[0];
rz(0.71022931) q[1];
sx q[1];
rz(-1.1980779) q[1];
sx q[1];
rz(-2.0517147) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37140815) q[0];
sx q[0];
rz(-0.84781269) q[0];
sx q[0];
rz(0.56130479) q[0];
rz(-pi) q[1];
rz(2.6925025) q[2];
sx q[2];
rz(-1.2586) q[2];
sx q[2];
rz(0.62836601) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8556887) q[1];
sx q[1];
rz(-0.82997417) q[1];
sx q[1];
rz(0.94372933) q[1];
x q[2];
rz(1.4185216) q[3];
sx q[3];
rz(-0.99968761) q[3];
sx q[3];
rz(-0.53779569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.39500427) q[2];
sx q[2];
rz(-2.2825664) q[2];
sx q[2];
rz(-1.0224226) q[2];
rz(-2.5452781) q[3];
sx q[3];
rz(-0.78323451) q[3];
sx q[3];
rz(-2.7885126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8280867) q[0];
sx q[0];
rz(-2.6441898) q[0];
sx q[0];
rz(1.5561546) q[0];
rz(1.4900788) q[1];
sx q[1];
rz(-0.45526344) q[1];
sx q[1];
rz(-2.1447287) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3123967) q[0];
sx q[0];
rz(-1.2252843) q[0];
sx q[0];
rz(1.4958044) q[0];
rz(-0.0016251621) q[2];
sx q[2];
rz(-1.3855943) q[2];
sx q[2];
rz(-2.8855326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0671419) q[1];
sx q[1];
rz(-1.07927) q[1];
sx q[1];
rz(-2.3105826) q[1];
x q[2];
rz(-2.2542265) q[3];
sx q[3];
rz(-2.4277504) q[3];
sx q[3];
rz(0.091370739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8263714) q[2];
sx q[2];
rz(-0.73306495) q[2];
sx q[2];
rz(2.9054387) q[2];
rz(-2.835623) q[3];
sx q[3];
rz(-1.1924084) q[3];
sx q[3];
rz(-1.4845622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9295101) q[0];
sx q[0];
rz(-2.4903553) q[0];
sx q[0];
rz(-2.4834852) q[0];
rz(-1.2492389) q[1];
sx q[1];
rz(-2.55195) q[1];
sx q[1];
rz(2.6499937) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015305681) q[0];
sx q[0];
rz(-1.1231866) q[0];
sx q[0];
rz(-2.5372895) q[0];
rz(-pi) q[1];
rz(-2.7986235) q[2];
sx q[2];
rz(-2.548647) q[2];
sx q[2];
rz(0.85747257) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.39855865) q[1];
sx q[1];
rz(-1.5614206) q[1];
sx q[1];
rz(-2.8658563) q[1];
rz(-2.059547) q[3];
sx q[3];
rz(-1.4982035) q[3];
sx q[3];
rz(-1.5047764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7813985) q[2];
sx q[2];
rz(-2.7935544) q[2];
sx q[2];
rz(1.4813102) q[2];
rz(2.4252452) q[3];
sx q[3];
rz(-0.64645386) q[3];
sx q[3];
rz(-2.9334478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79138712) q[0];
sx q[0];
rz(-1.5443784) q[0];
sx q[0];
rz(-2.7271893) q[0];
rz(-0.95175891) q[1];
sx q[1];
rz(-1.2928243) q[1];
sx q[1];
rz(-1.181319) q[1];
rz(2.1887171) q[2];
sx q[2];
rz(-0.97277736) q[2];
sx q[2];
rz(2.5680367) q[2];
rz(2.9708859) q[3];
sx q[3];
rz(-0.89912631) q[3];
sx q[3];
rz(2.5691433) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
