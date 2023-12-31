OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(-1.9987885) q[0];
sx q[0];
rz(1.2115275) q[0];
rz(2.9149574) q[1];
sx q[1];
rz(-1.5645138) q[1];
sx q[1];
rz(-0.29830631) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0878108) q[0];
sx q[0];
rz(-2.2349173) q[0];
sx q[0];
rz(-1.640663) q[0];
x q[1];
rz(1.9036129) q[2];
sx q[2];
rz(-1.6593854) q[2];
sx q[2];
rz(-2.7804136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5391985) q[1];
sx q[1];
rz(-2.1089923) q[1];
sx q[1];
rz(-2.5799275) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6739507) q[3];
sx q[3];
rz(-1.4872331) q[3];
sx q[3];
rz(-2.8148357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4937218) q[2];
sx q[2];
rz(-1.9335258) q[2];
sx q[2];
rz(0.99386627) q[2];
rz(0.99938756) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(0.58888155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64269972) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(-0.018218536) q[0];
rz(0.81623626) q[1];
sx q[1];
rz(-1.0304334) q[1];
sx q[1];
rz(-0.47168628) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90549201) q[0];
sx q[0];
rz(-0.24370757) q[0];
sx q[0];
rz(0.36867504) q[0];
x q[1];
rz(1.291044) q[2];
sx q[2];
rz(-2.0591759) q[2];
sx q[2];
rz(-2.5807057) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.7447409) q[1];
sx q[1];
rz(-0.62780118) q[1];
sx q[1];
rz(-2.6746034) q[1];
rz(1.7945812) q[3];
sx q[3];
rz(-0.61962485) q[3];
sx q[3];
rz(0.97298813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5919684) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(2.1726051) q[2];
rz(0.5747059) q[3];
sx q[3];
rz(-2.59022) q[3];
sx q[3];
rz(-1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4330924) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(2.95978) q[0];
rz(2.0388942) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(-1.4556494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91711125) q[0];
sx q[0];
rz(-0.79229504) q[0];
sx q[0];
rz(0.86865058) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4007912) q[2];
sx q[2];
rz(-1.891279) q[2];
sx q[2];
rz(-1.6522811) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1959343) q[1];
sx q[1];
rz(-0.12609005) q[1];
sx q[1];
rz(3.0581711) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5163243) q[3];
sx q[3];
rz(-1.4843974) q[3];
sx q[3];
rz(-0.82739917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.87749798) q[2];
sx q[2];
rz(-1.4870746) q[2];
sx q[2];
rz(0.96763119) q[2];
rz(-2.4140221) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(-0.23770604) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2847292) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(-0.63823429) q[0];
rz(2.0137285) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(-1.9086054) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15105948) q[0];
sx q[0];
rz(-2.6418243) q[0];
sx q[0];
rz(0.14453669) q[0];
rz(2.5523283) q[2];
sx q[2];
rz(-1.1009842) q[2];
sx q[2];
rz(-1.9643009) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7113263) q[1];
sx q[1];
rz(-1.3923936) q[1];
sx q[1];
rz(3.1049411) q[1];
rz(-pi) q[2];
rz(-1.8423555) q[3];
sx q[3];
rz(-2.4304667) q[3];
sx q[3];
rz(-1.7600876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91810742) q[2];
sx q[2];
rz(-1.8635668) q[2];
sx q[2];
rz(-2.3275862) q[2];
rz(2.0984086) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(-1.957318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84213132) q[0];
sx q[0];
rz(-1.786754) q[0];
sx q[0];
rz(-0.89170757) q[0];
rz(1.8978329) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(-0.2125425) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71056238) q[0];
sx q[0];
rz(-0.02598962) q[0];
sx q[0];
rz(2.2463069) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1754155) q[2];
sx q[2];
rz(-1.6024616) q[2];
sx q[2];
rz(-2.4633173) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5571641) q[1];
sx q[1];
rz(-2.2884528) q[1];
sx q[1];
rz(2.3889956) q[1];
rz(-pi) q[2];
rz(-2.6796474) q[3];
sx q[3];
rz(-0.89832234) q[3];
sx q[3];
rz(2.7355821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.033096878) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(0.5212211) q[2];
rz(1.3850348) q[3];
sx q[3];
rz(-1.228046) q[3];
sx q[3];
rz(-1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8591156) q[0];
sx q[0];
rz(-2.2127667) q[0];
sx q[0];
rz(2.916472) q[0];
rz(-1.3549995) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(-0.37757847) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29992732) q[0];
sx q[0];
rz(-2.8746434) q[0];
sx q[0];
rz(-1.4873234) q[0];
rz(-2.8541366) q[2];
sx q[2];
rz(-1.6171347) q[2];
sx q[2];
rz(3.1099144) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0203591) q[1];
sx q[1];
rz(-1.4660144) q[1];
sx q[1];
rz(2.0038414) q[1];
x q[2];
rz(3.1116629) q[3];
sx q[3];
rz(-2.3464591) q[3];
sx q[3];
rz(-0.078775725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1569415) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(-2.6605576) q[2];
rz(-0.40361079) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(0.31479442) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0890546) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(2.9470434) q[0];
rz(-0.21952195) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(-2.887168) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90760224) q[0];
sx q[0];
rz(-1.6798717) q[0];
sx q[0];
rz(1.8926189) q[0];
rz(-0.15433407) q[2];
sx q[2];
rz(-0.96417226) q[2];
sx q[2];
rz(-0.91147214) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9421778) q[1];
sx q[1];
rz(-0.53655469) q[1];
sx q[1];
rz(-2.5903827) q[1];
rz(2.0394578) q[3];
sx q[3];
rz(-1.161876) q[3];
sx q[3];
rz(-1.8917781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1640132) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(-1.3158201) q[2];
rz(-0.87604648) q[3];
sx q[3];
rz(-0.13893572) q[3];
sx q[3];
rz(-2.1379437) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2106237) q[0];
sx q[0];
rz(-0.3759149) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(0.82398206) q[1];
sx q[1];
rz(-0.95183698) q[1];
sx q[1];
rz(-1.5751858) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67714171) q[0];
sx q[0];
rz(-1.8832708) q[0];
sx q[0];
rz(-2.421445) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0303866) q[2];
sx q[2];
rz(-1.724913) q[2];
sx q[2];
rz(2.266778) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7944784) q[1];
sx q[1];
rz(-1.9273259) q[1];
sx q[1];
rz(1.1785248) q[1];
x q[2];
rz(-2.6353587) q[3];
sx q[3];
rz(-1.1971548) q[3];
sx q[3];
rz(0.95844275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.87551293) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(0.27754647) q[2];
rz(1.6905486) q[3];
sx q[3];
rz(-2.6896559) q[3];
sx q[3];
rz(1.014876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9797416) q[0];
sx q[0];
rz(-0.45270544) q[0];
sx q[0];
rz(-1.6850527) q[0];
rz(0.62943554) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(1.1368407) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0658543) q[0];
sx q[0];
rz(-2.542001) q[0];
sx q[0];
rz(-1.0435186) q[0];
rz(-pi) q[1];
rz(2.2898229) q[2];
sx q[2];
rz(-0.53817828) q[2];
sx q[2];
rz(-0.9966419) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.34502621) q[1];
sx q[1];
rz(-1.476164) q[1];
sx q[1];
rz(-1.252853) q[1];
rz(-pi) q[2];
rz(-2.4056899) q[3];
sx q[3];
rz(-1.8527485) q[3];
sx q[3];
rz(-1.499093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64951605) q[2];
sx q[2];
rz(-1.7121544) q[2];
sx q[2];
rz(1.0409522) q[2];
rz(3.1395636) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(0.31203976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.3025538) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(-0.94605207) q[0];
rz(2.229915) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(2.5295703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34110585) q[0];
sx q[0];
rz(-0.90011156) q[0];
sx q[0];
rz(2.1675046) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7962748) q[2];
sx q[2];
rz(-1.3441663) q[2];
sx q[2];
rz(0.40030865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.40572383) q[1];
sx q[1];
rz(-1.4475665) q[1];
sx q[1];
rz(-2.481639) q[1];
rz(-pi) q[2];
rz(0.53491433) q[3];
sx q[3];
rz(-1.8992918) q[3];
sx q[3];
rz(-1.2319777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0845906) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(2.4882312) q[2];
rz(0.35081321) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(-2.4408834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6012797) q[0];
sx q[0];
rz(-1.5355587) q[0];
sx q[0];
rz(-3.0328947) q[0];
rz(0.75469771) q[1];
sx q[1];
rz(-1.8221868) q[1];
sx q[1];
rz(1.6356161) q[1];
rz(-2.1379708) q[2];
sx q[2];
rz(-0.97273172) q[2];
sx q[2];
rz(1.6956971) q[2];
rz(1.3238974) q[3];
sx q[3];
rz(-0.3331475) q[3];
sx q[3];
rz(-3.0039136) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
