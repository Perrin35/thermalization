OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4133889) q[0];
sx q[0];
rz(-1.1336741) q[0];
sx q[0];
rz(1.5925621) q[0];
rz(1.6917317) q[1];
sx q[1];
rz(5.6258968) q[1];
sx q[1];
rz(13.110553) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.005851) q[0];
sx q[0];
rz(-0.07677456) q[0];
sx q[0];
rz(-2.6430921) q[0];
x q[1];
rz(0.10138301) q[2];
sx q[2];
rz(-1.5268491) q[2];
sx q[2];
rz(-1.6909388) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.7218329) q[1];
sx q[1];
rz(-0.83336035) q[1];
sx q[1];
rz(2.1268197) q[1];
x q[2];
rz(0.028624264) q[3];
sx q[3];
rz(-1.5531566) q[3];
sx q[3];
rz(-2.6580435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.071775285) q[2];
sx q[2];
rz(-1.2640307) q[2];
sx q[2];
rz(-1.3624181) q[2];
rz(-0.028256265) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(0.79022592) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64489275) q[0];
sx q[0];
rz(-2.4364478) q[0];
sx q[0];
rz(1.9702966) q[0];
rz(2.9303739) q[1];
sx q[1];
rz(-0.44208458) q[1];
sx q[1];
rz(1.7374977) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068721213) q[0];
sx q[0];
rz(-2.255548) q[0];
sx q[0];
rz(1.0147592) q[0];
x q[1];
rz(-0.30971576) q[2];
sx q[2];
rz(-2.6385912) q[2];
sx q[2];
rz(-0.75370698) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.029570015) q[1];
sx q[1];
rz(-0.28407541) q[1];
sx q[1];
rz(-0.082138852) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7066129) q[3];
sx q[3];
rz(-1.0944301) q[3];
sx q[3];
rz(1.7686896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8319548) q[2];
sx q[2];
rz(-1.4761304) q[2];
sx q[2];
rz(0.94397604) q[2];
rz(-2.5850463) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(-1.336162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8495162) q[0];
sx q[0];
rz(-2.3777666) q[0];
sx q[0];
rz(1.7180432) q[0];
rz(0.81958333) q[1];
sx q[1];
rz(-2.3312566) q[1];
sx q[1];
rz(-0.56366411) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99479988) q[0];
sx q[0];
rz(-2.1437862) q[0];
sx q[0];
rz(1.6550199) q[0];
rz(2.5416449) q[2];
sx q[2];
rz(-1.7581345) q[2];
sx q[2];
rz(1.7768605) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1395531) q[1];
sx q[1];
rz(-0.66322749) q[1];
sx q[1];
rz(1.0242277) q[1];
rz(1.1590957) q[3];
sx q[3];
rz(-0.99038306) q[3];
sx q[3];
rz(-1.8713649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6391969) q[2];
sx q[2];
rz(-1.1221308) q[2];
sx q[2];
rz(1.0602661) q[2];
rz(0.075573102) q[3];
sx q[3];
rz(-2.2836756) q[3];
sx q[3];
rz(-0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.32325) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(-2.9484205) q[0];
rz(-1.5441783) q[1];
sx q[1];
rz(-1.9924106) q[1];
sx q[1];
rz(0.65778041) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6000711) q[0];
sx q[0];
rz(-2.9642448) q[0];
sx q[0];
rz(-0.17139165) q[0];
rz(0.31099702) q[2];
sx q[2];
rz(-1.1592835) q[2];
sx q[2];
rz(2.3222773) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0187877) q[1];
sx q[1];
rz(-1.0446761) q[1];
sx q[1];
rz(-0.29463525) q[1];
x q[2];
rz(-1.2781906) q[3];
sx q[3];
rz(-1.8183823) q[3];
sx q[3];
rz(-2.2729371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0956991) q[2];
sx q[2];
rz(-2.7313488) q[2];
sx q[2];
rz(-2.0945385) q[2];
rz(1.0632769) q[3];
sx q[3];
rz(-1.9789109) q[3];
sx q[3];
rz(1.8803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3445774) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(1.0239333) q[0];
rz(-2.3539885) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(-2.5147298) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1889362) q[0];
sx q[0];
rz(-3.1075826) q[0];
sx q[0];
rz(1.0556428) q[0];
rz(-pi) q[1];
x q[1];
rz(2.824653) q[2];
sx q[2];
rz(-1.2845412) q[2];
sx q[2];
rz(-1.1405108) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6080731) q[1];
sx q[1];
rz(-2.7058209) q[1];
sx q[1];
rz(-0.79969745) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1781138) q[3];
sx q[3];
rz(-1.4564449) q[3];
sx q[3];
rz(-1.3254195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2287067) q[2];
sx q[2];
rz(-1.2330981) q[2];
sx q[2];
rz(-1.0181001) q[2];
rz(-0.61156887) q[3];
sx q[3];
rz(-1.3221778) q[3];
sx q[3];
rz(2.1900246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488895) q[0];
sx q[0];
rz(-1.3494116) q[0];
sx q[0];
rz(1.1992136) q[0];
rz(-1.9723643) q[1];
sx q[1];
rz(-2.0904082) q[1];
sx q[1];
rz(0.32454023) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31774662) q[0];
sx q[0];
rz(-2.4577603) q[0];
sx q[0];
rz(-0.86156396) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3338823) q[2];
sx q[2];
rz(-2.0732023) q[2];
sx q[2];
rz(-3.0903357) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.68795825) q[1];
sx q[1];
rz(-2.1046241) q[1];
sx q[1];
rz(-3.1127991) q[1];
x q[2];
rz(2.5472774) q[3];
sx q[3];
rz(-0.50110498) q[3];
sx q[3];
rz(-1.2831812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4679608) q[2];
sx q[2];
rz(-1.3053852) q[2];
sx q[2];
rz(1.8661873) q[2];
rz(-2.0868789) q[3];
sx q[3];
rz(-0.79189363) q[3];
sx q[3];
rz(-1.595343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6773029) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(-1.7328847) q[0];
rz(-0.45577058) q[1];
sx q[1];
rz(-0.20142889) q[1];
sx q[1];
rz(1.2021525) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8257608) q[0];
sx q[0];
rz(-1.9681265) q[0];
sx q[0];
rz(-2.2560918) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7289594) q[2];
sx q[2];
rz(-2.3473783) q[2];
sx q[2];
rz(-2.7318294) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0899883) q[1];
sx q[1];
rz(-3.0172536) q[1];
sx q[1];
rz(-2.405637) q[1];
x q[2];
rz(2.8771411) q[3];
sx q[3];
rz(-1.3925941) q[3];
sx q[3];
rz(-1.4101392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.000164) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(0.84623519) q[2];
rz(-0.49514654) q[3];
sx q[3];
rz(-2.4954002) q[3];
sx q[3];
rz(-1.600986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-3.0705567) q[0];
sx q[0];
rz(-1.2905916) q[0];
sx q[0];
rz(-2.0284247) q[0];
rz(-2.44599) q[1];
sx q[1];
rz(-0.38989392) q[1];
sx q[1];
rz(1.5323458) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42513645) q[0];
sx q[0];
rz(-1.3279337) q[0];
sx q[0];
rz(0.48432414) q[0];
x q[1];
rz(2.3959827) q[2];
sx q[2];
rz(-2.5157305) q[2];
sx q[2];
rz(-1.2127884) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6032226) q[1];
sx q[1];
rz(-0.45062989) q[1];
sx q[1];
rz(1.5249114) q[1];
rz(-pi) q[2];
rz(-2.9701482) q[3];
sx q[3];
rz(-0.60895863) q[3];
sx q[3];
rz(-0.82933784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1476851) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(2.2873986) q[2];
rz(1.2285852) q[3];
sx q[3];
rz(-1.5649786) q[3];
sx q[3];
rz(1.6555697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5937186) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(0.6643995) q[0];
rz(1.1876855) q[1];
sx q[1];
rz(-2.4408051) q[1];
sx q[1];
rz(1.857035) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8878471) q[0];
sx q[0];
rz(-0.77639025) q[0];
sx q[0];
rz(1.1573769) q[0];
x q[1];
rz(2.8753488) q[2];
sx q[2];
rz(-1.9387445) q[2];
sx q[2];
rz(-3.0014696) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.38332332) q[1];
sx q[1];
rz(-0.82793923) q[1];
sx q[1];
rz(-2.5187056) q[1];
rz(2.7515718) q[3];
sx q[3];
rz(-1.5961002) q[3];
sx q[3];
rz(2.7109773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61218843) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(-2.5908296) q[2];
rz(-2.4380056) q[3];
sx q[3];
rz(-2.1313322) q[3];
sx q[3];
rz(-1.6350869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7228912) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(-0.044145949) q[0];
rz(-1.6126397) q[1];
sx q[1];
rz(-1.9020558) q[1];
sx q[1];
rz(-0.77967656) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2432602) q[0];
sx q[0];
rz(-1.4295477) q[0];
sx q[0];
rz(0.54153533) q[0];
rz(-1.0829955) q[2];
sx q[2];
rz(-2.3239845) q[2];
sx q[2];
rz(-1.4873193) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.979094) q[1];
sx q[1];
rz(-1.3884123) q[1];
sx q[1];
rz(2.6136293) q[1];
rz(-pi) q[2];
rz(2.3603504) q[3];
sx q[3];
rz(-1.9133854) q[3];
sx q[3];
rz(2.4398746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.49446517) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(-1.54281) q[2];
rz(-2.4370082) q[3];
sx q[3];
rz(-1.6274118) q[3];
sx q[3];
rz(-0.012044756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.4298532) q[0];
sx q[0];
rz(-1.5562417) q[0];
sx q[0];
rz(-0.65162311) q[0];
rz(1.343887) q[1];
sx q[1];
rz(-1.6603036) q[1];
sx q[1];
rz(2.4659326) q[1];
rz(-0.30998183) q[2];
sx q[2];
rz(-1.9895171) q[2];
sx q[2];
rz(-2.4537357) q[2];
rz(1.4383437) q[3];
sx q[3];
rz(-0.9369029) q[3];
sx q[3];
rz(2.423219) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
