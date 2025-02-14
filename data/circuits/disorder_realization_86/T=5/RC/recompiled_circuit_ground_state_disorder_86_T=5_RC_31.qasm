OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3616537) q[0];
sx q[0];
rz(-3.0335463) q[0];
sx q[0];
rz(-1.894423) q[0];
rz(-0.78020686) q[1];
sx q[1];
rz(8.297774) q[1];
sx q[1];
rz(10.170593) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0184263) q[0];
sx q[0];
rz(-1.2300876) q[0];
sx q[0];
rz(2.9433204) q[0];
rz(-pi) q[1];
rz(-1.4100513) q[2];
sx q[2];
rz(-1.5622592) q[2];
sx q[2];
rz(-1.05913) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8565897) q[1];
sx q[1];
rz(-1.3599281) q[1];
sx q[1];
rz(-3.0233356) q[1];
rz(-pi) q[2];
rz(0.28663825) q[3];
sx q[3];
rz(-0.53813808) q[3];
sx q[3];
rz(1.9236444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9790393) q[2];
sx q[2];
rz(-0.59430846) q[2];
sx q[2];
rz(-1.7621367) q[2];
rz(0.72921324) q[3];
sx q[3];
rz(-2.3766434) q[3];
sx q[3];
rz(-2.2123857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5586229) q[0];
sx q[0];
rz(-2.0780777) q[0];
sx q[0];
rz(2.9003918) q[0];
rz(-0.25602117) q[1];
sx q[1];
rz(-2.119901) q[1];
sx q[1];
rz(-0.78114885) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22750073) q[0];
sx q[0];
rz(-2.1315398) q[0];
sx q[0];
rz(0.95277159) q[0];
x q[1];
rz(-2.5767869) q[2];
sx q[2];
rz(-1.5283268) q[2];
sx q[2];
rz(-0.04893411) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95805146) q[1];
sx q[1];
rz(-2.1528917) q[1];
sx q[1];
rz(-0.84611012) q[1];
rz(-pi) q[2];
rz(-3.097104) q[3];
sx q[3];
rz(-1.4334363) q[3];
sx q[3];
rz(3.1166164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8476734) q[2];
sx q[2];
rz(-1.5519698) q[2];
sx q[2];
rz(-1.5839362) q[2];
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
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3874338) q[0];
sx q[0];
rz(-1.4492946) q[0];
sx q[0];
rz(-0.19101983) q[0];
rz(-2.7905131) q[1];
sx q[1];
rz(-0.43126884) q[1];
sx q[1];
rz(1.692159) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3397749) q[0];
sx q[0];
rz(-1.9097985) q[0];
sx q[0];
rz(1.6247686) q[0];
rz(-pi) q[1];
rz(-1.1507785) q[2];
sx q[2];
rz(-1.1598829) q[2];
sx q[2];
rz(-2.0793629) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9026824) q[1];
sx q[1];
rz(-2.604831) q[1];
sx q[1];
rz(0.5384268) q[1];
x q[2];
rz(2.6552917) q[3];
sx q[3];
rz(-2.4913906) q[3];
sx q[3];
rz(0.40806684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8423975) q[2];
sx q[2];
rz(-1.7123875) q[2];
sx q[2];
rz(3.0136285) q[2];
rz(1.9757102) q[3];
sx q[3];
rz(-2.2539625) q[3];
sx q[3];
rz(0.33795801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2853476) q[0];
sx q[0];
rz(-2.0576394) q[0];
sx q[0];
rz(-1.6612843) q[0];
rz(-3.0604494) q[1];
sx q[1];
rz(-1.0973884) q[1];
sx q[1];
rz(0.88044423) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9362088) q[0];
sx q[0];
rz(-0.025553731) q[0];
sx q[0];
rz(-0.25746246) q[0];
x q[1];
rz(0.049186695) q[2];
sx q[2];
rz(-0.95251673) q[2];
sx q[2];
rz(-0.82727369) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9040981) q[1];
sx q[1];
rz(-1.4597055) q[1];
sx q[1];
rz(-2.6059125) q[1];
x q[2];
rz(-1.5945928) q[3];
sx q[3];
rz(-1.0131944) q[3];
sx q[3];
rz(-1.8228643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.60260281) q[2];
sx q[2];
rz(-2.4675641) q[2];
sx q[2];
rz(1.7146141) q[2];
rz(1.1566628) q[3];
sx q[3];
rz(-1.4256698) q[3];
sx q[3];
rz(-2.9837933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91931525) q[0];
sx q[0];
rz(-0.8150402) q[0];
sx q[0];
rz(0.94006938) q[0];
rz(-2.6433511) q[1];
sx q[1];
rz(-0.91612852) q[1];
sx q[1];
rz(-2.0169651) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8064708) q[0];
sx q[0];
rz(-2.2735519) q[0];
sx q[0];
rz(-0.051856144) q[0];
rz(-pi) q[1];
rz(-2.9369135) q[2];
sx q[2];
rz(-1.9726557) q[2];
sx q[2];
rz(-3.0167442) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7371205) q[1];
sx q[1];
rz(-2.0029161) q[1];
sx q[1];
rz(-2.4678556) q[1];
x q[2];
rz(1.5920958) q[3];
sx q[3];
rz(-0.7853342) q[3];
sx q[3];
rz(-1.9404279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.70790946) q[2];
sx q[2];
rz(-0.83432546) q[2];
sx q[2];
rz(1.679812) q[2];
rz(-0.24623571) q[3];
sx q[3];
rz(-0.88053954) q[3];
sx q[3];
rz(2.06854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6980625) q[0];
sx q[0];
rz(-1.9140697) q[0];
sx q[0];
rz(-0.45502934) q[0];
rz(0.21806923) q[1];
sx q[1];
rz(-2.2360305) q[1];
sx q[1];
rz(0.47854447) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45947507) q[0];
sx q[0];
rz(-2.2322901) q[0];
sx q[0];
rz(1.6949928) q[0];
rz(0.46314132) q[2];
sx q[2];
rz(-0.69370334) q[2];
sx q[2];
rz(0.59554539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6887321) q[1];
sx q[1];
rz(-2.0264027) q[1];
sx q[1];
rz(2.7745926) q[1];
x q[2];
rz(-0.68799893) q[3];
sx q[3];
rz(-1.08537) q[3];
sx q[3];
rz(-1.7958876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9559418) q[2];
sx q[2];
rz(-1.4943244) q[2];
sx q[2];
rz(2.4556124) q[2];
rz(1.8020804) q[3];
sx q[3];
rz(-1.9728262) q[3];
sx q[3];
rz(-1.3919938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5550391) q[0];
sx q[0];
rz(-1.3864484) q[0];
sx q[0];
rz(-2.6626124) q[0];
rz(2.2827177) q[1];
sx q[1];
rz(-2.5996467) q[1];
sx q[1];
rz(3.0368793) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090957932) q[0];
sx q[0];
rz(-0.43301157) q[0];
sx q[0];
rz(2.6754624) q[0];
x q[1];
rz(-1.9373879) q[2];
sx q[2];
rz(-2.0288003) q[2];
sx q[2];
rz(1.4187494) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5471902) q[1];
sx q[1];
rz(-2.5195751) q[1];
sx q[1];
rz(-1.9474496) q[1];
rz(1.0927622) q[3];
sx q[3];
rz(-1.6217188) q[3];
sx q[3];
rz(1.6855406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0844448) q[2];
sx q[2];
rz(-1.2556602) q[2];
sx q[2];
rz(0.88195938) q[2];
rz(0.070934892) q[3];
sx q[3];
rz(-1.8788012) q[3];
sx q[3];
rz(-0.96496636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.98081812) q[0];
sx q[0];
rz(-2.6698298) q[0];
sx q[0];
rz(-1.2534575) q[0];
rz(0.71022931) q[1];
sx q[1];
rz(-1.1980779) q[1];
sx q[1];
rz(1.089878) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37140815) q[0];
sx q[0];
rz(-0.84781269) q[0];
sx q[0];
rz(-0.56130479) q[0];
rz(1.9148255) q[2];
sx q[2];
rz(-1.9967177) q[2];
sx q[2];
rz(-2.3461208) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0347669) q[1];
sx q[1];
rz(-2.2112101) q[1];
sx q[1];
rz(0.57042112) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4185216) q[3];
sx q[3];
rz(-2.141905) q[3];
sx q[3];
rz(2.603797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7465884) q[2];
sx q[2];
rz(-2.2825664) q[2];
sx q[2];
rz(-2.11917) q[2];
rz(0.59631452) q[3];
sx q[3];
rz(-2.3583581) q[3];
sx q[3];
rz(2.7885126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31350598) q[0];
sx q[0];
rz(-2.6441898) q[0];
sx q[0];
rz(1.585438) q[0];
rz(1.4900788) q[1];
sx q[1];
rz(-2.6863292) q[1];
sx q[1];
rz(2.1447287) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23295924) q[0];
sx q[0];
rz(-1.6413488) q[0];
sx q[0];
rz(-0.34640991) q[0];
x q[1];
rz(-1.385594) q[2];
sx q[2];
rz(-1.5723937) q[2];
sx q[2];
rz(-1.8271556) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.026406757) q[1];
sx q[1];
rz(-0.86167012) q[1];
sx q[1];
rz(-0.89964189) q[1];
x q[2];
rz(0.50052754) q[3];
sx q[3];
rz(-2.1033035) q[3];
sx q[3];
rz(-0.73126572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8263714) q[2];
sx q[2];
rz(-0.73306495) q[2];
sx q[2];
rz(0.23615393) q[2];
rz(-2.835623) q[3];
sx q[3];
rz(-1.1924084) q[3];
sx q[3];
rz(1.6570305) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9295101) q[0];
sx q[0];
rz(-0.65123737) q[0];
sx q[0];
rz(-0.65810743) q[0];
rz(-1.2492389) q[1];
sx q[1];
rz(-0.58964261) q[1];
sx q[1];
rz(-2.6499937) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.126287) q[0];
sx q[0];
rz(-1.1231866) q[0];
sx q[0];
rz(0.60430312) q[0];
rz(0.56546338) q[2];
sx q[2];
rz(-1.3817543) q[2];
sx q[2];
rz(0.42542419) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39855865) q[1];
sx q[1];
rz(-1.5614206) q[1];
sx q[1];
rz(0.27573632) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0594143) q[3];
sx q[3];
rz(-2.058147) q[3];
sx q[3];
rz(0.10457071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3601941) q[2];
sx q[2];
rz(-0.34803826) q[2];
sx q[2];
rz(1.6602824) q[2];
rz(2.4252452) q[3];
sx q[3];
rz(-2.4951388) q[3];
sx q[3];
rz(-0.20814482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3502055) q[0];
sx q[0];
rz(-1.5972142) q[0];
sx q[0];
rz(0.41440339) q[0];
rz(2.1898337) q[1];
sx q[1];
rz(-1.2928243) q[1];
sx q[1];
rz(-1.181319) q[1];
rz(-2.4454115) q[2];
sx q[2];
rz(-2.0700818) q[2];
sx q[2];
rz(-2.5250057) q[2];
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
