OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0754492) q[0];
sx q[0];
rz(-1.0439405) q[0];
sx q[0];
rz(-3.1312842) q[0];
rz(-0.97996867) q[1];
sx q[1];
rz(-1.6966532) q[1];
sx q[1];
rz(0.57656062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0703556) q[0];
sx q[0];
rz(-0.33215392) q[0];
sx q[0];
rz(0.37918703) q[0];
rz(-pi) q[1];
rz(-2.2498807) q[2];
sx q[2];
rz(-2.0713965) q[2];
sx q[2];
rz(-0.5989738) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.7652055) q[1];
sx q[1];
rz(-0.55843267) q[1];
sx q[1];
rz(0.33513481) q[1];
rz(-2.303894) q[3];
sx q[3];
rz(-1.695096) q[3];
sx q[3];
rz(1.4364786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.52446857) q[2];
sx q[2];
rz(-1.6813797) q[2];
sx q[2];
rz(-1.8560393) q[2];
rz(-1.589795) q[3];
sx q[3];
rz(-1.0714622) q[3];
sx q[3];
rz(-2.0901399) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0153506) q[0];
sx q[0];
rz(-1.9168357) q[0];
sx q[0];
rz(-0.24573627) q[0];
rz(-1.0579146) q[1];
sx q[1];
rz(-1.3163047) q[1];
sx q[1];
rz(-2.8071075) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34821389) q[0];
sx q[0];
rz(-2.0248027) q[0];
sx q[0];
rz(0.55310849) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7203015) q[2];
sx q[2];
rz(-0.80688804) q[2];
sx q[2];
rz(0.56299201) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1002754) q[1];
sx q[1];
rz(-1.154049) q[1];
sx q[1];
rz(-2.7303425) q[1];
x q[2];
rz(0.59774996) q[3];
sx q[3];
rz(-1.808262) q[3];
sx q[3];
rz(-1.5721842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1841396) q[2];
sx q[2];
rz(-2.2440971) q[2];
sx q[2];
rz(-1.755836) q[2];
rz(-0.98006788) q[3];
sx q[3];
rz(-2.6762784) q[3];
sx q[3];
rz(0.050962713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362713) q[0];
sx q[0];
rz(-0.956981) q[0];
sx q[0];
rz(-1.3901688) q[0];
rz(-2.7888489) q[1];
sx q[1];
rz(-2.2309062) q[1];
sx q[1];
rz(0.62044755) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7095637) q[0];
sx q[0];
rz(-1.480446) q[0];
sx q[0];
rz(1.6455151) q[0];
x q[1];
rz(-3.0376126) q[2];
sx q[2];
rz(-2.1562088) q[2];
sx q[2];
rz(-2.3217161) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3603649) q[1];
sx q[1];
rz(-1.2457799) q[1];
sx q[1];
rz(-0.31045352) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6246715) q[3];
sx q[3];
rz(-0.94603387) q[3];
sx q[3];
rz(0.79510915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0188401) q[2];
sx q[2];
rz(-0.41308013) q[2];
sx q[2];
rz(-1.2472461) q[2];
rz(2.9774169) q[3];
sx q[3];
rz(-1.434606) q[3];
sx q[3];
rz(2.5119787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4267047) q[0];
sx q[0];
rz(-2.1737104) q[0];
sx q[0];
rz(-1.2870652) q[0];
rz(2.5413068) q[1];
sx q[1];
rz(-1.3515819) q[1];
sx q[1];
rz(0.90369019) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1290734) q[0];
sx q[0];
rz(-3.0154069) q[0];
sx q[0];
rz(0.55380765) q[0];
rz(-pi) q[1];
rz(0.37995423) q[2];
sx q[2];
rz(-1.6340874) q[2];
sx q[2];
rz(1.5587057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43634448) q[1];
sx q[1];
rz(-1.7503993) q[1];
sx q[1];
rz(-2.8840547) q[1];
x q[2];
rz(0.41208668) q[3];
sx q[3];
rz(-1.8672322) q[3];
sx q[3];
rz(-0.81338289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6220182) q[2];
sx q[2];
rz(-0.833424) q[2];
sx q[2];
rz(-2.3681417) q[2];
rz(-1.2889688) q[3];
sx q[3];
rz(-2.6123612) q[3];
sx q[3];
rz(-2.4066063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89161038) q[0];
sx q[0];
rz(-2.1365428) q[0];
sx q[0];
rz(2.6494359) q[0];
rz(0.53897578) q[1];
sx q[1];
rz(-1.4424126) q[1];
sx q[1];
rz(-1.0999365) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.625222) q[0];
sx q[0];
rz(-1.3110135) q[0];
sx q[0];
rz(0.35210877) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1308168) q[2];
sx q[2];
rz(-2.1824565) q[2];
sx q[2];
rz(-1.7623368) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.089842794) q[1];
sx q[1];
rz(-0.81498346) q[1];
sx q[1];
rz(-3.0137193) q[1];
rz(-pi) q[2];
rz(-0.57657974) q[3];
sx q[3];
rz(-0.55906536) q[3];
sx q[3];
rz(2.7311529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.34053549) q[2];
sx q[2];
rz(-0.33581442) q[2];
sx q[2];
rz(-2.578793) q[2];
rz(0.69989145) q[3];
sx q[3];
rz(-1.2908582) q[3];
sx q[3];
rz(0.25115299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7971147) q[0];
sx q[0];
rz(-2.4176702) q[0];
sx q[0];
rz(2.7864454) q[0];
rz(1.2190602) q[1];
sx q[1];
rz(-1.9025758) q[1];
sx q[1];
rz(-0.452279) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4773524) q[0];
sx q[0];
rz(-1.9684122) q[0];
sx q[0];
rz(-2.2009322) q[0];
rz(-0.81461774) q[2];
sx q[2];
rz(-0.53485188) q[2];
sx q[2];
rz(-1.8823106) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4411583) q[1];
sx q[1];
rz(-1.5467073) q[1];
sx q[1];
rz(1.248293) q[1];
x q[2];
rz(-1.4283435) q[3];
sx q[3];
rz(-0.94678426) q[3];
sx q[3];
rz(-0.54572661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5111115) q[2];
sx q[2];
rz(-2.6177572) q[2];
sx q[2];
rz(0.13709489) q[2];
rz(1.5591722) q[3];
sx q[3];
rz(-1.1538006) q[3];
sx q[3];
rz(-2.3649575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1906076) q[0];
sx q[0];
rz(-2.3793716) q[0];
sx q[0];
rz(-1.5451587) q[0];
rz(0.51721382) q[1];
sx q[1];
rz(-1.1511753) q[1];
sx q[1];
rz(0.98701611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1008489) q[0];
sx q[0];
rz(-0.13253875) q[0];
sx q[0];
rz(-2.9805471) q[0];
rz(-pi) q[1];
rz(2.0275063) q[2];
sx q[2];
rz(-0.65908495) q[2];
sx q[2];
rz(2.6353419) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4728022) q[1];
sx q[1];
rz(-1.3061532) q[1];
sx q[1];
rz(1.8084779) q[1];
x q[2];
rz(-1.2521947) q[3];
sx q[3];
rz(-1.4907537) q[3];
sx q[3];
rz(2.6032676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78642693) q[2];
sx q[2];
rz(-1.0309018) q[2];
sx q[2];
rz(3.1332341) q[2];
rz(2.9110294) q[3];
sx q[3];
rz(-1.9356666) q[3];
sx q[3];
rz(-1.6647388) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1281328) q[0];
sx q[0];
rz(-3.1210493) q[0];
sx q[0];
rz(1.531456) q[0];
rz(1.6186591) q[1];
sx q[1];
rz(-1.5956655) q[1];
sx q[1];
rz(-1.5787554) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42199907) q[0];
sx q[0];
rz(-2.5223456) q[0];
sx q[0];
rz(0.25281711) q[0];
rz(-0.35308102) q[2];
sx q[2];
rz(-2.7432979) q[2];
sx q[2];
rz(2.9495267) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7111499) q[1];
sx q[1];
rz(-2.4949269) q[1];
sx q[1];
rz(2.4025687) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2501026) q[3];
sx q[3];
rz(-0.9281635) q[3];
sx q[3];
rz(-1.927141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3030887) q[2];
sx q[2];
rz(-1.1469301) q[2];
sx q[2];
rz(-3.1040891) q[2];
rz(2.904902) q[3];
sx q[3];
rz(-0.37875566) q[3];
sx q[3];
rz(1.8481351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8732052) q[0];
sx q[0];
rz(-2.3516042) q[0];
sx q[0];
rz(2.068212) q[0];
rz(2.7997596) q[1];
sx q[1];
rz(-2.5624202) q[1];
sx q[1];
rz(-0.62172186) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9455915) q[0];
sx q[0];
rz(-0.47179963) q[0];
sx q[0];
rz(0.5990754) q[0];
rz(-2.4814168) q[2];
sx q[2];
rz(-1.4731506) q[2];
sx q[2];
rz(2.4843189) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7555862) q[1];
sx q[1];
rz(-1.8749798) q[1];
sx q[1];
rz(2.6261397) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6584431) q[3];
sx q[3];
rz(-1.4306127) q[3];
sx q[3];
rz(-2.6722081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6217893) q[2];
sx q[2];
rz(-0.18814627) q[2];
sx q[2];
rz(-2.5780047) q[2];
rz(-1.8300736) q[3];
sx q[3];
rz(-1.2389641) q[3];
sx q[3];
rz(-0.81609503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9594864) q[0];
sx q[0];
rz(-0.86903787) q[0];
sx q[0];
rz(2.2667789) q[0];
rz(-1.4080217) q[1];
sx q[1];
rz(-0.66649109) q[1];
sx q[1];
rz(-2.5061238) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1615026) q[0];
sx q[0];
rz(-2.3005565) q[0];
sx q[0];
rz(1.1671216) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0032015) q[2];
sx q[2];
rz(-0.48241189) q[2];
sx q[2];
rz(2.8495827) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8732152) q[1];
sx q[1];
rz(-2.4627373) q[1];
sx q[1];
rz(0.89880235) q[1];
rz(-pi) q[2];
rz(1.0069153) q[3];
sx q[3];
rz(-2.0479408) q[3];
sx q[3];
rz(1.7686219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.94329876) q[2];
sx q[2];
rz(-2.0512927) q[2];
sx q[2];
rz(2.3401006) q[2];
rz(1.21579) q[3];
sx q[3];
rz(-2.2299168) q[3];
sx q[3];
rz(-0.041778684) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6112678) q[0];
sx q[0];
rz(-1.4401191) q[0];
sx q[0];
rz(-2.8339207) q[0];
rz(1.2904185) q[1];
sx q[1];
rz(-0.94480521) q[1];
sx q[1];
rz(0.31029846) q[1];
rz(2.5038638) q[2];
sx q[2];
rz(-0.22508937) q[2];
sx q[2];
rz(-1.3300016) q[2];
rz(-1.7054059) q[3];
sx q[3];
rz(-0.91840875) q[3];
sx q[3];
rz(2.8847532) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
