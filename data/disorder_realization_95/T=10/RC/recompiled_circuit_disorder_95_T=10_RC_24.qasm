OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.99958217) q[0];
sx q[0];
rz(5.2810623) q[0];
sx q[0];
rz(5.3856344) q[0];
rz(2.907213) q[1];
sx q[1];
rz(-2.8657764) q[1];
sx q[1];
rz(-2.0770567) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8911154) q[0];
sx q[0];
rz(-0.63397898) q[0];
sx q[0];
rz(1.6888213) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35383309) q[2];
sx q[2];
rz(-0.96828038) q[2];
sx q[2];
rz(-1.8941855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7251016) q[1];
sx q[1];
rz(-0.93826586) q[1];
sx q[1];
rz(-1.3593258) q[1];
rz(-pi) q[2];
rz(1.3863871) q[3];
sx q[3];
rz(-2.3204436) q[3];
sx q[3];
rz(-2.8781995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66427461) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(-0.73195362) q[2];
rz(-0.96015635) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(-1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17523781) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(2.2251341) q[0];
rz(2.6610999) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(0.8786456) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4247596) q[0];
sx q[0];
rz(-0.75100198) q[0];
sx q[0];
rz(-2.1483833) q[0];
rz(2.7777113) q[2];
sx q[2];
rz(-0.65081396) q[2];
sx q[2];
rz(1.770307) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52777427) q[1];
sx q[1];
rz(-1.6442181) q[1];
sx q[1];
rz(-2.7223177) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8804178) q[3];
sx q[3];
rz(-1.505758) q[3];
sx q[3];
rz(1.0950973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8615222) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(-0.64727616) q[2];
rz(2.9679126) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(-0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6140401) q[0];
sx q[0];
rz(-2.5019167) q[0];
sx q[0];
rz(2.8955984) q[0];
rz(1.4100769) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(-2.0203967) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5593528) q[0];
sx q[0];
rz(-1.2057349) q[0];
sx q[0];
rz(0.75612005) q[0];
x q[1];
rz(-1.9636743) q[2];
sx q[2];
rz(-1.5117206) q[2];
sx q[2];
rz(1.7977561) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1998636) q[1];
sx q[1];
rz(-0.45048303) q[1];
sx q[1];
rz(0.73168879) q[1];
rz(-pi) q[2];
rz(1.4533914) q[3];
sx q[3];
rz(-1.0798287) q[3];
sx q[3];
rz(1.1743869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7401509) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(-2.1901954) q[2];
rz(-0.65008632) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(-2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51405108) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(-0.78805584) q[0];
rz(2.0846941) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(0.11638164) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43198904) q[0];
sx q[0];
rz(-2.5076712) q[0];
sx q[0];
rz(-3.0984127) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94159796) q[2];
sx q[2];
rz(-2.5430352) q[2];
sx q[2];
rz(-1.4231921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0105646) q[1];
sx q[1];
rz(-0.66337913) q[1];
sx q[1];
rz(0.69372155) q[1];
rz(3.0451123) q[3];
sx q[3];
rz(-2.6498142) q[3];
sx q[3];
rz(1.552856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5300166) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(-0.25203618) q[2];
rz(0.37825545) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(-2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(0.53260032) q[0];
rz(-1.4415007) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(1.9929569) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.900433) q[0];
sx q[0];
rz(-0.94928375) q[0];
sx q[0];
rz(-3.050699) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2448036) q[2];
sx q[2];
rz(-1.9325581) q[2];
sx q[2];
rz(2.0828431) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8685638) q[1];
sx q[1];
rz(-1.3820573) q[1];
sx q[1];
rz(2.7815458) q[1];
rz(-pi) q[2];
rz(-2.4162021) q[3];
sx q[3];
rz(-1.2687506) q[3];
sx q[3];
rz(-2.6186752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0107515) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(-2.1095236) q[2];
rz(-0.71470913) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(-3.1177974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(-0.39598879) q[0];
rz(-1.4453325) q[1];
sx q[1];
rz(-1.6991801) q[1];
sx q[1];
rz(2.9352303) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8800031) q[0];
sx q[0];
rz(-2.3149009) q[0];
sx q[0];
rz(-1.4522626) q[0];
rz(-0.66531078) q[2];
sx q[2];
rz(-2.1449001) q[2];
sx q[2];
rz(2.2523508) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.57950912) q[1];
sx q[1];
rz(-1.675203) q[1];
sx q[1];
rz(0.097192055) q[1];
rz(1.0422802) q[3];
sx q[3];
rz(-2.7665666) q[3];
sx q[3];
rz(3.022775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6017194) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(0.27077857) q[2];
rz(-2.3932636) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(-2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2449743) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(-0.57624972) q[0];
rz(-0.17310625) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(-1.9304088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21748397) q[0];
sx q[0];
rz(-0.78290126) q[0];
sx q[0];
rz(-0.80362513) q[0];
rz(-pi) q[1];
rz(-1.7089825) q[2];
sx q[2];
rz(-2.0565363) q[2];
sx q[2];
rz(1.1952343) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32528977) q[1];
sx q[1];
rz(-1.6386697) q[1];
sx q[1];
rz(-1.7273278) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1738449) q[3];
sx q[3];
rz(-1.3021886) q[3];
sx q[3];
rz(2.8318162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(0.024519196) q[2];
rz(-2.426614) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(-1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.0054758469) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(2.1210282) q[0];
rz(0.15469805) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(-1.221009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97604254) q[0];
sx q[0];
rz(-1.3001469) q[0];
sx q[0];
rz(0.793215) q[0];
rz(-0.74322015) q[2];
sx q[2];
rz(-2.9033702) q[2];
sx q[2];
rz(-1.5526349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.78302661) q[1];
sx q[1];
rz(-1.4604124) q[1];
sx q[1];
rz(-3.0049938) q[1];
x q[2];
rz(0.33631781) q[3];
sx q[3];
rz(-0.85018966) q[3];
sx q[3];
rz(-1.99828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28875479) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(0.70518804) q[2];
rz(-0.60338902) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(-1.5195297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.110638) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(-3.0045793) q[0];
rz(0.6048454) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(3.0158214) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9916358) q[0];
sx q[0];
rz(-1.7466674) q[0];
sx q[0];
rz(1.3791023) q[0];
rz(-pi) q[1];
rz(-2.8081886) q[2];
sx q[2];
rz(-1.7575022) q[2];
sx q[2];
rz(-2.8508027) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8217433) q[1];
sx q[1];
rz(-0.84801596) q[1];
sx q[1];
rz(-0.38139947) q[1];
x q[2];
rz(2.845876) q[3];
sx q[3];
rz(-1.7641281) q[3];
sx q[3];
rz(-0.28702345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13828364) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(-0.70927817) q[2];
rz(0.62018958) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(-0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(-1.2003157) q[0];
rz(2.8740846) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(-1.1402003) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70779078) q[0];
sx q[0];
rz(-1.8330488) q[0];
sx q[0];
rz(3.0387525) q[0];
rz(-pi) q[1];
rz(-1.1085547) q[2];
sx q[2];
rz(-1.3369563) q[2];
sx q[2];
rz(-3.0548981) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.66265857) q[1];
sx q[1];
rz(-1.3824029) q[1];
sx q[1];
rz(1.6608095) q[1];
rz(-1.2530008) q[3];
sx q[3];
rz(-2.3386049) q[3];
sx q[3];
rz(-3.0190937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(2.2429788) q[2];
rz(3.1344154) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(-1.7858508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89467775) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(2.6208411) q[1];
sx q[1];
rz(-1.755935) q[1];
sx q[1];
rz(-1.4204949) q[1];
rz(0.76277914) q[2];
sx q[2];
rz(-2.503958) q[2];
sx q[2];
rz(-2.5867953) q[2];
rz(0.58017147) q[3];
sx q[3];
rz(-0.67065722) q[3];
sx q[3];
rz(-1.8967659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
