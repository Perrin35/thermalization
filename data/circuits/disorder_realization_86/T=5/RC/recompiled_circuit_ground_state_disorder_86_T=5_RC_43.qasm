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
rz(2.3613858) q[1];
sx q[1];
rz(-2.0145887) q[1];
sx q[1];
rz(-2.3957774) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1231664) q[0];
sx q[0];
rz(-1.9115051) q[0];
sx q[0];
rz(-0.19827224) q[0];
rz(1.4100513) q[2];
sx q[2];
rz(-1.5622592) q[2];
sx q[2];
rz(1.05913) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76872003) q[1];
sx q[1];
rz(-0.24133397) q[1];
sx q[1];
rz(-1.0670426) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7379825) q[3];
sx q[3];
rz(-2.0847528) q[3];
sx q[3];
rz(1.5486167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9790393) q[2];
sx q[2];
rz(-2.5472842) q[2];
sx q[2];
rz(-1.7621367) q[2];
rz(-0.72921324) q[3];
sx q[3];
rz(-0.76494923) q[3];
sx q[3];
rz(-2.2123857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58296975) q[0];
sx q[0];
rz(-1.0635149) q[0];
sx q[0];
rz(-2.9003918) q[0];
rz(-2.8855715) q[1];
sx q[1];
rz(-1.0216917) q[1];
sx q[1];
rz(-0.78114885) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.985551) q[0];
sx q[0];
rz(-0.80901481) q[0];
sx q[0];
rz(-2.3963905) q[0];
x q[1];
rz(-0.56480572) q[2];
sx q[2];
rz(-1.5283268) q[2];
sx q[2];
rz(0.04893411) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9731739) q[1];
sx q[1];
rz(-0.89514625) q[1];
sx q[1];
rz(-0.78898095) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2595176) q[3];
sx q[3];
rz(-2.9972509) q[3];
sx q[3];
rz(0.33931574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.29391924) q[2];
sx q[2];
rz(-1.5519698) q[2];
sx q[2];
rz(1.5839362) q[2];
rz(0.0047575792) q[3];
sx q[3];
rz(-0.59499756) q[3];
sx q[3];
rz(-0.34576542) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7541589) q[0];
sx q[0];
rz(-1.4492946) q[0];
sx q[0];
rz(2.9505728) q[0];
rz(2.7905131) q[1];
sx q[1];
rz(-0.43126884) q[1];
sx q[1];
rz(1.4494337) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75101484) q[0];
sx q[0];
rz(-1.6216941) q[0];
sx q[0];
rz(-2.8021332) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1507785) q[2];
sx q[2];
rz(-1.1598829) q[2];
sx q[2];
rz(2.0793629) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8061132) q[1];
sx q[1];
rz(-1.3054784) q[1];
sx q[1];
rz(0.47228864) q[1];
rz(-pi) q[2];
rz(-1.2292858) q[3];
sx q[3];
rz(-1.0060898) q[3];
sx q[3];
rz(-2.963394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8423975) q[2];
sx q[2];
rz(-1.7123875) q[2];
sx q[2];
rz(-3.0136285) q[2];
rz(-1.9757102) q[3];
sx q[3];
rz(-2.2539625) q[3];
sx q[3];
rz(2.8036346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8562451) q[0];
sx q[0];
rz(-1.0839533) q[0];
sx q[0];
rz(1.6612843) q[0];
rz(-0.081143204) q[1];
sx q[1];
rz(-1.0973884) q[1];
sx q[1];
rz(2.2611484) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5187982) q[0];
sx q[0];
rz(-1.5642903) q[0];
sx q[0];
rz(3.1168808) q[0];
rz(1.6398076) q[2];
sx q[2];
rz(-2.5216148) q[2];
sx q[2];
rz(2.2295956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8740011) q[1];
sx q[1];
rz(-2.1028215) q[1];
sx q[1];
rz(1.6997972) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5577291) q[3];
sx q[3];
rz(-1.5506049) q[3];
sx q[3];
rz(2.9021183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.60260281) q[2];
sx q[2];
rz(-2.4675641) q[2];
sx q[2];
rz(1.7146141) q[2];
rz(1.1566628) q[3];
sx q[3];
rz(-1.7159228) q[3];
sx q[3];
rz(-0.15779933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(2.2222774) q[0];
sx q[0];
rz(-0.8150402) q[0];
sx q[0];
rz(0.94006938) q[0];
rz(2.6433511) q[1];
sx q[1];
rz(-2.2254641) q[1];
sx q[1];
rz(1.1246276) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2021411) q[0];
sx q[0];
rz(-1.6103585) q[0];
sx q[0];
rz(2.2742154) q[0];
rz(-pi) q[1];
rz(-0.20467918) q[2];
sx q[2];
rz(-1.9726557) q[2];
sx q[2];
rz(3.0167442) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4044721) q[1];
sx q[1];
rz(-2.0029161) q[1];
sx q[1];
rz(2.4678556) q[1];
rz(-pi) q[2];
rz(-0.78557555) q[3];
sx q[3];
rz(-1.5858558) q[3];
sx q[3];
rz(-2.756898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70790946) q[2];
sx q[2];
rz(-2.3072672) q[2];
sx q[2];
rz(-1.4617807) q[2];
rz(-2.8953569) q[3];
sx q[3];
rz(-0.88053954) q[3];
sx q[3];
rz(-2.06854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4435302) q[0];
sx q[0];
rz(-1.9140697) q[0];
sx q[0];
rz(-0.45502934) q[0];
rz(-2.9235234) q[1];
sx q[1];
rz(-2.2360305) q[1];
sx q[1];
rz(-2.6630482) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25898257) q[0];
sx q[0];
rz(-0.6713258) q[0];
sx q[0];
rz(0.15780003) q[0];
rz(-pi) q[1];
rz(-0.63963525) q[2];
sx q[2];
rz(-1.2811077) q[2];
sx q[2];
rz(-1.7997431) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8763833) q[1];
sx q[1];
rz(-2.5647616) q[1];
sx q[1];
rz(0.93872197) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6932186) q[3];
sx q[3];
rz(-2.3230334) q[3];
sx q[3];
rz(2.4002176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9559418) q[2];
sx q[2];
rz(-1.4943244) q[2];
sx q[2];
rz(-2.4556124) q[2];
rz(1.8020804) q[3];
sx q[3];
rz(-1.9728262) q[3];
sx q[3];
rz(-1.3919938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5550391) q[0];
sx q[0];
rz(-1.3864484) q[0];
sx q[0];
rz(2.6626124) q[0];
rz(-0.85887495) q[1];
sx q[1];
rz(-0.54194599) q[1];
sx q[1];
rz(0.10471334) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59706748) q[0];
sx q[0];
rz(-1.1865718) q[0];
sx q[0];
rz(-1.3659507) q[0];
x q[1];
rz(-0.62868406) q[2];
sx q[2];
rz(-0.57839823) q[2];
sx q[2];
rz(-2.4378928) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.59440246) q[1];
sx q[1];
rz(-2.5195751) q[1];
sx q[1];
rz(-1.1941431) q[1];
x q[2];
rz(1.4604578) q[3];
sx q[3];
rz(-0.48052999) q[3];
sx q[3];
rz(-0.21267173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0844448) q[2];
sx q[2];
rz(-1.8859325) q[2];
sx q[2];
rz(-0.88195938) q[2];
rz(-0.070934892) q[3];
sx q[3];
rz(-1.2627914) q[3];
sx q[3];
rz(2.1766263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98081812) q[0];
sx q[0];
rz(-2.6698298) q[0];
sx q[0];
rz(1.8881352) q[0];
rz(-2.4313633) q[1];
sx q[1];
rz(-1.9435147) q[1];
sx q[1];
rz(2.0517147) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1313497) q[0];
sx q[0];
rz(-2.2585218) q[0];
sx q[0];
rz(-1.0279845) q[0];
rz(-pi) q[1];
rz(0.44909018) q[2];
sx q[2];
rz(-1.8829926) q[2];
sx q[2];
rz(-2.5132266) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3115499) q[1];
sx q[1];
rz(-1.122992) q[1];
sx q[1];
rz(2.2954659) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5651781) q[3];
sx q[3];
rz(-1.6987598) q[3];
sx q[3];
rz(2.0258249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.39500427) q[2];
sx q[2];
rz(-2.2825664) q[2];
sx q[2];
rz(-2.11917) q[2];
rz(-0.59631452) q[3];
sx q[3];
rz(-0.78323451) q[3];
sx q[3];
rz(-0.35308009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8280867) q[0];
sx q[0];
rz(-2.6441898) q[0];
sx q[0];
rz(1.585438) q[0];
rz(1.4900788) q[1];
sx q[1];
rz(-0.45526344) q[1];
sx q[1];
rz(0.99686399) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23295924) q[0];
sx q[0];
rz(-1.5002439) q[0];
sx q[0];
rz(2.7951827) q[0];
rz(-1.5794706) q[2];
sx q[2];
rz(-0.18520912) q[2];
sx q[2];
rz(-0.26488525) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1151859) q[1];
sx q[1];
rz(-0.86167012) q[1];
sx q[1];
rz(-0.89964189) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97932696) q[3];
sx q[3];
rz(-1.1445657) q[3];
sx q[3];
rz(-1.1104079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.31522125) q[2];
sx q[2];
rz(-0.73306495) q[2];
sx q[2];
rz(-0.23615393) q[2];
rz(-0.30596966) q[3];
sx q[3];
rz(-1.1924084) q[3];
sx q[3];
rz(-1.6570305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21208256) q[0];
sx q[0];
rz(-0.65123737) q[0];
sx q[0];
rz(-2.4834852) q[0];
rz(-1.8923538) q[1];
sx q[1];
rz(-0.58964261) q[1];
sx q[1];
rz(2.6499937) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8458873) q[0];
sx q[0];
rz(-2.108556) q[0];
sx q[0];
rz(2.0989492) q[0];
x q[1];
rz(-0.34296918) q[2];
sx q[2];
rz(-0.59294564) q[2];
sx q[2];
rz(0.85747257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1695849) q[1];
sx q[1];
rz(-1.8465202) q[1];
sx q[1];
rz(1.5610525) q[1];
rz(1.4171353) q[3];
sx q[3];
rz(-2.6479122) q[3];
sx q[3];
rz(0.069531893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7813985) q[2];
sx q[2];
rz(-0.34803826) q[2];
sx q[2];
rz(1.4813102) q[2];
rz(0.71634746) q[3];
sx q[3];
rz(-2.4951388) q[3];
sx q[3];
rz(0.20814482) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.70474456) q[2];
sx q[2];
rz(-0.8316883) q[2];
sx q[2];
rz(-1.4744454) q[2];
rz(-0.17070676) q[3];
sx q[3];
rz(-0.89912631) q[3];
sx q[3];
rz(2.5691433) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
