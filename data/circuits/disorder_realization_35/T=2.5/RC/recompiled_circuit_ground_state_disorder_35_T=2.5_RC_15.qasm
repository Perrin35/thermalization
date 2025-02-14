OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2090981) q[0];
sx q[0];
rz(-3.0734835) q[0];
sx q[0];
rz(-0.77686247) q[0];
rz(-1.3070973) q[1];
sx q[1];
rz(-0.96944648) q[1];
sx q[1];
rz(-1.8196655) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8213733) q[0];
sx q[0];
rz(-2.2801274) q[0];
sx q[0];
rz(-0.18289645) q[0];
rz(-pi) q[1];
rz(-1.2250617) q[2];
sx q[2];
rz(-2.0421197) q[2];
sx q[2];
rz(0.87475785) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3275431) q[1];
sx q[1];
rz(-0.40990489) q[1];
sx q[1];
rz(1.6442714) q[1];
rz(1.8613937) q[3];
sx q[3];
rz(-1.3826256) q[3];
sx q[3];
rz(-3.1332603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0778568) q[2];
sx q[2];
rz(-1.7813762) q[2];
sx q[2];
rz(0.35749164) q[2];
rz(2.5743971) q[3];
sx q[3];
rz(-1.4128069) q[3];
sx q[3];
rz(0.43958694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1883063) q[0];
sx q[0];
rz(-2.8585241) q[0];
sx q[0];
rz(-1.333492) q[0];
rz(-0.4370583) q[1];
sx q[1];
rz(-2.5407365) q[1];
sx q[1];
rz(2.3016047) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6463722) q[0];
sx q[0];
rz(-0.53421181) q[0];
sx q[0];
rz(1.87508) q[0];
rz(1.7425547) q[2];
sx q[2];
rz(-2.1585625) q[2];
sx q[2];
rz(-1.5422271) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31038302) q[1];
sx q[1];
rz(-2.134517) q[1];
sx q[1];
rz(-2.4604164) q[1];
rz(-pi) q[2];
rz(-0.038957314) q[3];
sx q[3];
rz(-2.1202728) q[3];
sx q[3];
rz(-2.6770858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4240894) q[2];
sx q[2];
rz(-1.4428029) q[2];
sx q[2];
rz(1.7355512) q[2];
rz(-2.9794335) q[3];
sx q[3];
rz(-1.5806961) q[3];
sx q[3];
rz(-2.6774008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22299396) q[0];
sx q[0];
rz(-0.080568947) q[0];
sx q[0];
rz(-2.719847) q[0];
rz(-2.2987507) q[1];
sx q[1];
rz(-1.755736) q[1];
sx q[1];
rz(0.27580321) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3746412) q[0];
sx q[0];
rz(-1.8948484) q[0];
sx q[0];
rz(1.3270017) q[0];
x q[1];
rz(0.47068542) q[2];
sx q[2];
rz(-1.7690725) q[2];
sx q[2];
rz(-2.7338365) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2308064) q[1];
sx q[1];
rz(-1.8911929) q[1];
sx q[1];
rz(-1.6256534) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4546605) q[3];
sx q[3];
rz(-0.7051055) q[3];
sx q[3];
rz(0.0034611387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80321035) q[2];
sx q[2];
rz(-0.79757491) q[2];
sx q[2];
rz(-0.69022834) q[2];
rz(2.7817182) q[3];
sx q[3];
rz(-1.3709603) q[3];
sx q[3];
rz(-0.38351044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91367078) q[0];
sx q[0];
rz(-1.0424732) q[0];
sx q[0];
rz(0.87164718) q[0];
rz(-2.7323885) q[1];
sx q[1];
rz(-2.6458461) q[1];
sx q[1];
rz(-0.31235487) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.271628) q[0];
sx q[0];
rz(-1.4532491) q[0];
sx q[0];
rz(-2.6948476) q[0];
rz(-pi) q[1];
rz(-2.6957507) q[2];
sx q[2];
rz(-0.22805691) q[2];
sx q[2];
rz(-0.40846172) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4620381) q[1];
sx q[1];
rz(-1.3674539) q[1];
sx q[1];
rz(0.36851369) q[1];
rz(-pi) q[2];
rz(-2.8866037) q[3];
sx q[3];
rz(-2.4772128) q[3];
sx q[3];
rz(0.57306108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13545869) q[2];
sx q[2];
rz(-1.4486855) q[2];
sx q[2];
rz(0.82497605) q[2];
rz(-1.7168335) q[3];
sx q[3];
rz(-2.8202839) q[3];
sx q[3];
rz(-2.7062374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7216126) q[0];
sx q[0];
rz(-1.5897607) q[0];
sx q[0];
rz(2.2671674) q[0];
rz(-2.8127316) q[1];
sx q[1];
rz(-1.7560274) q[1];
sx q[1];
rz(1.0985451) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.335137) q[0];
sx q[0];
rz(-2.7493434) q[0];
sx q[0];
rz(-2.8492938) q[0];
rz(-pi) q[1];
rz(-0.83580612) q[2];
sx q[2];
rz(-1.5975738) q[2];
sx q[2];
rz(1.4876897) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85531536) q[1];
sx q[1];
rz(-0.5868656) q[1];
sx q[1];
rz(-1.0081968) q[1];
rz(-pi) q[2];
rz(-0.27733488) q[3];
sx q[3];
rz(-0.99858741) q[3];
sx q[3];
rz(-0.25857833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.141779) q[2];
sx q[2];
rz(-1.9580611) q[2];
sx q[2];
rz(2.4675274) q[2];
rz(1.4874124) q[3];
sx q[3];
rz(-1.6964648) q[3];
sx q[3];
rz(2.8580247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11033002) q[0];
sx q[0];
rz(-1.0227579) q[0];
sx q[0];
rz(1.7559825) q[0];
rz(2.022187) q[1];
sx q[1];
rz(-1.4675354) q[1];
sx q[1];
rz(-1.3853692) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.114217) q[0];
sx q[0];
rz(-1.3189335) q[0];
sx q[0];
rz(-0.92433874) q[0];
rz(-2.643061) q[2];
sx q[2];
rz(-1.4953002) q[2];
sx q[2];
rz(-0.99270644) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1277571) q[1];
sx q[1];
rz(-2.7633939) q[1];
sx q[1];
rz(-0.54658039) q[1];
x q[2];
rz(-0.22988152) q[3];
sx q[3];
rz(-0.15860367) q[3];
sx q[3];
rz(1.6386709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0783656) q[2];
sx q[2];
rz(-2.207022) q[2];
sx q[2];
rz(-2.7654977) q[2];
rz(1.1345351) q[3];
sx q[3];
rz(-2.7781656) q[3];
sx q[3];
rz(2.334972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019808708) q[0];
sx q[0];
rz(-2.5406295) q[0];
sx q[0];
rz(-0.10251775) q[0];
rz(2.5566697) q[1];
sx q[1];
rz(-2.1939317) q[1];
sx q[1];
rz(1.9580511) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4393504) q[0];
sx q[0];
rz(-2.9493414) q[0];
sx q[0];
rz(-0.73701145) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4790484) q[2];
sx q[2];
rz(-1.8383611) q[2];
sx q[2];
rz(2.2217158) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9283824) q[1];
sx q[1];
rz(-2.2021535) q[1];
sx q[1];
rz(-0.31698314) q[1];
x q[2];
rz(-0.46681301) q[3];
sx q[3];
rz(-1.8100421) q[3];
sx q[3];
rz(-1.5249263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3465603) q[2];
sx q[2];
rz(-0.9298032) q[2];
sx q[2];
rz(-2.9417876) q[2];
rz(-1.0605158) q[3];
sx q[3];
rz(-0.54148713) q[3];
sx q[3];
rz(3.0919891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.878433) q[0];
sx q[0];
rz(-1.6596154) q[0];
sx q[0];
rz(-0.31295452) q[0];
rz(2.1921659) q[1];
sx q[1];
rz(-1.3525617) q[1];
sx q[1];
rz(-1.4535646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1412555) q[0];
sx q[0];
rz(-0.20788684) q[0];
sx q[0];
rz(2.6882386) q[0];
rz(-0.61954478) q[2];
sx q[2];
rz(-1.5063553) q[2];
sx q[2];
rz(-0.10945129) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9447228) q[1];
sx q[1];
rz(-0.88366449) q[1];
sx q[1];
rz(1.2792759) q[1];
x q[2];
rz(-1.231915) q[3];
sx q[3];
rz(-1.1383071) q[3];
sx q[3];
rz(2.5755455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7877385) q[2];
sx q[2];
rz(-0.9757897) q[2];
sx q[2];
rz(1.3346416) q[2];
rz(1.8169962) q[3];
sx q[3];
rz(-0.66671222) q[3];
sx q[3];
rz(1.2142396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1052833) q[0];
sx q[0];
rz(-2.5258625) q[0];
sx q[0];
rz(2.4435254) q[0];
rz(-0.87567323) q[1];
sx q[1];
rz(-2.0798637) q[1];
sx q[1];
rz(2.7811513) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9746285) q[0];
sx q[0];
rz(-0.95050838) q[0];
sx q[0];
rz(1.0247158) q[0];
rz(2.9827732) q[2];
sx q[2];
rz(-0.38714001) q[2];
sx q[2];
rz(2.5617245) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7061758) q[1];
sx q[1];
rz(-0.98505322) q[1];
sx q[1];
rz(0.46153586) q[1];
rz(-pi) q[2];
rz(-1.7385941) q[3];
sx q[3];
rz(-2.1133964) q[3];
sx q[3];
rz(2.8436273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91161072) q[2];
sx q[2];
rz(-1.4260099) q[2];
sx q[2];
rz(0.81725517) q[2];
rz(-0.83003712) q[3];
sx q[3];
rz(-2.6006112) q[3];
sx q[3];
rz(-0.219492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-2.5068186) q[0];
sx q[0];
rz(-2.2932678) q[0];
sx q[0];
rz(3.1095374) q[0];
rz(1.3196779) q[1];
sx q[1];
rz(-1.3751043) q[1];
sx q[1];
rz(-0.65418902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20248948) q[0];
sx q[0];
rz(-1.3746885) q[0];
sx q[0];
rz(1.6926426) q[0];
rz(0.78829576) q[2];
sx q[2];
rz(-1.239778) q[2];
sx q[2];
rz(0.85110474) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54683751) q[1];
sx q[1];
rz(-1.0887301) q[1];
sx q[1];
rz(2.099311) q[1];
x q[2];
rz(-2.8579162) q[3];
sx q[3];
rz(-0.95796889) q[3];
sx q[3];
rz(-0.046663849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1303945) q[2];
sx q[2];
rz(-1.0426714) q[2];
sx q[2];
rz(1.0065669) q[2];
rz(-2.448163) q[3];
sx q[3];
rz(-1.8208241) q[3];
sx q[3];
rz(1.8600195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4257767) q[0];
sx q[0];
rz(-0.67831138) q[0];
sx q[0];
rz(-0.074445733) q[0];
rz(0.88059942) q[1];
sx q[1];
rz(-2.0279299) q[1];
sx q[1];
rz(0.50150064) q[1];
rz(2.2702552) q[2];
sx q[2];
rz(-1.9819145) q[2];
sx q[2];
rz(-2.869538) q[2];
rz(0.44742031) q[3];
sx q[3];
rz(-0.74759737) q[3];
sx q[3];
rz(2.916593) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
