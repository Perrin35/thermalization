OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.35535204) q[0];
sx q[0];
rz(-0.30682895) q[0];
sx q[0];
rz(0.29875779) q[0];
rz(-0.6660676) q[1];
sx q[1];
rz(-0.63036418) q[1];
sx q[1];
rz(1.4730374) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9763223) q[0];
sx q[0];
rz(-2.3403694) q[0];
sx q[0];
rz(-3.0972605) q[0];
rz(-pi) q[1];
x q[1];
rz(0.017919964) q[2];
sx q[2];
rz(-1.8715053) q[2];
sx q[2];
rz(1.9422046) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7466429) q[1];
sx q[1];
rz(-0.87374765) q[1];
sx q[1];
rz(0.66956981) q[1];
rz(-pi) q[2];
rz(-2.8979635) q[3];
sx q[3];
rz(-2.8578161) q[3];
sx q[3];
rz(2.711913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.95639688) q[2];
sx q[2];
rz(-0.38072017) q[2];
sx q[2];
rz(-1.8308651) q[2];
rz(0.091536097) q[3];
sx q[3];
rz(-2.5444578) q[3];
sx q[3];
rz(-2.6105647) q[3];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0700584) q[0];
sx q[0];
rz(-1.0907084) q[0];
sx q[0];
rz(-2.6614905) q[0];
rz(-1.9942888) q[1];
sx q[1];
rz(-2.5105748) q[1];
sx q[1];
rz(0.70153418) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.498256) q[0];
sx q[0];
rz(-2.069888) q[0];
sx q[0];
rz(-2.3168111) q[0];
rz(0.417493) q[2];
sx q[2];
rz(-2.001363) q[2];
sx q[2];
rz(-0.21829641) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5614723) q[1];
sx q[1];
rz(-1.1271203) q[1];
sx q[1];
rz(0.89718141) q[1];
x q[2];
rz(0.63752709) q[3];
sx q[3];
rz(-1.8211604) q[3];
sx q[3];
rz(-1.375578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7445765) q[2];
sx q[2];
rz(-1.959356) q[2];
sx q[2];
rz(-1.5783295) q[2];
rz(-2.8619316) q[3];
sx q[3];
rz(-0.71664387) q[3];
sx q[3];
rz(-0.40951148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-0.79769832) q[0];
sx q[0];
rz(-2.7293623) q[0];
sx q[0];
rz(1.718234) q[0];
rz(-3.0176945) q[1];
sx q[1];
rz(-0.84404498) q[1];
sx q[1];
rz(-1.5843676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9950969) q[0];
sx q[0];
rz(-2.2687758) q[0];
sx q[0];
rz(0.3905889) q[0];
rz(0.84752797) q[2];
sx q[2];
rz(-1.28994) q[2];
sx q[2];
rz(-2.6658863) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.023097087) q[1];
sx q[1];
rz(-1.4916991) q[1];
sx q[1];
rz(1.817607) q[1];
rz(-pi) q[2];
rz(-2.4067114) q[3];
sx q[3];
rz(-2.2874444) q[3];
sx q[3];
rz(-1.7295966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.065217) q[2];
sx q[2];
rz(-2.4880444) q[2];
sx q[2];
rz(-0.93488133) q[2];
rz(0.30385083) q[3];
sx q[3];
rz(-0.6128208) q[3];
sx q[3];
rz(2.2936308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76871753) q[0];
sx q[0];
rz(-0.77943742) q[0];
sx q[0];
rz(-0.26707643) q[0];
rz(-3.0067387) q[1];
sx q[1];
rz(-0.57661533) q[1];
sx q[1];
rz(2.3042302) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7202111) q[0];
sx q[0];
rz(-1.0307587) q[0];
sx q[0];
rz(-2.8296986) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5576511) q[2];
sx q[2];
rz(-2.2211233) q[2];
sx q[2];
rz(0.038829858) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0935614) q[1];
sx q[1];
rz(-1.8037027) q[1];
sx q[1];
rz(-1.3711817) q[1];
x q[2];
rz(0.046810026) q[3];
sx q[3];
rz(-1.6322281) q[3];
sx q[3];
rz(-2.2991179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5460633) q[2];
sx q[2];
rz(-0.17543051) q[2];
sx q[2];
rz(2.9626633) q[2];
rz(1.1692125) q[3];
sx q[3];
rz(-1.7763276) q[3];
sx q[3];
rz(-0.21807142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1201852) q[0];
sx q[0];
rz(-2.6299801) q[0];
sx q[0];
rz(-0.87274337) q[0];
rz(-1.1574289) q[1];
sx q[1];
rz(-2.2667784) q[1];
sx q[1];
rz(-0.016955888) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8166148) q[0];
sx q[0];
rz(-0.8362174) q[0];
sx q[0];
rz(-1.5543544) q[0];
rz(-pi) q[1];
rz(-2.5418607) q[2];
sx q[2];
rz(-1.5625192) q[2];
sx q[2];
rz(-0.042009609) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.03093623) q[1];
sx q[1];
rz(-0.82582322) q[1];
sx q[1];
rz(0.40542545) q[1];
rz(-1.2797375) q[3];
sx q[3];
rz(-1.3123056) q[3];
sx q[3];
rz(2.2973826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0852647) q[2];
sx q[2];
rz(-1.7076098) q[2];
sx q[2];
rz(0.31923527) q[2];
rz(-0.6790092) q[3];
sx q[3];
rz(-1.7947936) q[3];
sx q[3];
rz(2.3398248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7769258) q[0];
sx q[0];
rz(-2.8209782) q[0];
sx q[0];
rz(-2.4989682) q[0];
rz(1.7721666) q[1];
sx q[1];
rz(-0.6232999) q[1];
sx q[1];
rz(-2.1646037) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95648051) q[0];
sx q[0];
rz(-1.6688884) q[0];
sx q[0];
rz(-3.0177659) q[0];
x q[1];
rz(-2.0989292) q[2];
sx q[2];
rz(-2.0181877) q[2];
sx q[2];
rz(1.3716979) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9639324) q[1];
sx q[1];
rz(-2.2626082) q[1];
sx q[1];
rz(3.1236137) q[1];
x q[2];
rz(-0.98501916) q[3];
sx q[3];
rz(-1.5259229) q[3];
sx q[3];
rz(0.45811996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52648181) q[2];
sx q[2];
rz(-1.4075764) q[2];
sx q[2];
rz(1.9132445) q[2];
rz(0.86722106) q[3];
sx q[3];
rz(-2.7286178) q[3];
sx q[3];
rz(-2.373608) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41384554) q[0];
sx q[0];
rz(-2.9242046) q[0];
sx q[0];
rz(0.17284285) q[0];
rz(0.8051644) q[1];
sx q[1];
rz(-0.54904896) q[1];
sx q[1];
rz(3.0056312) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522536) q[0];
sx q[0];
rz(-2.2704008) q[0];
sx q[0];
rz(0.1564349) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4589086) q[2];
sx q[2];
rz(-0.75190836) q[2];
sx q[2];
rz(-0.4900527) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.48366212) q[1];
sx q[1];
rz(-1.6854981) q[1];
sx q[1];
rz(1.5899659) q[1];
rz(-pi) q[2];
rz(0.32973955) q[3];
sx q[3];
rz(-1.6981594) q[3];
sx q[3];
rz(-1.691526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.92676306) q[2];
sx q[2];
rz(-1.5462993) q[2];
sx q[2];
rz(-0.1845486) q[2];
rz(-0.18937011) q[3];
sx q[3];
rz(-0.53818494) q[3];
sx q[3];
rz(-0.93809938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9823343) q[0];
sx q[0];
rz(-2.8026447) q[0];
sx q[0];
rz(0.92639297) q[0];
rz(3.1023846) q[1];
sx q[1];
rz(-0.67752939) q[1];
sx q[1];
rz(1.0367941) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3711972) q[0];
sx q[0];
rz(-1.3950893) q[0];
sx q[0];
rz(1.2716588) q[0];
rz(-pi) q[1];
rz(0.43208684) q[2];
sx q[2];
rz(-1.7361904) q[2];
sx q[2];
rz(1.337932) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7651789) q[1];
sx q[1];
rz(-2.661099) q[1];
sx q[1];
rz(-2.5960467) q[1];
x q[2];
rz(1.3831656) q[3];
sx q[3];
rz(-1.1867282) q[3];
sx q[3];
rz(1.7373067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.89256531) q[2];
sx q[2];
rz(-2.0165063) q[2];
sx q[2];
rz(0.79891515) q[2];
rz(0.51698452) q[3];
sx q[3];
rz(-0.40701443) q[3];
sx q[3];
rz(-0.5504722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7901881) q[0];
sx q[0];
rz(-3.1365972) q[0];
sx q[0];
rz(1.5193526) q[0];
rz(1.5478569) q[1];
sx q[1];
rz(-2.3212815) q[1];
sx q[1];
rz(2.4511852) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.658118) q[0];
sx q[0];
rz(-0.42040792) q[0];
sx q[0];
rz(-0.57256521) q[0];
x q[1];
rz(2.3423455) q[2];
sx q[2];
rz(-1.4991781) q[2];
sx q[2];
rz(1.3992753) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4286297) q[1];
sx q[1];
rz(-2.3979514) q[1];
sx q[1];
rz(-0.78674591) q[1];
rz(-pi) q[2];
rz(0.14491187) q[3];
sx q[3];
rz(-2.1050801) q[3];
sx q[3];
rz(2.2458959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4385628) q[2];
sx q[2];
rz(-2.5327693) q[2];
sx q[2];
rz(-0.969886) q[2];
rz(0.25740933) q[3];
sx q[3];
rz(-0.22203797) q[3];
sx q[3];
rz(-1.359587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53512204) q[0];
sx q[0];
rz(-0.73000014) q[0];
sx q[0];
rz(0.085513376) q[0];
rz(2.8657148) q[1];
sx q[1];
rz(-0.62092263) q[1];
sx q[1];
rz(-1.9524908) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71578854) q[0];
sx q[0];
rz(-2.5809408) q[0];
sx q[0];
rz(2.406757) q[0];
rz(-pi) q[1];
rz(3.1136572) q[2];
sx q[2];
rz(-2.8314674) q[2];
sx q[2];
rz(-1.5903697) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4170749) q[1];
sx q[1];
rz(-0.41890422) q[1];
sx q[1];
rz(-1.8663097) q[1];
rz(-pi) q[2];
rz(0.25190763) q[3];
sx q[3];
rz(-0.35874507) q[3];
sx q[3];
rz(-2.9210966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6659866) q[2];
sx q[2];
rz(-0.87916547) q[2];
sx q[2];
rz(2.6574668) q[2];
rz(-0.13949805) q[3];
sx q[3];
rz(-0.77553427) q[3];
sx q[3];
rz(0.16105306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.465268) q[0];
sx q[0];
rz(-1.9857255) q[0];
sx q[0];
rz(2.3406512) q[0];
rz(1.4576661) q[1];
sx q[1];
rz(-1.6438345) q[1];
sx q[1];
rz(1.8297292) q[1];
rz(0.2659145) q[2];
sx q[2];
rz(-2.4752046) q[2];
sx q[2];
rz(2.2023329) q[2];
rz(-0.030293754) q[3];
sx q[3];
rz(-2.3843063) q[3];
sx q[3];
rz(1.6142308) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
