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
rz(0.41843721) q[0];
sx q[0];
rz(5.3199407) q[0];
sx q[0];
rz(9.6286019) q[0];
rz(1.0526429) q[1];
sx q[1];
rz(3.9349603) q[1];
sx q[1];
rz(10.959672) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14102916) q[0];
sx q[0];
rz(-1.805843) q[0];
sx q[0];
rz(-2.0268593) q[0];
rz(-1.1821177) q[2];
sx q[2];
rz(-0.76180327) q[2];
sx q[2];
rz(2.3488059) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2633367) q[1];
sx q[1];
rz(-1.9227826) q[1];
sx q[1];
rz(-0.89603591) q[1];
rz(-0.58310572) q[3];
sx q[3];
rz(-2.8994377) q[3];
sx q[3];
rz(0.21447578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7242929) q[2];
sx q[2];
rz(-0.55880004) q[2];
sx q[2];
rz(-2.144045) q[2];
rz(0.7558465) q[3];
sx q[3];
rz(-1.7852802) q[3];
sx q[3];
rz(2.1006179) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7928829) q[0];
sx q[0];
rz(-0.097276874) q[0];
sx q[0];
rz(-1.8541699) q[0];
rz(-0.48311326) q[1];
sx q[1];
rz(-1.0419507) q[1];
sx q[1];
rz(-1.9226673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3588803) q[0];
sx q[0];
rz(-2.8887889) q[0];
sx q[0];
rz(2.2532399) q[0];
rz(-1.7886247) q[2];
sx q[2];
rz(-2.2309897) q[2];
sx q[2];
rz(-0.24581395) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2698631) q[1];
sx q[1];
rz(-2.8262889) q[1];
sx q[1];
rz(-1.0878776) q[1];
rz(1.1514436) q[3];
sx q[3];
rz(-0.73691982) q[3];
sx q[3];
rz(-1.1518971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7634742) q[2];
sx q[2];
rz(-0.066630445) q[2];
sx q[2];
rz(2.5197869) q[2];
rz(-0.63831896) q[3];
sx q[3];
rz(-0.74898762) q[3];
sx q[3];
rz(2.2973255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8353552) q[0];
sx q[0];
rz(-2.4995646) q[0];
sx q[0];
rz(0.21632347) q[0];
rz(-2.5430039) q[1];
sx q[1];
rz(-1.8363771) q[1];
sx q[1];
rz(-2.6187706) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.030634681) q[0];
sx q[0];
rz(-1.8930264) q[0];
sx q[0];
rz(-2.8186174) q[0];
rz(-0.36575138) q[2];
sx q[2];
rz(-2.7257082) q[2];
sx q[2];
rz(1.5378086) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.712341) q[1];
sx q[1];
rz(-2.0852226) q[1];
sx q[1];
rz(1.3144668) q[1];
x q[2];
rz(0.82562311) q[3];
sx q[3];
rz(-2.7047727) q[3];
sx q[3];
rz(0.88479155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.181695) q[2];
sx q[2];
rz(-1.2744224) q[2];
sx q[2];
rz(-1.6240906) q[2];
rz(2.6767139) q[3];
sx q[3];
rz(-2.2113694) q[3];
sx q[3];
rz(-1.5215065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14438039) q[0];
sx q[0];
rz(-0.388044) q[0];
sx q[0];
rz(1.602518) q[0];
rz(-2.1082361) q[1];
sx q[1];
rz(-0.76834232) q[1];
sx q[1];
rz(1.296952) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1259114) q[0];
sx q[0];
rz(-0.75443479) q[0];
sx q[0];
rz(1.8777281) q[0];
rz(-pi) q[1];
rz(-1.4960852) q[2];
sx q[2];
rz(-2.4444408) q[2];
sx q[2];
rz(-2.8139204) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.32758157) q[1];
sx q[1];
rz(-1.1594698) q[1];
sx q[1];
rz(-0.41366215) q[1];
rz(-pi) q[2];
rz(-2.6565032) q[3];
sx q[3];
rz(-0.61044932) q[3];
sx q[3];
rz(1.7318673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3492655) q[2];
sx q[2];
rz(-0.63750625) q[2];
sx q[2];
rz(2.0959334) q[2];
rz(0.21444923) q[3];
sx q[3];
rz(-0.72434536) q[3];
sx q[3];
rz(-1.9946056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3359208) q[0];
sx q[0];
rz(-1.8593973) q[0];
sx q[0];
rz(-2.6222498) q[0];
rz(-0.94201159) q[1];
sx q[1];
rz(-2.8603034) q[1];
sx q[1];
rz(-1.5325783) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62270861) q[0];
sx q[0];
rz(-1.9340589) q[0];
sx q[0];
rz(-0.69464442) q[0];
rz(-pi) q[1];
rz(-0.66560676) q[2];
sx q[2];
rz(-2.9133248) q[2];
sx q[2];
rz(2.1537202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0377887) q[1];
sx q[1];
rz(-2.6409864) q[1];
sx q[1];
rz(-1.5816825) q[1];
rz(-pi) q[2];
rz(0.10974613) q[3];
sx q[3];
rz(-1.6648431) q[3];
sx q[3];
rz(1.4274244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.64451009) q[2];
sx q[2];
rz(-0.77748674) q[2];
sx q[2];
rz(-0.28437781) q[2];
rz(-2.583368) q[3];
sx q[3];
rz(-1.3642718) q[3];
sx q[3];
rz(-2.8936581) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068950653) q[0];
sx q[0];
rz(-0.41256368) q[0];
sx q[0];
rz(-0.64055842) q[0];
rz(1.6259646) q[1];
sx q[1];
rz(-1.1311572) q[1];
sx q[1];
rz(-0.21387771) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.36422) q[0];
sx q[0];
rz(-1.3427444) q[0];
sx q[0];
rz(1.8814726) q[0];
rz(-pi) q[1];
rz(-1.1919695) q[2];
sx q[2];
rz(-0.86244273) q[2];
sx q[2];
rz(1.5713991) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2474413) q[1];
sx q[1];
rz(-0.94235984) q[1];
sx q[1];
rz(3.0559866) q[1];
x q[2];
rz(3.0701105) q[3];
sx q[3];
rz(-2.0897647) q[3];
sx q[3];
rz(-1.9148358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3512257) q[2];
sx q[2];
rz(-0.48206097) q[2];
sx q[2];
rz(-0.17769979) q[2];
rz(1.7443582) q[3];
sx q[3];
rz(-1.1763562) q[3];
sx q[3];
rz(2.7241404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0419256) q[0];
sx q[0];
rz(-2.5302027) q[0];
sx q[0];
rz(-1.1055111) q[0];
rz(0.28327495) q[1];
sx q[1];
rz(-2.6073644) q[1];
sx q[1];
rz(1.4615321) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4257129) q[0];
sx q[0];
rz(-1.0277896) q[0];
sx q[0];
rz(2.3433861) q[0];
rz(-pi) q[1];
rz(0.92723691) q[2];
sx q[2];
rz(-0.74802784) q[2];
sx q[2];
rz(-3.1035556) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0431598) q[1];
sx q[1];
rz(-2.5665356) q[1];
sx q[1];
rz(-2.9270323) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8959778) q[3];
sx q[3];
rz(-1.0187314) q[3];
sx q[3];
rz(-2.0887041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0561515) q[2];
sx q[2];
rz(-1.5458919) q[2];
sx q[2];
rz(-0.20480569) q[2];
rz(1.0497931) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(-0.36627305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7751223) q[0];
sx q[0];
rz(-2.6752052) q[0];
sx q[0];
rz(-3.108016) q[0];
rz(1.4749984) q[1];
sx q[1];
rz(-2.3971403) q[1];
sx q[1];
rz(1.144145) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1468723) q[0];
sx q[0];
rz(-1.084274) q[0];
sx q[0];
rz(-1.542132) q[0];
rz(-1.894925) q[2];
sx q[2];
rz(-1.6648544) q[2];
sx q[2];
rz(-2.6878074) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.543964) q[1];
sx q[1];
rz(-1.9457327) q[1];
sx q[1];
rz(0.72245325) q[1];
rz(-0.28067067) q[3];
sx q[3];
rz(-1.9037399) q[3];
sx q[3];
rz(-1.5047634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7555776) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(1.8899274) q[2];
rz(-1.3517316) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(-0.98021093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1397454) q[0];
sx q[0];
rz(-2.503105) q[0];
sx q[0];
rz(-1.8852604) q[0];
rz(-0.095257692) q[1];
sx q[1];
rz(-2.2525411) q[1];
sx q[1];
rz(-0.5307861) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19671104) q[0];
sx q[0];
rz(-1.3232348) q[0];
sx q[0];
rz(0.16522466) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31802788) q[2];
sx q[2];
rz(-1.6075168) q[2];
sx q[2];
rz(1.6317692) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4640749) q[1];
sx q[1];
rz(-1.4428268) q[1];
sx q[1];
rz(-2.900338) q[1];
rz(1.4559326) q[3];
sx q[3];
rz(-2.0146166) q[3];
sx q[3];
rz(-0.97568363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.329616) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(-2.6007593) q[2];
rz(2.3234308) q[3];
sx q[3];
rz(-1.6988138) q[3];
sx q[3];
rz(-0.63280672) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5166017) q[0];
sx q[0];
rz(-1.3636959) q[0];
sx q[0];
rz(0.3987819) q[0];
rz(-1.7575691) q[1];
sx q[1];
rz(-2.3868491) q[1];
sx q[1];
rz(0.049364518) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2731664) q[0];
sx q[0];
rz(-0.12463649) q[0];
sx q[0];
rz(-1.0602555) q[0];
rz(0.10714679) q[2];
sx q[2];
rz(-1.8320036) q[2];
sx q[2];
rz(0.10814737) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80222622) q[1];
sx q[1];
rz(-0.58103937) q[1];
sx q[1];
rz(-1.6520348) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7743763) q[3];
sx q[3];
rz(-2.1927811) q[3];
sx q[3];
rz(-2.3443215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.12386879) q[2];
sx q[2];
rz(-1.8327291) q[2];
sx q[2];
rz(1.6309942) q[2];
rz(-0.92783582) q[3];
sx q[3];
rz(-1.3414914) q[3];
sx q[3];
rz(-1.235289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8265726) q[0];
sx q[0];
rz(-2.6980504) q[0];
sx q[0];
rz(2.0499688) q[0];
rz(-1.0646959) q[1];
sx q[1];
rz(-0.5263435) q[1];
sx q[1];
rz(1.975504) q[1];
rz(-2.0666368) q[2];
sx q[2];
rz(-1.1100162) q[2];
sx q[2];
rz(0.77592862) q[2];
rz(-2.8228685) q[3];
sx q[3];
rz(-1.0464077) q[3];
sx q[3];
rz(2.9628021) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
