OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.14804949) q[0];
sx q[0];
rz(-0.43819675) q[0];
sx q[0];
rz(0.73102695) q[0];
rz(-2.464715) q[1];
sx q[1];
rz(2.4426065) q[1];
sx q[1];
rz(12.01233) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95931388) q[0];
sx q[0];
rz(-2.4327299) q[0];
sx q[0];
rz(2.3793263) q[0];
rz(-pi) q[1];
rz(-1.117716) q[2];
sx q[2];
rz(-0.52342285) q[2];
sx q[2];
rz(-2.7414309) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3202127) q[1];
sx q[1];
rz(-2.2578635) q[1];
sx q[1];
rz(-3.0979731) q[1];
rz(-0.88907928) q[3];
sx q[3];
rz(-1.9535716) q[3];
sx q[3];
rz(1.0600404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76250166) q[2];
sx q[2];
rz(-1.5753626) q[2];
sx q[2];
rz(2.9762034) q[2];
rz(-2.9108544) q[3];
sx q[3];
rz(-0.24300353) q[3];
sx q[3];
rz(-0.86133426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6623401) q[0];
sx q[0];
rz(-2.8133744) q[0];
sx q[0];
rz(1.7829371) q[0];
rz(-1.5232297) q[1];
sx q[1];
rz(-2.3871469) q[1];
sx q[1];
rz(0.84017909) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62094394) q[0];
sx q[0];
rz(-2.0387702) q[0];
sx q[0];
rz(-2.8415989) q[0];
rz(-pi) q[1];
rz(0.86960881) q[2];
sx q[2];
rz(-1.879862) q[2];
sx q[2];
rz(0.78534195) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3749373) q[1];
sx q[1];
rz(-1.8647653) q[1];
sx q[1];
rz(-0.73237441) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9884913) q[3];
sx q[3];
rz(-0.77147064) q[3];
sx q[3];
rz(0.16530802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.75300616) q[2];
sx q[2];
rz(-1.9084088) q[2];
sx q[2];
rz(-2.2071655) q[2];
rz(1.2855444) q[3];
sx q[3];
rz(-0.779874) q[3];
sx q[3];
rz(0.88750315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11338209) q[0];
sx q[0];
rz(-2.1206355) q[0];
sx q[0];
rz(-1.5895948) q[0];
rz(-1.3681083) q[1];
sx q[1];
rz(-1.869447) q[1];
sx q[1];
rz(1.1205463) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5277953) q[0];
sx q[0];
rz(-0.58301914) q[0];
sx q[0];
rz(0.9436508) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0650313) q[2];
sx q[2];
rz(-2.6496844) q[2];
sx q[2];
rz(3.0658403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6419598) q[1];
sx q[1];
rz(-1.4781535) q[1];
sx q[1];
rz(-0.53643061) q[1];
rz(-1.9777704) q[3];
sx q[3];
rz(-1.455179) q[3];
sx q[3];
rz(-1.8979929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0997194) q[2];
sx q[2];
rz(-2.9283044) q[2];
sx q[2];
rz(-2.8495157) q[2];
rz(1.9000351) q[3];
sx q[3];
rz(-1.440666) q[3];
sx q[3];
rz(-0.45274538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1555772) q[0];
sx q[0];
rz(-0.41322511) q[0];
sx q[0];
rz(-0.43169942) q[0];
rz(-0.15402928) q[1];
sx q[1];
rz(-1.5715716) q[1];
sx q[1];
rz(1.0607176) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5397388) q[0];
sx q[0];
rz(-2.665747) q[0];
sx q[0];
rz(-0.81625508) q[0];
x q[1];
rz(1.3260822) q[2];
sx q[2];
rz(-1.7438466) q[2];
sx q[2];
rz(1.7987378) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1338443) q[1];
sx q[1];
rz(-2.7415726) q[1];
sx q[1];
rz(0.96537867) q[1];
rz(-pi) q[2];
rz(-0.19017724) q[3];
sx q[3];
rz(-1.7253644) q[3];
sx q[3];
rz(3.0319253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3920307) q[2];
sx q[2];
rz(-1.531484) q[2];
sx q[2];
rz(2.5677666) q[2];
rz(3.0841893) q[3];
sx q[3];
rz(-1.4797689) q[3];
sx q[3];
rz(0.81006947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5350128) q[0];
sx q[0];
rz(-2.8856394) q[0];
sx q[0];
rz(2.8888597) q[0];
rz(-0.17929721) q[1];
sx q[1];
rz(-1.7959692) q[1];
sx q[1];
rz(-1.7326694) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64302639) q[0];
sx q[0];
rz(-1.4966949) q[0];
sx q[0];
rz(-1.3331494) q[0];
x q[1];
rz(-1.9581355) q[2];
sx q[2];
rz(-1.172386) q[2];
sx q[2];
rz(-0.18926316) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.682413) q[1];
sx q[1];
rz(-2.3485567) q[1];
sx q[1];
rz(-2.9694818) q[1];
x q[2];
rz(0.5770251) q[3];
sx q[3];
rz(-0.60363704) q[3];
sx q[3];
rz(1.9848422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0130284) q[2];
sx q[2];
rz(-1.7594254) q[2];
sx q[2];
rz(-1.9405091) q[2];
rz(2.3342093) q[3];
sx q[3];
rz(-1.3599334) q[3];
sx q[3];
rz(2.0090296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(0.078868911) q[0];
sx q[0];
rz(-1.6428592) q[0];
sx q[0];
rz(-2.3003182) q[0];
rz(1.2048644) q[1];
sx q[1];
rz(-1.3179904) q[1];
sx q[1];
rz(-2.1119609) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1618274) q[0];
sx q[0];
rz(-1.0603535) q[0];
sx q[0];
rz(-0.011988601) q[0];
rz(-pi) q[1];
rz(2.6961521) q[2];
sx q[2];
rz(-2.5800319) q[2];
sx q[2];
rz(0.5024006) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8334044) q[1];
sx q[1];
rz(-2.9627315) q[1];
sx q[1];
rz(-0.043828242) q[1];
x q[2];
rz(-2.5978686) q[3];
sx q[3];
rz(-2.7910376) q[3];
sx q[3];
rz(-1.9026827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.131375) q[2];
sx q[2];
rz(-2.3848332) q[2];
sx q[2];
rz(-0.49794623) q[2];
rz(-2.61854) q[3];
sx q[3];
rz(-1.4191351) q[3];
sx q[3];
rz(-3.0150692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2987357) q[0];
sx q[0];
rz(-3.0099478) q[0];
sx q[0];
rz(-0.81197062) q[0];
rz(-2.2506591) q[1];
sx q[1];
rz(-1.1633326) q[1];
sx q[1];
rz(-0.99937159) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1645714) q[0];
sx q[0];
rz(-2.0334822) q[0];
sx q[0];
rz(2.0624731) q[0];
rz(0.096316263) q[2];
sx q[2];
rz(-1.4681946) q[2];
sx q[2];
rz(1.4136537) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68518406) q[1];
sx q[1];
rz(-1.0899836) q[1];
sx q[1];
rz(0.64688869) q[1];
rz(-2.8511397) q[3];
sx q[3];
rz(-2.2155665) q[3];
sx q[3];
rz(-1.6381263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.823395) q[2];
sx q[2];
rz(-2.2237033) q[2];
sx q[2];
rz(-1.8434175) q[2];
rz(-0.038481742) q[3];
sx q[3];
rz(-2.0352071) q[3];
sx q[3];
rz(-1.6160256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0865974) q[0];
sx q[0];
rz(-2.0107858) q[0];
sx q[0];
rz(0.31563345) q[0];
rz(1.9075958) q[1];
sx q[1];
rz(-2.4256568) q[1];
sx q[1];
rz(-1.2589781) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0961279) q[0];
sx q[0];
rz(-0.79497951) q[0];
sx q[0];
rz(0.40277512) q[0];
x q[1];
rz(-1.3975436) q[2];
sx q[2];
rz(-2.3286741) q[2];
sx q[2];
rz(-0.11101857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1728744) q[1];
sx q[1];
rz(-0.58191381) q[1];
sx q[1];
rz(-0.7919148) q[1];
x q[2];
rz(-0.057825967) q[3];
sx q[3];
rz(-3.0396121) q[3];
sx q[3];
rz(-1.1241871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.2235609) q[2];
sx q[2];
rz(-0.74350244) q[2];
sx q[2];
rz(-3.1136759) q[2];
rz(0.74808407) q[3];
sx q[3];
rz(-1.7223765) q[3];
sx q[3];
rz(2.9885651) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7276723) q[0];
sx q[0];
rz(-2.1509009) q[0];
sx q[0];
rz(2.682611) q[0];
rz(-2.0052295) q[1];
sx q[1];
rz(-1.4308948) q[1];
sx q[1];
rz(-2.9232025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1735575) q[0];
sx q[0];
rz(-1.5941248) q[0];
sx q[0];
rz(-1.5874046) q[0];
rz(-pi) q[1];
rz(1.2480466) q[2];
sx q[2];
rz(-0.40229978) q[2];
sx q[2];
rz(2.2108271) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23109197) q[1];
sx q[1];
rz(-1.0808966) q[1];
sx q[1];
rz(-1.7375768) q[1];
x q[2];
rz(0.13111643) q[3];
sx q[3];
rz(-1.9558441) q[3];
sx q[3];
rz(-2.7484886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8454664) q[2];
sx q[2];
rz(-2.8870236) q[2];
sx q[2];
rz(2.0938342) q[2];
rz(-0.74538499) q[3];
sx q[3];
rz(-1.5630009) q[3];
sx q[3];
rz(-1.6489702) q[3];
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
rz(0.02136852) q[0];
sx q[0];
rz(-0.57149082) q[0];
sx q[0];
rz(-2.7987203) q[0];
rz(-2.951237) q[1];
sx q[1];
rz(-1.0089259) q[1];
sx q[1];
rz(2.6284133) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8489055) q[0];
sx q[0];
rz(-0.93528895) q[0];
sx q[0];
rz(0.19750316) q[0];
rz(-pi) q[1];
rz(-2.1634899) q[2];
sx q[2];
rz(-0.58916559) q[2];
sx q[2];
rz(-1.9798673) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3575705) q[1];
sx q[1];
rz(-2.5139132) q[1];
sx q[1];
rz(0.72369544) q[1];
rz(-pi) q[2];
rz(1.1499728) q[3];
sx q[3];
rz(-1.7003891) q[3];
sx q[3];
rz(1.1269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.43789431) q[2];
sx q[2];
rz(-0.63592211) q[2];
sx q[2];
rz(-0.61100125) q[2];
rz(0.25887394) q[3];
sx q[3];
rz(-1.9301819) q[3];
sx q[3];
rz(1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7543058) q[0];
sx q[0];
rz(-1.5504693) q[0];
sx q[0];
rz(-1.5691527) q[0];
rz(0.98987956) q[1];
sx q[1];
rz(-1.2073333) q[1];
sx q[1];
rz(-2.5055199) q[1];
rz(-1.5237332) q[2];
sx q[2];
rz(-0.18885352) q[2];
sx q[2];
rz(0.39299008) q[2];
rz(-2.4425735) q[3];
sx q[3];
rz(-1.6215646) q[3];
sx q[3];
rz(-0.42536964) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
