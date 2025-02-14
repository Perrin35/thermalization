OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6736421) q[0];
sx q[0];
rz(-1.8704432) q[0];
sx q[0];
rz(-0.76572642) q[0];
rz(-0.39628059) q[1];
sx q[1];
rz(-0.09274617) q[1];
sx q[1];
rz(0.5527817) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46229773) q[0];
sx q[0];
rz(-0.57749417) q[0];
sx q[0];
rz(-0.80120835) q[0];
x q[1];
rz(-0.67019318) q[2];
sx q[2];
rz(-0.92447058) q[2];
sx q[2];
rz(-2.0984416) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1172792) q[1];
sx q[1];
rz(-1.8872042) q[1];
sx q[1];
rz(-2.1561619) q[1];
rz(-pi) q[2];
rz(2.1215274) q[3];
sx q[3];
rz(-2.4006776) q[3];
sx q[3];
rz(-1.7307841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7736241) q[2];
sx q[2];
rz(-1.5017193) q[2];
sx q[2];
rz(3.0527414) q[2];
rz(2.7446274) q[3];
sx q[3];
rz(-1.0395972) q[3];
sx q[3];
rz(1.4575492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6615768) q[0];
sx q[0];
rz(-0.50742298) q[0];
sx q[0];
rz(-2.2870824) q[0];
rz(-0.62141934) q[1];
sx q[1];
rz(-0.3365376) q[1];
sx q[1];
rz(1.052676) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7490279) q[0];
sx q[0];
rz(-1.6727007) q[0];
sx q[0];
rz(3.0046731) q[0];
x q[1];
rz(-0.41469736) q[2];
sx q[2];
rz(-1.1578348) q[2];
sx q[2];
rz(0.1916445) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.32089511) q[1];
sx q[1];
rz(-1.7716494) q[1];
sx q[1];
rz(2.3066375) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0016453513) q[3];
sx q[3];
rz(-0.96685322) q[3];
sx q[3];
rz(1.920948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4738327) q[2];
sx q[2];
rz(-2.2809873) q[2];
sx q[2];
rz(-2.1916126) q[2];
rz(-1.0085603) q[3];
sx q[3];
rz(-0.63298321) q[3];
sx q[3];
rz(0.27957255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94473332) q[0];
sx q[0];
rz(-1.9254528) q[0];
sx q[0];
rz(1.3470294) q[0];
rz(1.0579717) q[1];
sx q[1];
rz(-2.8949013) q[1];
sx q[1];
rz(0.11014858) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5025217) q[0];
sx q[0];
rz(-2.5701231) q[0];
sx q[0];
rz(2.1854464) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0827237) q[2];
sx q[2];
rz(-2.9721568) q[2];
sx q[2];
rz(-1.2537341) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0855838) q[1];
sx q[1];
rz(-0.97404503) q[1];
sx q[1];
rz(1.4419773) q[1];
rz(-pi) q[2];
rz(-2.0581117) q[3];
sx q[3];
rz(-1.1705453) q[3];
sx q[3];
rz(2.6292173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9614253) q[2];
sx q[2];
rz(-1.578293) q[2];
sx q[2];
rz(-0.041672826) q[2];
rz(-0.20798072) q[3];
sx q[3];
rz(-2.9198923) q[3];
sx q[3];
rz(2.8878133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5623077) q[0];
sx q[0];
rz(-1.8203745) q[0];
sx q[0];
rz(-2.1606309) q[0];
rz(0.43308577) q[1];
sx q[1];
rz(-1.4739477) q[1];
sx q[1];
rz(1.513419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.604852) q[0];
sx q[0];
rz(-1.3910595) q[0];
sx q[0];
rz(2.7374703) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89940091) q[2];
sx q[2];
rz(-1.5601306) q[2];
sx q[2];
rz(0.15029066) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95483741) q[1];
sx q[1];
rz(-1.673497) q[1];
sx q[1];
rz(1.0726733) q[1];
rz(2.8582005) q[3];
sx q[3];
rz(-0.95855721) q[3];
sx q[3];
rz(-1.3536842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5341586) q[2];
sx q[2];
rz(-2.6223493) q[2];
sx q[2];
rz(-2.329211) q[2];
rz(1.8937998) q[3];
sx q[3];
rz(-1.2500074) q[3];
sx q[3];
rz(-1.8947424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8060551) q[0];
sx q[0];
rz(-2.1195109) q[0];
sx q[0];
rz(-1.6988276) q[0];
rz(0.45817786) q[1];
sx q[1];
rz(-1.6280326) q[1];
sx q[1];
rz(-2.3853669) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0201705) q[0];
sx q[0];
rz(-0.64682942) q[0];
sx q[0];
rz(1.3117164) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5014072) q[2];
sx q[2];
rz(-2.019503) q[2];
sx q[2];
rz(-1.2011091) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9994558) q[1];
sx q[1];
rz(-0.161006) q[1];
sx q[1];
rz(-2.5672037) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2268158) q[3];
sx q[3];
rz(-1.5369253) q[3];
sx q[3];
rz(0.86465166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.68236399) q[2];
sx q[2];
rz(-1.3882779) q[2];
sx q[2];
rz(-0.96561042) q[2];
rz(2.5718578) q[3];
sx q[3];
rz(-2.0163048) q[3];
sx q[3];
rz(-1.6857356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4023034) q[0];
sx q[0];
rz(-1.1981413) q[0];
sx q[0];
rz(1.385561) q[0];
rz(0.7631453) q[1];
sx q[1];
rz(-1.6940073) q[1];
sx q[1];
rz(0.10173434) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0245314) q[0];
sx q[0];
rz(-2.5145338) q[0];
sx q[0];
rz(-1.7282906) q[0];
rz(-3.0252803) q[2];
sx q[2];
rz(-0.38517932) q[2];
sx q[2];
rz(-2.0907837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.26120034) q[1];
sx q[1];
rz(-1.8414547) q[1];
sx q[1];
rz(-3.0780869) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2564628) q[3];
sx q[3];
rz(-0.99381522) q[3];
sx q[3];
rz(-2.5654405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.10852854) q[2];
sx q[2];
rz(-1.6988924) q[2];
sx q[2];
rz(2.4326883) q[2];
rz(2.3809643) q[3];
sx q[3];
rz(-1.1516738) q[3];
sx q[3];
rz(2.2061548) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0660504) q[0];
sx q[0];
rz(-1.1533371) q[0];
sx q[0];
rz(0.72108889) q[0];
rz(-0.9779633) q[1];
sx q[1];
rz(-2.5826192) q[1];
sx q[1];
rz(-1.5616547) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2943731) q[0];
sx q[0];
rz(-0.10589639) q[0];
sx q[0];
rz(-0.76283331) q[0];
rz(2.1660837) q[2];
sx q[2];
rz(-1.3689318) q[2];
sx q[2];
rz(1.5306115) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9813919) q[1];
sx q[1];
rz(-2.1213648) q[1];
sx q[1];
rz(2.0345313) q[1];
rz(2.7851289) q[3];
sx q[3];
rz(-0.93178669) q[3];
sx q[3];
rz(1.654341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.28479031) q[2];
sx q[2];
rz(-2.6140116) q[2];
sx q[2];
rz(-1.6214726) q[2];
rz(-2.704845) q[3];
sx q[3];
rz(-2.2468061) q[3];
sx q[3];
rz(0.53957087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0590416) q[0];
sx q[0];
rz(-2.3024547) q[0];
sx q[0];
rz(0.82373291) q[0];
rz(-2.0027022) q[1];
sx q[1];
rz(-1.7887807) q[1];
sx q[1];
rz(1.1762071) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19768077) q[0];
sx q[0];
rz(-1.0689387) q[0];
sx q[0];
rz(-3.0905064) q[0];
x q[1];
rz(-2.5912639) q[2];
sx q[2];
rz(-0.68745733) q[2];
sx q[2];
rz(-0.94151173) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2580793) q[1];
sx q[1];
rz(-0.23830676) q[1];
sx q[1];
rz(0.041447354) q[1];
rz(-pi) q[2];
rz(-0.59112143) q[3];
sx q[3];
rz(-1.6038461) q[3];
sx q[3];
rz(-0.37632468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1463683) q[2];
sx q[2];
rz(-2.5547042) q[2];
sx q[2];
rz(-2.5541019) q[2];
rz(2.7557709) q[3];
sx q[3];
rz(-0.84857517) q[3];
sx q[3];
rz(-2.2390656) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1189482) q[0];
sx q[0];
rz(-2.1653439) q[0];
sx q[0];
rz(-2.5318085) q[0];
rz(-1.3439641) q[1];
sx q[1];
rz(-1.1577497) q[1];
sx q[1];
rz(-2.9885898) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1399779) q[0];
sx q[0];
rz(-1.5511594) q[0];
sx q[0];
rz(-2.6756712) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.010092334) q[2];
sx q[2];
rz(-2.3463425) q[2];
sx q[2];
rz(-2.681596) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9003332) q[1];
sx q[1];
rz(-1.5894842) q[1];
sx q[1];
rz(-0.92872031) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2629635) q[3];
sx q[3];
rz(-1.9464916) q[3];
sx q[3];
rz(3.0082448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1028221) q[2];
sx q[2];
rz(-1.5800579) q[2];
sx q[2];
rz(-2.5193396) q[2];
rz(-1.2004987) q[3];
sx q[3];
rz(-1.7222722) q[3];
sx q[3];
rz(0.57435575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95933652) q[0];
sx q[0];
rz(-0.065956235) q[0];
sx q[0];
rz(0.35520735) q[0];
rz(1.1025053) q[1];
sx q[1];
rz(-3.1246154) q[1];
sx q[1];
rz(-0.65974081) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44910494) q[0];
sx q[0];
rz(-2.4665255) q[0];
sx q[0];
rz(0.27691845) q[0];
rz(-pi) q[1];
rz(1.8546472) q[2];
sx q[2];
rz(-0.42869332) q[2];
sx q[2];
rz(-2.5108166) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7463508) q[1];
sx q[1];
rz(-1.244427) q[1];
sx q[1];
rz(-2.8119948) q[1];
rz(-pi) q[2];
rz(0.89307488) q[3];
sx q[3];
rz(-1.5297946) q[3];
sx q[3];
rz(-1.6779016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.79992646) q[2];
sx q[2];
rz(-3.0604) q[2];
sx q[2];
rz(0.1989092) q[2];
rz(2.343446) q[3];
sx q[3];
rz(-1.1856368) q[3];
sx q[3];
rz(1.2822436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.2220919) q[0];
sx q[0];
rz(-1.6035447) q[0];
sx q[0];
rz(-1.0731687) q[0];
rz(-2.3566698) q[1];
sx q[1];
rz(-0.85826086) q[1];
sx q[1];
rz(0.8716743) q[1];
rz(-1.6941407) q[2];
sx q[2];
rz(-1.5236749) q[2];
sx q[2];
rz(2.9070367) q[2];
rz(2.4724805) q[3];
sx q[3];
rz(-1.3539066) q[3];
sx q[3];
rz(2.2216968) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
