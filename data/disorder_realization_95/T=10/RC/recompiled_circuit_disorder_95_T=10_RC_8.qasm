OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(0.89755091) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(2.0770567) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39660397) q[0];
sx q[0];
rz(-0.94192266) q[0];
sx q[0];
rz(3.0552342) q[0];
x q[1];
rz(-1.1041553) q[2];
sx q[2];
rz(-2.4541514) q[2];
sx q[2];
rz(-1.825037) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0734288) q[1];
sx q[1];
rz(-0.66232077) q[1];
sx q[1];
rz(-0.27889241) q[1];
x q[2];
rz(-0.1944794) q[3];
sx q[3];
rz(-0.76768657) q[3];
sx q[3];
rz(-0.0038113468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.66427461) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(-0.73195362) q[2];
rz(2.1814363) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17523781) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(0.91645855) q[0];
rz(2.6610999) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(-0.8786456) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71683305) q[0];
sx q[0];
rz(-0.75100198) q[0];
sx q[0];
rz(0.99320937) q[0];
x q[1];
rz(-2.5231045) q[2];
sx q[2];
rz(-1.3534708) q[2];
sx q[2];
rz(-0.49371142) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0658873) q[1];
sx q[1];
rz(-1.1527219) q[1];
sx q[1];
rz(1.4904406) q[1];
rz(-1.8804178) q[3];
sx q[3];
rz(-1.6358346) q[3];
sx q[3];
rz(-2.0464954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2800704) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(0.64727616) q[2];
rz(-0.17368008) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(-2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(1.1211959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66378731) q[0];
sx q[0];
rz(-0.87513798) q[0];
sx q[0];
rz(1.0870766) q[0];
rz(-0.063935117) q[2];
sx q[2];
rz(-1.1786412) q[2];
sx q[2];
rz(2.8901697) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0908302) q[1];
sx q[1];
rz(-1.275626) q[1];
sx q[1];
rz(-2.7961618) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6882012) q[3];
sx q[3];
rz(-2.061764) q[3];
sx q[3];
rz(-1.9672058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7401509) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(-2.1901954) q[2];
rz(-0.65008632) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6275416) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(0.78805584) q[0];
rz(-1.0568985) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(0.11638164) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1040092) q[0];
sx q[0];
rz(-1.5452256) q[0];
sx q[0];
rz(-2.5081162) q[0];
rz(2.0747244) q[2];
sx q[2];
rz(-1.908784) q[2];
sx q[2];
rz(-2.4525814) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0105646) q[1];
sx q[1];
rz(-2.4782135) q[1];
sx q[1];
rz(-2.4478711) q[1];
rz(-pi) q[2];
rz(-1.5192401) q[3];
sx q[3];
rz(-1.0815074) q[3];
sx q[3];
rz(-1.4435022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6115761) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(-2.8895565) q[2];
rz(2.7633372) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(-2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5761121) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(-0.53260032) q[0];
rz(-1.4415007) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(1.1486357) q[1];
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
rz(2.1158475) q[2];
sx q[2];
rz(-2.3902241) q[2];
sx q[2];
rz(-0.92912208) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7733113) q[1];
sx q[1];
rz(-1.9241663) q[1];
sx q[1];
rz(-1.7721304) q[1];
rz(-2.4162021) q[3];
sx q[3];
rz(-1.2687506) q[3];
sx q[3];
rz(-2.6186752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0107515) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(-2.1095236) q[2];
rz(2.4268835) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(-0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(2.7456039) q[0];
rz(-1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(2.9352303) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2287184) q[0];
sx q[0];
rz(-1.4836856) q[0];
sx q[0];
rz(0.74761439) q[0];
x q[1];
rz(0.88271898) q[2];
sx q[2];
rz(-2.1157017) q[2];
sx q[2];
rz(2.8628652) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3316162) q[1];
sx q[1];
rz(-0.14252256) q[1];
sx q[1];
rz(-0.82377164) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9456375) q[3];
sx q[3];
rz(-1.8926419) q[3];
sx q[3];
rz(2.4623507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6017194) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(-0.27077857) q[2];
rz(0.74832908) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(-1.07871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.8966184) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(2.5653429) q[0];
rz(0.17310625) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(1.9304088) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9528708) q[0];
sx q[0];
rz(-1.0591918) q[0];
sx q[0];
rz(0.94922025) q[0];
x q[1];
rz(-2.6518875) q[2];
sx q[2];
rz(-1.6929132) q[2];
sx q[2];
rz(2.7011938) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3019575) q[1];
sx q[1];
rz(-0.17050276) q[1];
sx q[1];
rz(-1.1595999) q[1];
rz(-pi) q[2];
rz(2.022642) q[3];
sx q[3];
rz(-2.4882836) q[3];
sx q[3];
rz(-2.2484231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8048191) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(0.024519196) q[2];
rz(0.71497861) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(-1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1361168) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(-1.0205644) q[0];
rz(2.9868946) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(1.9205836) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97604254) q[0];
sx q[0];
rz(-1.8414458) q[0];
sx q[0];
rz(-0.793215) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74322015) q[2];
sx q[2];
rz(-2.9033702) q[2];
sx q[2];
rz(-1.5526349) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.11215969) q[1];
sx q[1];
rz(-0.1754079) q[1];
sx q[1];
rz(2.4584241) q[1];
x q[2];
rz(1.9302759) q[3];
sx q[3];
rz(-0.7822789) q[3];
sx q[3];
rz(-0.6560916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(-0.70518804) q[2];
rz(0.60338902) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(1.5195297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.110638) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(-0.13701339) q[0];
rz(0.6048454) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(-3.0158214) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9916358) q[0];
sx q[0];
rz(-1.3949252) q[0];
sx q[0];
rz(-1.3791023) q[0];
rz(-pi) q[1];
rz(2.8081886) q[2];
sx q[2];
rz(-1.7575022) q[2];
sx q[2];
rz(-0.29078996) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6313435) q[1];
sx q[1];
rz(-1.8537087) q[1];
sx q[1];
rz(0.81088539) q[1];
rz(-2.5500507) q[3];
sx q[3];
rz(-0.35174832) q[3];
sx q[3];
rz(-1.2951617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.003309) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(-0.70927817) q[2];
rz(-2.5214031) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(2.9848849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9897292) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(-1.9412769) q[0];
rz(0.26750803) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(2.0013924) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4338019) q[0];
sx q[0];
rz(-1.8330488) q[0];
sx q[0];
rz(0.10284013) q[0];
x q[1];
rz(2.0613725) q[2];
sx q[2];
rz(-2.6274101) q[2];
sx q[2];
rz(1.0487923) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0298437) q[1];
sx q[1];
rz(-2.9330301) q[1];
sx q[1];
rz(-2.7010121) q[1];
rz(-pi) q[2];
rz(2.8285698) q[3];
sx q[3];
rz(-0.81837294) q[3];
sx q[3];
rz(2.5767874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(-0.89861384) q[2];
rz(0.0071772655) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.89467775) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(2.6208411) q[1];
sx q[1];
rz(-1.755935) q[1];
sx q[1];
rz(-1.4204949) q[1];
rz(-2.043914) q[2];
sx q[2];
rz(-2.0156779) q[2];
sx q[2];
rz(1.426522) q[2];
rz(2.5614212) q[3];
sx q[3];
rz(-2.4709354) q[3];
sx q[3];
rz(1.2448268) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
