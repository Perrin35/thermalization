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
rz(-2.3209651) q[0];
sx q[0];
rz(-2.4790915) q[0];
sx q[0];
rz(-2.878046) q[0];
rz(-2.0343556) q[1];
sx q[1];
rz(-1.2821481) q[1];
sx q[1];
rz(-0.67556226) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5731544) q[0];
sx q[0];
rz(-2.805925) q[0];
sx q[0];
rz(-1.5511309) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58568875) q[2];
sx q[2];
rz(-1.0694475) q[2];
sx q[2];
rz(0.96282178) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5225884) q[1];
sx q[1];
rz(-1.1824885) q[1];
sx q[1];
rz(-0.16242811) q[1];
x q[2];
rz(0.3278424) q[3];
sx q[3];
rz(-2.050542) q[3];
sx q[3];
rz(-2.858851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88966173) q[2];
sx q[2];
rz(-1.2222922) q[2];
sx q[2];
rz(-2.7737854) q[2];
rz(0.31428549) q[3];
sx q[3];
rz(-2.8511484) q[3];
sx q[3];
rz(2.9634326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80938584) q[0];
sx q[0];
rz(-0.089501373) q[0];
sx q[0];
rz(-1.7820763) q[0];
rz(1.8571732) q[1];
sx q[1];
rz(-1.1651243) q[1];
sx q[1];
rz(-3.0899835) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70003874) q[0];
sx q[0];
rz(-1.3897898) q[0];
sx q[0];
rz(-1.8405528) q[0];
rz(2.935595) q[2];
sx q[2];
rz(-1.9340746) q[2];
sx q[2];
rz(1.3000803) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9549841) q[1];
sx q[1];
rz(-2.767635) q[1];
sx q[1];
rz(-1.5692406) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4447024) q[3];
sx q[3];
rz(-2.5514388) q[3];
sx q[3];
rz(-1.7916288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.08903683) q[2];
sx q[2];
rz(-0.26036981) q[2];
sx q[2];
rz(0.54711771) q[2];
rz(1.268092) q[3];
sx q[3];
rz(-1.3601466) q[3];
sx q[3];
rz(-2.8030296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38076213) q[0];
sx q[0];
rz(-1.0365423) q[0];
sx q[0];
rz(1.8007675) q[0];
rz(-2.2876168) q[1];
sx q[1];
rz(-1.8546591) q[1];
sx q[1];
rz(0.30240789) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0722457) q[0];
sx q[0];
rz(-0.69127842) q[0];
sx q[0];
rz(-1.1290347) q[0];
rz(-pi) q[1];
rz(-1.339389) q[2];
sx q[2];
rz(-1.1303899) q[2];
sx q[2];
rz(-2.6803859) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.21061347) q[1];
sx q[1];
rz(-2.180778) q[1];
sx q[1];
rz(0.90175866) q[1];
x q[2];
rz(-2.0539927) q[3];
sx q[3];
rz(-2.5977166) q[3];
sx q[3];
rz(-2.9915006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2053908) q[2];
sx q[2];
rz(-1.8855636) q[2];
sx q[2];
rz(0.54226792) q[2];
rz(2.1451456) q[3];
sx q[3];
rz(-0.22765972) q[3];
sx q[3];
rz(-3.0136133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.19313136) q[0];
sx q[0];
rz(-1.6574991) q[0];
sx q[0];
rz(1.2633854) q[0];
rz(0.43426934) q[1];
sx q[1];
rz(-0.77249211) q[1];
sx q[1];
rz(-1.6901406) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5562486) q[0];
sx q[0];
rz(-2.5894269) q[0];
sx q[0];
rz(-0.57741965) q[0];
x q[1];
rz(2.5927438) q[2];
sx q[2];
rz(-1.9611231) q[2];
sx q[2];
rz(3.0496719) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2512379) q[1];
sx q[1];
rz(-1.4571937) q[1];
sx q[1];
rz(3.0556728) q[1];
rz(-pi) q[2];
rz(0.85335387) q[3];
sx q[3];
rz(-2.382016) q[3];
sx q[3];
rz(2.4977998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0749391) q[2];
sx q[2];
rz(-2.1982919) q[2];
sx q[2];
rz(-0.88610506) q[2];
rz(0.0063627176) q[3];
sx q[3];
rz(-1.3402091) q[3];
sx q[3];
rz(1.1191012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5804382) q[0];
sx q[0];
rz(-2.6851324) q[0];
sx q[0];
rz(1.891267) q[0];
rz(2.9087032) q[1];
sx q[1];
rz(-0.75585514) q[1];
sx q[1];
rz(0.09045352) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4912495) q[0];
sx q[0];
rz(-1.8199662) q[0];
sx q[0];
rz(0.12604524) q[0];
rz(1.8402507) q[2];
sx q[2];
rz(-1.1134992) q[2];
sx q[2];
rz(-2.2895165) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.0098086987) q[1];
sx q[1];
rz(-2.4854641) q[1];
sx q[1];
rz(-3.0387278) q[1];
rz(1.2430322) q[3];
sx q[3];
rz(-1.48945) q[3];
sx q[3];
rz(1.7880754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.49454409) q[2];
sx q[2];
rz(-2.9477951) q[2];
sx q[2];
rz(-1.9336644) q[2];
rz(-0.9211933) q[3];
sx q[3];
rz(-0.84417206) q[3];
sx q[3];
rz(-2.6827961) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0334466) q[0];
sx q[0];
rz(-1.7195846) q[0];
sx q[0];
rz(-2.6500927) q[0];
rz(-2.4132701) q[1];
sx q[1];
rz(-2.0305384) q[1];
sx q[1];
rz(-2.947015) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96724961) q[0];
sx q[0];
rz(-1.3493269) q[0];
sx q[0];
rz(-2.1465149) q[0];
rz(-0.80912239) q[2];
sx q[2];
rz(-2.1702907) q[2];
sx q[2];
rz(0.15669151) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66524551) q[1];
sx q[1];
rz(-2.0788296) q[1];
sx q[1];
rz(0.82583896) q[1];
x q[2];
rz(-1.803502) q[3];
sx q[3];
rz(-0.92706087) q[3];
sx q[3];
rz(-0.88730592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7877833) q[2];
sx q[2];
rz(-1.0156735) q[2];
sx q[2];
rz(-0.0018399012) q[2];
rz(1.6483824) q[3];
sx q[3];
rz(-2.4108678) q[3];
sx q[3];
rz(-0.28436896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2581185) q[0];
sx q[0];
rz(-0.26171568) q[0];
sx q[0];
rz(-2.0765685) q[0];
rz(-0.26271543) q[1];
sx q[1];
rz(-2.5698667) q[1];
sx q[1];
rz(-2.6782742) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7119448) q[0];
sx q[0];
rz(-1.4139626) q[0];
sx q[0];
rz(0.98008396) q[0];
rz(0.72603308) q[2];
sx q[2];
rz(-1.9872857) q[2];
sx q[2];
rz(-1.287578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1541248) q[1];
sx q[1];
rz(-0.5041669) q[1];
sx q[1];
rz(-2.8567863) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3635395) q[3];
sx q[3];
rz(-0.1175783) q[3];
sx q[3];
rz(-2.1685947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.00039438417) q[2];
sx q[2];
rz(-1.980915) q[2];
sx q[2];
rz(-2.311643) q[2];
rz(-2.0937894) q[3];
sx q[3];
rz(-2.6959097) q[3];
sx q[3];
rz(-0.1786264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038427) q[0];
sx q[0];
rz(-2.1512845) q[0];
sx q[0];
rz(0.75482279) q[0];
rz(0.38914514) q[1];
sx q[1];
rz(-1.4578338) q[1];
sx q[1];
rz(0.39774242) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75466418) q[0];
sx q[0];
rz(-2.7542186) q[0];
sx q[0];
rz(2.3088916) q[0];
rz(-pi) q[1];
rz(-2.7806578) q[2];
sx q[2];
rz(-0.58141469) q[2];
sx q[2];
rz(-1.5370528) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.81521846) q[1];
sx q[1];
rz(-0.37719101) q[1];
sx q[1];
rz(0.81881028) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3944309) q[3];
sx q[3];
rz(-2.0748169) q[3];
sx q[3];
rz(2.6256743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8336746) q[2];
sx q[2];
rz(-0.29611823) q[2];
sx q[2];
rz(2.8274242) q[2];
rz(-1.2416154) q[3];
sx q[3];
rz(-1.5528468) q[3];
sx q[3];
rz(-1.2392905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.817713) q[0];
sx q[0];
rz(-2.4535024) q[0];
sx q[0];
rz(-2.896198) q[0];
rz(1.2303526) q[1];
sx q[1];
rz(-3.0407258) q[1];
sx q[1];
rz(2.6255677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2986261) q[0];
sx q[0];
rz(-0.43382513) q[0];
sx q[0];
rz(2.6515657) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74036171) q[2];
sx q[2];
rz(-0.90351331) q[2];
sx q[2];
rz(1.1297392) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.302205) q[1];
sx q[1];
rz(-2.6714462) q[1];
sx q[1];
rz(-0.66580982) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13652674) q[3];
sx q[3];
rz(-1.5093922) q[3];
sx q[3];
rz(1.9207973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5838991) q[2];
sx q[2];
rz(-1.9198753) q[2];
sx q[2];
rz(-1.0830797) q[2];
rz(0.22100581) q[3];
sx q[3];
rz(-2.998304) q[3];
sx q[3];
rz(-2.3451282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1204656) q[0];
sx q[0];
rz(-3.0607304) q[0];
sx q[0];
rz(3.0057111) q[0];
rz(1.4295626) q[1];
sx q[1];
rz(-2.0202961) q[1];
sx q[1];
rz(-2.8543465) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3843975) q[0];
sx q[0];
rz(-1.383019) q[0];
sx q[0];
rz(1.8191992) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9075228) q[2];
sx q[2];
rz(-2.8904466) q[2];
sx q[2];
rz(2.2876571) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4212227) q[1];
sx q[1];
rz(-2.1541641) q[1];
sx q[1];
rz(-1.3066584) q[1];
rz(-pi) q[2];
rz(1.997825) q[3];
sx q[3];
rz(-2.5134519) q[3];
sx q[3];
rz(2.028156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.91934943) q[2];
sx q[2];
rz(-1.9025849) q[2];
sx q[2];
rz(0.83402056) q[2];
rz(-0.30759865) q[3];
sx q[3];
rz(-0.37720507) q[3];
sx q[3];
rz(-2.8789177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.96497) q[0];
sx q[0];
rz(-1.8360092) q[0];
sx q[0];
rz(-1.0259884) q[0];
rz(2.7550244) q[1];
sx q[1];
rz(-1.3810806) q[1];
sx q[1];
rz(2.5234533) q[1];
rz(2.8564524) q[2];
sx q[2];
rz(-0.44024537) q[2];
sx q[2];
rz(2.5633752) q[2];
rz(-1.3597506) q[3];
sx q[3];
rz(-1.4283709) q[3];
sx q[3];
rz(-1.4621468) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
