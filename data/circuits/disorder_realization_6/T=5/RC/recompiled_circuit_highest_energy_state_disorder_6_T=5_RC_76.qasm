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
rz(1.037984) q[0];
sx q[0];
rz(4.8259566) q[0];
sx q[0];
rz(11.2136) q[0];
rz(1.1846722) q[1];
sx q[1];
rz(-0.67263043) q[1];
sx q[1];
rz(1.5157359) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88600776) q[0];
sx q[0];
rz(-2.8681462) q[0];
sx q[0];
rz(-0.81583623) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0667384) q[2];
sx q[2];
rz(-2.3124394) q[2];
sx q[2];
rz(3.0665891) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.63187088) q[1];
sx q[1];
rz(-2.0364822) q[1];
sx q[1];
rz(1.9655955) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42077371) q[3];
sx q[3];
rz(-1.3075365) q[3];
sx q[3];
rz(-0.0076310633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1379913) q[2];
sx q[2];
rz(-3.0333952) q[2];
sx q[2];
rz(0.18599621) q[2];
rz(3.1229535) q[3];
sx q[3];
rz(-1.2942856) q[3];
sx q[3];
rz(2.0712461) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5587392) q[0];
sx q[0];
rz(-2.5753729) q[0];
sx q[0];
rz(-2.8232316) q[0];
rz(-0.63703498) q[1];
sx q[1];
rz(-2.36167) q[1];
sx q[1];
rz(0.46924082) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1701654) q[0];
sx q[0];
rz(-1.3999456) q[0];
sx q[0];
rz(2.9396482) q[0];
rz(-pi) q[1];
rz(3.0576433) q[2];
sx q[2];
rz(-1.1962039) q[2];
sx q[2];
rz(-2.2427246) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0035867) q[1];
sx q[1];
rz(-2.7452861) q[1];
sx q[1];
rz(0.62162865) q[1];
x q[2];
rz(-1.5777228) q[3];
sx q[3];
rz(-0.68994656) q[3];
sx q[3];
rz(-1.6405218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8239173) q[2];
sx q[2];
rz(-0.21024148) q[2];
sx q[2];
rz(-0.69327411) q[2];
rz(1.3200101) q[3];
sx q[3];
rz(-1.8157248) q[3];
sx q[3];
rz(2.0461953) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0005118) q[0];
sx q[0];
rz(-0.63758272) q[0];
sx q[0];
rz(2.6620423) q[0];
rz(2.6374822) q[1];
sx q[1];
rz(-1.9934374) q[1];
sx q[1];
rz(3.0832916) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2445912) q[0];
sx q[0];
rz(-2.069733) q[0];
sx q[0];
rz(1.2759491) q[0];
x q[1];
rz(2.53251) q[2];
sx q[2];
rz(-2.3533965) q[2];
sx q[2];
rz(2.4969375) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7229798) q[1];
sx q[1];
rz(-1.449391) q[1];
sx q[1];
rz(2.1239807) q[1];
x q[2];
rz(0.4543484) q[3];
sx q[3];
rz(-0.55255167) q[3];
sx q[3];
rz(-2.1228028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38640675) q[2];
sx q[2];
rz(-1.1815973) q[2];
sx q[2];
rz(3.1043261) q[2];
rz(1.6218328) q[3];
sx q[3];
rz(-2.683679) q[3];
sx q[3];
rz(-0.38145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.195049) q[0];
sx q[0];
rz(-0.46849546) q[0];
sx q[0];
rz(0.066548912) q[0];
rz(1.6171803) q[1];
sx q[1];
rz(-1.2434554) q[1];
sx q[1];
rz(0.7349416) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3108705) q[0];
sx q[0];
rz(-1.8025711) q[0];
sx q[0];
rz(-0.59344296) q[0];
rz(2.8822172) q[2];
sx q[2];
rz(-0.7470567) q[2];
sx q[2];
rz(0.41058985) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1030324) q[1];
sx q[1];
rz(-1.8502697) q[1];
sx q[1];
rz(-2.7725262) q[1];
x q[2];
rz(-2.9959277) q[3];
sx q[3];
rz(-0.85547912) q[3];
sx q[3];
rz(3.0112034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.28896004) q[2];
sx q[2];
rz(-1.5639049) q[2];
sx q[2];
rz(0.10870474) q[2];
rz(0.92681256) q[3];
sx q[3];
rz(-0.97065297) q[3];
sx q[3];
rz(-1.5638117) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0842593) q[0];
sx q[0];
rz(-2.5045392) q[0];
sx q[0];
rz(2.0507226) q[0];
rz(2.5690761) q[1];
sx q[1];
rz(-1.1895836) q[1];
sx q[1];
rz(-0.82430878) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0181247) q[0];
sx q[0];
rz(-1.195692) q[0];
sx q[0];
rz(-0.25183046) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9277497) q[2];
sx q[2];
rz(-1.130419) q[2];
sx q[2];
rz(-1.3470896) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.2356154) q[1];
sx q[1];
rz(-2.4433377) q[1];
sx q[1];
rz(-0.92963107) q[1];
rz(-pi) q[2];
rz(2.7359928) q[3];
sx q[3];
rz(-1.1685017) q[3];
sx q[3];
rz(1.9674439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.480964) q[2];
sx q[2];
rz(-0.97233665) q[2];
sx q[2];
rz(2.8503897) q[2];
rz(2.5045942) q[3];
sx q[3];
rz(-1.6773418) q[3];
sx q[3];
rz(1.252482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3785192) q[0];
sx q[0];
rz(-2.851649) q[0];
sx q[0];
rz(3.0105403) q[0];
rz(-1.9746926) q[1];
sx q[1];
rz(-0.98375541) q[1];
sx q[1];
rz(0.022620591) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2065434) q[0];
sx q[0];
rz(-0.66842043) q[0];
sx q[0];
rz(2.5517625) q[0];
x q[1];
rz(-0.51949595) q[2];
sx q[2];
rz(-0.24458376) q[2];
sx q[2];
rz(-1.0102538) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.80161051) q[1];
sx q[1];
rz(-1.2622854) q[1];
sx q[1];
rz(2.1478081) q[1];
rz(-pi) q[2];
rz(-0.87225391) q[3];
sx q[3];
rz(-0.24980751) q[3];
sx q[3];
rz(-1.8810617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.268198) q[2];
sx q[2];
rz(-2.4384273) q[2];
sx q[2];
rz(2.0568636) q[2];
rz(-1.1449413) q[3];
sx q[3];
rz(-1.8549253) q[3];
sx q[3];
rz(2.1196608) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5893843) q[0];
sx q[0];
rz(-1.6227868) q[0];
sx q[0];
rz(2.6680706) q[0];
rz(-1.2208968) q[1];
sx q[1];
rz(-0.41243204) q[1];
sx q[1];
rz(-1.268505) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.761577) q[0];
sx q[0];
rz(-0.76671874) q[0];
sx q[0];
rz(-0.79464998) q[0];
rz(-pi) q[1];
rz(2.2008197) q[2];
sx q[2];
rz(-0.93859276) q[2];
sx q[2];
rz(-0.8559627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2655609) q[1];
sx q[1];
rz(-2.4596493) q[1];
sx q[1];
rz(-2.4443786) q[1];
rz(-pi) q[2];
rz(-1.394519) q[3];
sx q[3];
rz(-0.44197861) q[3];
sx q[3];
rz(-0.37261841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1641757) q[2];
sx q[2];
rz(-1.5693393) q[2];
sx q[2];
rz(0.24641985) q[2];
rz(-2.1988846) q[3];
sx q[3];
rz(-1.0859414) q[3];
sx q[3];
rz(0.76343083) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3281658) q[0];
sx q[0];
rz(-1.8712217) q[0];
sx q[0];
rz(-0.060591977) q[0];
rz(1.5024441) q[1];
sx q[1];
rz(-1.363058) q[1];
sx q[1];
rz(-1.5477808) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12822026) q[0];
sx q[0];
rz(-0.0039847535) q[0];
sx q[0];
rz(0.26040034) q[0];
x q[1];
rz(-2.5663788) q[2];
sx q[2];
rz(-1.4123431) q[2];
sx q[2];
rz(2.2666933) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.036435) q[1];
sx q[1];
rz(-2.1686825) q[1];
sx q[1];
rz(0.5374686) q[1];
x q[2];
rz(-2.4259858) q[3];
sx q[3];
rz(-1.2168435) q[3];
sx q[3];
rz(2.7277814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8760406) q[2];
sx q[2];
rz(-0.92090845) q[2];
sx q[2];
rz(1.6843686) q[2];
rz(3.0217116) q[3];
sx q[3];
rz(-2.9370152) q[3];
sx q[3];
rz(2.1129107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2940755) q[0];
sx q[0];
rz(-2.1471922) q[0];
sx q[0];
rz(1.1353528) q[0];
rz(0.10336939) q[1];
sx q[1];
rz(-1.7366948) q[1];
sx q[1];
rz(0.9476544) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4554868) q[0];
sx q[0];
rz(-0.64160937) q[0];
sx q[0];
rz(3.0296586) q[0];
rz(-pi) q[1];
rz(-2.728711) q[2];
sx q[2];
rz(-0.92501174) q[2];
sx q[2];
rz(1.4476763) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7458742) q[1];
sx q[1];
rz(-1.8367136) q[1];
sx q[1];
rz(2.2661792) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2470948) q[3];
sx q[3];
rz(-0.94493659) q[3];
sx q[3];
rz(-1.8602399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4244298) q[2];
sx q[2];
rz(-0.62817502) q[2];
sx q[2];
rz(-0.90432811) q[2];
rz(-1.8782015) q[3];
sx q[3];
rz(-1.9046013) q[3];
sx q[3];
rz(2.625107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.271027) q[0];
sx q[0];
rz(-2.3834383) q[0];
sx q[0];
rz(2.1562321) q[0];
rz(-2.789978) q[1];
sx q[1];
rz(-2.1814587) q[1];
sx q[1];
rz(0.34686372) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0766633) q[0];
sx q[0];
rz(-2.4682625) q[0];
sx q[0];
rz(-1.2162186) q[0];
x q[1];
rz(-0.11254452) q[2];
sx q[2];
rz(-1.1273317) q[2];
sx q[2];
rz(-2.7751768) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.092791768) q[1];
sx q[1];
rz(-1.6597224) q[1];
sx q[1];
rz(-2.4503178) q[1];
rz(-1.204416) q[3];
sx q[3];
rz(-1.9854331) q[3];
sx q[3];
rz(2.7759027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88006192) q[2];
sx q[2];
rz(-0.71892771) q[2];
sx q[2];
rz(3.0066709) q[2];
rz(3.0123582) q[3];
sx q[3];
rz(-2.7415469) q[3];
sx q[3];
rz(2.1028886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018910949) q[0];
sx q[0];
rz(-1.2553348) q[0];
sx q[0];
rz(-3.0154764) q[0];
rz(-0.46161721) q[1];
sx q[1];
rz(-1.7414265) q[1];
sx q[1];
rz(0.19207676) q[1];
rz(2.1264524) q[2];
sx q[2];
rz(-0.75426741) q[2];
sx q[2];
rz(-1.5744899) q[2];
rz(0.25855385) q[3];
sx q[3];
rz(-1.7333442) q[3];
sx q[3];
rz(2.6556591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
