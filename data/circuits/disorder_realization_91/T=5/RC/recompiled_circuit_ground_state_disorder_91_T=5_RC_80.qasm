OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.087091669) q[0];
sx q[0];
rz(-2.2890685) q[0];
sx q[0];
rz(2.8791715) q[0];
rz(2.8334795) q[1];
sx q[1];
rz(-2.3051655) q[1];
sx q[1];
rz(-0.0094553789) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.981009) q[0];
sx q[0];
rz(-2.0045337) q[0];
sx q[0];
rz(2.1061312) q[0];
x q[1];
rz(0.043536206) q[2];
sx q[2];
rz(-1.5015209) q[2];
sx q[2];
rz(-2.2486698) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6188742) q[1];
sx q[1];
rz(-2.1565003) q[1];
sx q[1];
rz(-2.5453091) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49008781) q[3];
sx q[3];
rz(-2.2667829) q[3];
sx q[3];
rz(2.0139607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0414163) q[2];
sx q[2];
rz(-0.88465038) q[2];
sx q[2];
rz(-0.46119383) q[2];
rz(-0.50208107) q[3];
sx q[3];
rz(-2.3177948) q[3];
sx q[3];
rz(-1.4095149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13040386) q[0];
sx q[0];
rz(-1.6809502) q[0];
sx q[0];
rz(-2.9050264) q[0];
rz(-2.323281) q[1];
sx q[1];
rz(-1.2032443) q[1];
sx q[1];
rz(-2.1990105) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078694669) q[0];
sx q[0];
rz(-1.601614) q[0];
sx q[0];
rz(-0.0027058733) q[0];
rz(-0.06748345) q[2];
sx q[2];
rz(-1.4500256) q[2];
sx q[2];
rz(-1.6492594) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.67737867) q[1];
sx q[1];
rz(-2.3056952) q[1];
sx q[1];
rz(-1.5688495) q[1];
rz(-pi) q[2];
rz(0.97901042) q[3];
sx q[3];
rz(-2.5755582) q[3];
sx q[3];
rz(-2.4477959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4054823) q[2];
sx q[2];
rz(-0.63850275) q[2];
sx q[2];
rz(2.45641) q[2];
rz(0.80125609) q[3];
sx q[3];
rz(-1.6359436) q[3];
sx q[3];
rz(1.8906458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41246978) q[0];
sx q[0];
rz(-0.47054371) q[0];
sx q[0];
rz(0.98010081) q[0];
rz(-3.0209172) q[1];
sx q[1];
rz(-1.1680892) q[1];
sx q[1];
rz(1.4792222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2187933) q[0];
sx q[0];
rz(-0.81577071) q[0];
sx q[0];
rz(1.7936617) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0329887) q[2];
sx q[2];
rz(-2.0304907) q[2];
sx q[2];
rz(2.5407956) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9514346) q[1];
sx q[1];
rz(-0.010316545) q[1];
sx q[1];
rz(-2.6134788) q[1];
rz(-pi) q[2];
rz(-2.6605786) q[3];
sx q[3];
rz(-0.64669197) q[3];
sx q[3];
rz(2.8493136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4270619) q[2];
sx q[2];
rz(-1.4947299) q[2];
sx q[2];
rz(-0.18043268) q[2];
rz(0.42029542) q[3];
sx q[3];
rz(-0.97729483) q[3];
sx q[3];
rz(-0.54496566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.9412398) q[0];
sx q[0];
rz(-1.350133) q[0];
sx q[0];
rz(3.0923162) q[0];
rz(0.73206466) q[1];
sx q[1];
rz(-0.8492291) q[1];
sx q[1];
rz(-1.3759618) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9889209) q[0];
sx q[0];
rz(-1.3878569) q[0];
sx q[0];
rz(1.4519917) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3591197) q[2];
sx q[2];
rz(-2.3251901) q[2];
sx q[2];
rz(-2.8704081) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.84280864) q[1];
sx q[1];
rz(-2.7453303) q[1];
sx q[1];
rz(2.9664459) q[1];
rz(-pi) q[2];
rz(-0.64582156) q[3];
sx q[3];
rz(-2.1649515) q[3];
sx q[3];
rz(-0.37120968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.10776082) q[2];
sx q[2];
rz(-1.3904479) q[2];
sx q[2];
rz(-0.66229171) q[2];
rz(-1.8108588) q[3];
sx q[3];
rz(-0.8225421) q[3];
sx q[3];
rz(-1.8432157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0996284) q[0];
sx q[0];
rz(-2.8350416) q[0];
sx q[0];
rz(2.1339259) q[0];
rz(0.10534605) q[1];
sx q[1];
rz(-0.65709972) q[1];
sx q[1];
rz(2.9109921) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5916633) q[0];
sx q[0];
rz(-1.9029034) q[0];
sx q[0];
rz(-1.9426533) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2097539) q[2];
sx q[2];
rz(-1.7954383) q[2];
sx q[2];
rz(1.063907) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19781348) q[1];
sx q[1];
rz(-0.18584419) q[1];
sx q[1];
rz(-0.4798823) q[1];
rz(0.19959656) q[3];
sx q[3];
rz(-0.65465876) q[3];
sx q[3];
rz(-2.9195291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0838123) q[2];
sx q[2];
rz(-2.1746706) q[2];
sx q[2];
rz(-2.9528565) q[2];
rz(1.2356637) q[3];
sx q[3];
rz(-2.6344447) q[3];
sx q[3];
rz(-1.88131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.1756303) q[0];
sx q[0];
rz(-1.9251134) q[0];
sx q[0];
rz(-0.29773444) q[0];
rz(3.0689902) q[1];
sx q[1];
rz(-2.4563792) q[1];
sx q[1];
rz(-2.4593478) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84083623) q[0];
sx q[0];
rz(-1.2781236) q[0];
sx q[0];
rz(-2.2038644) q[0];
x q[1];
rz(2.2043292) q[2];
sx q[2];
rz(-1.6906977) q[2];
sx q[2];
rz(-1.0064358) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.678735) q[1];
sx q[1];
rz(-1.917352) q[1];
sx q[1];
rz(-2.9842675) q[1];
x q[2];
rz(2.7704084) q[3];
sx q[3];
rz(-2.1371578) q[3];
sx q[3];
rz(1.0584099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5785152) q[2];
sx q[2];
rz(-2.7882521) q[2];
sx q[2];
rz(-3.0541218) q[2];
rz(2.2443917) q[3];
sx q[3];
rz(-1.8950491) q[3];
sx q[3];
rz(-2.7520666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0422269) q[0];
sx q[0];
rz(-1.1023738) q[0];
sx q[0];
rz(-1.9788096) q[0];
rz(0.21576628) q[1];
sx q[1];
rz(-2.3240604) q[1];
sx q[1];
rz(-0.93528265) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1638086) q[0];
sx q[0];
rz(-0.80923128) q[0];
sx q[0];
rz(0.33235312) q[0];
rz(-0.14752702) q[2];
sx q[2];
rz(-1.3152244) q[2];
sx q[2];
rz(-2.7856261) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.062309655) q[1];
sx q[1];
rz(-0.78719646) q[1];
sx q[1];
rz(1.4592378) q[1];
rz(-pi) q[2];
rz(-1.2452081) q[3];
sx q[3];
rz(-1.8703484) q[3];
sx q[3];
rz(1.3249782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0995348) q[2];
sx q[2];
rz(-0.65379405) q[2];
sx q[2];
rz(0.81134861) q[2];
rz(-0.30284303) q[3];
sx q[3];
rz(-1.1648014) q[3];
sx q[3];
rz(-0.6306878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52043668) q[0];
sx q[0];
rz(-2.8148837) q[0];
sx q[0];
rz(-0.69123554) q[0];
rz(-0.65903819) q[1];
sx q[1];
rz(-1.6233416) q[1];
sx q[1];
rz(-0.59115994) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6499929) q[0];
sx q[0];
rz(-2.8584456) q[0];
sx q[0];
rz(1.5294019) q[0];
rz(-pi) q[1];
rz(-1.0760087) q[2];
sx q[2];
rz(-2.587834) q[2];
sx q[2];
rz(-0.29386917) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4860079) q[1];
sx q[1];
rz(-2.297245) q[1];
sx q[1];
rz(-0.4419216) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.05686111) q[3];
sx q[3];
rz(-1.2788075) q[3];
sx q[3];
rz(2.9546471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9968694) q[2];
sx q[2];
rz(-0.99088061) q[2];
sx q[2];
rz(1.4608176) q[2];
rz(-2.3878494) q[3];
sx q[3];
rz(-1.6921348) q[3];
sx q[3];
rz(-2.6678705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6853471) q[0];
sx q[0];
rz(-1.101475) q[0];
sx q[0];
rz(-0.22612485) q[0];
rz(1.6926758) q[1];
sx q[1];
rz(-0.96335226) q[1];
sx q[1];
rz(1.0015782) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152214) q[0];
sx q[0];
rz(-2.5363521) q[0];
sx q[0];
rz(2.8855118) q[0];
rz(-1.1923157) q[2];
sx q[2];
rz(-0.90736249) q[2];
sx q[2];
rz(-2.0326234) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7563419) q[1];
sx q[1];
rz(-1.4651734) q[1];
sx q[1];
rz(-2.9526605) q[1];
rz(-pi) q[2];
rz(1.8705192) q[3];
sx q[3];
rz(-2.2329438) q[3];
sx q[3];
rz(-1.3817092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61233026) q[2];
sx q[2];
rz(-0.77028209) q[2];
sx q[2];
rz(-2.9840577) q[2];
rz(-2.9869288) q[3];
sx q[3];
rz(-1.1636584) q[3];
sx q[3];
rz(-0.4304339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73151076) q[0];
sx q[0];
rz(-1.2393351) q[0];
sx q[0];
rz(0.70247689) q[0];
rz(-0.53600535) q[1];
sx q[1];
rz(-0.82967007) q[1];
sx q[1];
rz(-1.747267) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4634906) q[0];
sx q[0];
rz(-0.16568389) q[0];
sx q[0];
rz(1.7278306) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.093981103) q[2];
sx q[2];
rz(-1.4620442) q[2];
sx q[2];
rz(-2.5861339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.179217) q[1];
sx q[1];
rz(-2.0834647) q[1];
sx q[1];
rz(2.8151399) q[1];
x q[2];
rz(-1.570618) q[3];
sx q[3];
rz(-1.1552253) q[3];
sx q[3];
rz(-0.10997575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.006762) q[2];
sx q[2];
rz(-2.4595021) q[2];
sx q[2];
rz(-1.6949867) q[2];
rz(-0.26099482) q[3];
sx q[3];
rz(-2.1311396) q[3];
sx q[3];
rz(1.235777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21988729) q[0];
sx q[0];
rz(-0.6548665) q[0];
sx q[0];
rz(-1.0153216) q[0];
rz(1.075853) q[1];
sx q[1];
rz(-1.8668108) q[1];
sx q[1];
rz(1.6783953) q[1];
rz(-1.6816611) q[2];
sx q[2];
rz(-1.1959497) q[2];
sx q[2];
rz(-0.49899286) q[2];
rz(2.0113284) q[3];
sx q[3];
rz(-2.2838099) q[3];
sx q[3];
rz(1.6259307) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
