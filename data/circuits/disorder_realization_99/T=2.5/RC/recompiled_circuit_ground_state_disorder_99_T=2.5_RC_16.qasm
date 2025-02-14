OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34243256) q[0];
sx q[0];
rz(2.7437796) q[0];
sx q[0];
rz(7.5510511) q[0];
rz(-0.59029382) q[1];
sx q[1];
rz(-0.78650147) q[1];
sx q[1];
rz(-1.2193349) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0578227) q[0];
sx q[0];
rz(-1.738213) q[0];
sx q[0];
rz(2.1517702) q[0];
x q[1];
rz(0.54852416) q[2];
sx q[2];
rz(-1.7146352) q[2];
sx q[2];
rz(1.583969) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.089069033) q[1];
sx q[1];
rz(-2.7929651) q[1];
sx q[1];
rz(-2.2050956) q[1];
rz(-1.8699617) q[3];
sx q[3];
rz(-2.4623021) q[3];
sx q[3];
rz(-0.75377476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5775602) q[2];
sx q[2];
rz(-1.5215678) q[2];
sx q[2];
rz(1.7731898) q[2];
rz(0.91156256) q[3];
sx q[3];
rz(-1.7715958) q[3];
sx q[3];
rz(3.1172359) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0381222) q[0];
sx q[0];
rz(-0.39762527) q[0];
sx q[0];
rz(3.0225515) q[0];
rz(-1.1301522) q[1];
sx q[1];
rz(-2.1739013) q[1];
sx q[1];
rz(-1.4281323) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9790065) q[0];
sx q[0];
rz(-1.9784728) q[0];
sx q[0];
rz(-2.5752221) q[0];
rz(-1.9910286) q[2];
sx q[2];
rz(-1.2017851) q[2];
sx q[2];
rz(-2.8431161) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.94257689) q[1];
sx q[1];
rz(-1.7046832) q[1];
sx q[1];
rz(-1.1034637) q[1];
rz(-0.30562206) q[3];
sx q[3];
rz(-1.247331) q[3];
sx q[3];
rz(-2.4572069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6156562) q[2];
sx q[2];
rz(-1.6776513) q[2];
sx q[2];
rz(0.76457912) q[2];
rz(2.6584451) q[3];
sx q[3];
rz(-1.8254447) q[3];
sx q[3];
rz(2.29276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1044384) q[0];
sx q[0];
rz(-0.26532441) q[0];
sx q[0];
rz(-1.6695439) q[0];
rz(-0.56023359) q[1];
sx q[1];
rz(-1.3872223) q[1];
sx q[1];
rz(1.701042) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36828968) q[0];
sx q[0];
rz(-1.76923) q[0];
sx q[0];
rz(-1.3318568) q[0];
rz(-pi) q[1];
rz(-1.5650827) q[2];
sx q[2];
rz(-2.2457321) q[2];
sx q[2];
rz(0.99944464) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5270823) q[1];
sx q[1];
rz(-2.0918814) q[1];
sx q[1];
rz(-2.7481542) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7974924) q[3];
sx q[3];
rz(-1.8009225) q[3];
sx q[3];
rz(-3.074444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4855839) q[2];
sx q[2];
rz(-1.6501082) q[2];
sx q[2];
rz(-2.7991925) q[2];
rz(-1.7229236) q[3];
sx q[3];
rz(-2.4125621) q[3];
sx q[3];
rz(0.5595783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(2.2320084) q[0];
sx q[0];
rz(-2.7689458) q[0];
sx q[0];
rz(1.5268071) q[0];
rz(-2.3341663) q[1];
sx q[1];
rz(-1.5680983) q[1];
sx q[1];
rz(1.2015013) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75185173) q[0];
sx q[0];
rz(-2.3065615) q[0];
sx q[0];
rz(0.40919183) q[0];
rz(2.4929659) q[2];
sx q[2];
rz(-1.0264215) q[2];
sx q[2];
rz(2.8971162) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2953342) q[1];
sx q[1];
rz(-2.3148119) q[1];
sx q[1];
rz(1.0215525) q[1];
x q[2];
rz(2.8852984) q[3];
sx q[3];
rz(-1.6209037) q[3];
sx q[3];
rz(0.20370349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3379007) q[2];
sx q[2];
rz(-2.1521229) q[2];
sx q[2];
rz(1.2376415) q[2];
rz(-1.3112274) q[3];
sx q[3];
rz(-1.5179736) q[3];
sx q[3];
rz(1.1782014) q[3];
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
rz(1.118498) q[0];
sx q[0];
rz(-1.7685522) q[0];
sx q[0];
rz(2.232724) q[0];
rz(-2.2591649) q[1];
sx q[1];
rz(-2.4274223) q[1];
sx q[1];
rz(-0.73208255) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1650989) q[0];
sx q[0];
rz(-1.1015176) q[0];
sx q[0];
rz(-0.83013541) q[0];
x q[1];
rz(-2.8516099) q[2];
sx q[2];
rz(-0.78943832) q[2];
sx q[2];
rz(-1.0944927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4914507) q[1];
sx q[1];
rz(-2.53778) q[1];
sx q[1];
rz(0.73765386) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77712147) q[3];
sx q[3];
rz(-2.2279871) q[3];
sx q[3];
rz(1.4025721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0714134) q[2];
sx q[2];
rz(-0.81255239) q[2];
sx q[2];
rz(1.8522235) q[2];
rz(-1.4687126) q[3];
sx q[3];
rz(-1.9332935) q[3];
sx q[3];
rz(-0.41951352) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5656723) q[0];
sx q[0];
rz(-1.9603632) q[0];
sx q[0];
rz(2.0400203) q[0];
rz(-2.9167602) q[1];
sx q[1];
rz(-1.79554) q[1];
sx q[1];
rz(-2.1563931) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.824788) q[0];
sx q[0];
rz(-1.7508306) q[0];
sx q[0];
rz(2.2145674) q[0];
rz(-1.9198138) q[2];
sx q[2];
rz(-2.1563081) q[2];
sx q[2];
rz(1.3427918) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8781153) q[1];
sx q[1];
rz(-2.2190071) q[1];
sx q[1];
rz(-0.99261673) q[1];
rz(0.52805488) q[3];
sx q[3];
rz(-1.5925358) q[3];
sx q[3];
rz(0.69247144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3605986) q[2];
sx q[2];
rz(-2.4807319) q[2];
sx q[2];
rz(1.9471656) q[2];
rz(3.0418975) q[3];
sx q[3];
rz(-1.7710779) q[3];
sx q[3];
rz(2.1626507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4286956) q[0];
sx q[0];
rz(-2.6850057) q[0];
sx q[0];
rz(-2.0275443) q[0];
rz(-1.7440589) q[1];
sx q[1];
rz(-1.3797398) q[1];
sx q[1];
rz(-2.8278415) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6614051) q[0];
sx q[0];
rz(-1.5208986) q[0];
sx q[0];
rz(-1.8636835) q[0];
rz(-2.9767562) q[2];
sx q[2];
rz(-1.7007174) q[2];
sx q[2];
rz(-0.31229737) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.113768) q[1];
sx q[1];
rz(-1.8701487) q[1];
sx q[1];
rz(-2.37948) q[1];
rz(-pi) q[2];
rz(2.4558732) q[3];
sx q[3];
rz(-1.7758533) q[3];
sx q[3];
rz(-3.0577537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.067001192) q[2];
sx q[2];
rz(-0.40412298) q[2];
sx q[2];
rz(-0.60603777) q[2];
rz(-2.0634985) q[3];
sx q[3];
rz(-0.86229101) q[3];
sx q[3];
rz(1.3585453) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17230497) q[0];
sx q[0];
rz(-0.91738874) q[0];
sx q[0];
rz(1.1267927) q[0];
rz(0.91398319) q[1];
sx q[1];
rz(-0.4387478) q[1];
sx q[1];
rz(-2.9235358) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870616) q[0];
sx q[0];
rz(-2.0045223) q[0];
sx q[0];
rz(-1.0591255) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3534691) q[2];
sx q[2];
rz(-2.1160853) q[2];
sx q[2];
rz(0.1696378) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7331024) q[1];
sx q[1];
rz(-0.6576076) q[1];
sx q[1];
rz(-1.7773184) q[1];
rz(2.2770834) q[3];
sx q[3];
rz(-0.1801404) q[3];
sx q[3];
rz(0.064398191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6286991) q[2];
sx q[2];
rz(-0.67535526) q[2];
sx q[2];
rz(-0.2891573) q[2];
rz(2.1469927) q[3];
sx q[3];
rz(-1.1887487) q[3];
sx q[3];
rz(2.6836256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7116123) q[0];
sx q[0];
rz(-1.0884322) q[0];
sx q[0];
rz(-0.76643884) q[0];
rz(-1.537716) q[1];
sx q[1];
rz(-1.8285373) q[1];
sx q[1];
rz(-2.2235353) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9452451) q[0];
sx q[0];
rz(-2.1736988) q[0];
sx q[0];
rz(0.57749282) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4337344) q[2];
sx q[2];
rz(-0.22032693) q[2];
sx q[2];
rz(-0.79263055) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2868216) q[1];
sx q[1];
rz(-2.3769551) q[1];
sx q[1];
rz(0.36061339) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10580801) q[3];
sx q[3];
rz(-1.322016) q[3];
sx q[3];
rz(0.63423587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90091577) q[2];
sx q[2];
rz(-2.675481) q[2];
sx q[2];
rz(-0.051699836) q[2];
rz(2.2895571) q[3];
sx q[3];
rz(-1.7107191) q[3];
sx q[3];
rz(-2.8813072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8155415) q[0];
sx q[0];
rz(-1.4651848) q[0];
sx q[0];
rz(-3.0503804) q[0];
rz(-2.3444029) q[1];
sx q[1];
rz(-1.7557095) q[1];
sx q[1];
rz(0.43112722) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9190311) q[0];
sx q[0];
rz(-3.0662144) q[0];
sx q[0];
rz(2.4576709) q[0];
rz(-pi) q[1];
rz(2.2195039) q[2];
sx q[2];
rz(-1.6493634) q[2];
sx q[2];
rz(-2.9804413) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9111002) q[1];
sx q[1];
rz(-0.5020895) q[1];
sx q[1];
rz(-2.6761495) q[1];
rz(-pi) q[2];
rz(-2.1608716) q[3];
sx q[3];
rz(-2.4517165) q[3];
sx q[3];
rz(-1.4341314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5690696) q[2];
sx q[2];
rz(-1.5745682) q[2];
sx q[2];
rz(-1.8756867) q[2];
rz(-1.2931394) q[3];
sx q[3];
rz(-1.061941) q[3];
sx q[3];
rz(-0.40614793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-3.1311998) q[0];
sx q[0];
rz(-2.4507903) q[0];
sx q[0];
rz(-2.3660085) q[0];
rz(0.56397437) q[1];
sx q[1];
rz(-2.0717944) q[1];
sx q[1];
rz(2.0800128) q[1];
rz(1.3971267) q[2];
sx q[2];
rz(-0.43890719) q[2];
sx q[2];
rz(-0.70499805) q[2];
rz(-0.98714491) q[3];
sx q[3];
rz(-2.5844283) q[3];
sx q[3];
rz(1.9682932) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
