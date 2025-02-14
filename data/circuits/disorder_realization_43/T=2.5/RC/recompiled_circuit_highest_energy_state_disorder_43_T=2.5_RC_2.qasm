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
rz(0.046959538) q[0];
sx q[0];
rz(1.5927915) q[0];
sx q[0];
rz(11.172063) q[0];
rz(-0.46051639) q[1];
sx q[1];
rz(4.4098201) q[1];
sx q[1];
rz(8.9486651) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9481422) q[0];
sx q[0];
rz(-2.1102437) q[0];
sx q[0];
rz(2.6953117) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23525518) q[2];
sx q[2];
rz(-1.3089184) q[2];
sx q[2];
rz(2.1507598) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.75501498) q[1];
sx q[1];
rz(-1.5534895) q[1];
sx q[1];
rz(-0.14825082) q[1];
x q[2];
rz(-1.5067817) q[3];
sx q[3];
rz(-1.598017) q[3];
sx q[3];
rz(2.343319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32690471) q[2];
sx q[2];
rz(-1.9196332) q[2];
sx q[2];
rz(1.8587221) q[2];
rz(3.0021216) q[3];
sx q[3];
rz(-1.3160416) q[3];
sx q[3];
rz(2.4563346) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4271456) q[0];
sx q[0];
rz(-2.7834263) q[0];
sx q[0];
rz(-2.8156679) q[0];
rz(-2.4320995) q[1];
sx q[1];
rz(-2.8804417) q[1];
sx q[1];
rz(0.66169468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5241338) q[0];
sx q[0];
rz(-1.5781998) q[0];
sx q[0];
rz(-1.5803158) q[0];
rz(2.2108033) q[2];
sx q[2];
rz(-2.402596) q[2];
sx q[2];
rz(1.3698824) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.41971808) q[1];
sx q[1];
rz(-1.9456753) q[1];
sx q[1];
rz(-2.8931395) q[1];
x q[2];
rz(0.68485188) q[3];
sx q[3];
rz(-2.0461444) q[3];
sx q[3];
rz(2.1641987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1208531) q[2];
sx q[2];
rz(-1.8031305) q[2];
sx q[2];
rz(0.2600812) q[2];
rz(2.2777879) q[3];
sx q[3];
rz(-1.6983906) q[3];
sx q[3];
rz(-1.9728194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.23673713) q[0];
sx q[0];
rz(-2.3598292) q[0];
sx q[0];
rz(0.28808638) q[0];
rz(-1.9347363) q[1];
sx q[1];
rz(-2.0286045) q[1];
sx q[1];
rz(-2.3634214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85976493) q[0];
sx q[0];
rz(-0.71548235) q[0];
sx q[0];
rz(-1.2802109) q[0];
rz(0.065041754) q[2];
sx q[2];
rz(-0.75924004) q[2];
sx q[2];
rz(-2.0652079) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6699804) q[1];
sx q[1];
rz(-2.0005352) q[1];
sx q[1];
rz(2.5915036) q[1];
rz(-2.7500626) q[3];
sx q[3];
rz(-0.96116306) q[3];
sx q[3];
rz(-3.0972598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8975767) q[2];
sx q[2];
rz(-1.1268758) q[2];
sx q[2];
rz(2.8934532) q[2];
rz(1.0330307) q[3];
sx q[3];
rz(-1.6525729) q[3];
sx q[3];
rz(-1.2584794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4397864) q[0];
sx q[0];
rz(-0.93455625) q[0];
sx q[0];
rz(2.5820861) q[0];
rz(1.4065546) q[1];
sx q[1];
rz(-2.1969257) q[1];
sx q[1];
rz(-1.9962126) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9713719) q[0];
sx q[0];
rz(-1.6165501) q[0];
sx q[0];
rz(-0.45661033) q[0];
rz(-2.4907095) q[2];
sx q[2];
rz(-1.1090298) q[2];
sx q[2];
rz(-0.24591638) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8662474) q[1];
sx q[1];
rz(-1.2762349) q[1];
sx q[1];
rz(0.7742851) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5942176) q[3];
sx q[3];
rz(-2.5716647) q[3];
sx q[3];
rz(-0.83372926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2650602) q[2];
sx q[2];
rz(-2.127485) q[2];
sx q[2];
rz(-1.8939135) q[2];
rz(-2.0606591) q[3];
sx q[3];
rz(-1.7526046) q[3];
sx q[3];
rz(1.8814253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0697698) q[0];
sx q[0];
rz(-2.6520196) q[0];
sx q[0];
rz(2.3967632) q[0];
rz(-1.7597594) q[1];
sx q[1];
rz(-2.0227183) q[1];
sx q[1];
rz(3.0435069) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25135219) q[0];
sx q[0];
rz(-0.47516631) q[0];
sx q[0];
rz(0.71936468) q[0];
rz(0.44386835) q[2];
sx q[2];
rz(-1.8763246) q[2];
sx q[2];
rz(1.533184) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91589663) q[1];
sx q[1];
rz(-2.4476123) q[1];
sx q[1];
rz(-3.1129614) q[1];
rz(-pi) q[2];
rz(-0.10048203) q[3];
sx q[3];
rz(-2.049787) q[3];
sx q[3];
rz(0.27083957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.44094008) q[2];
sx q[2];
rz(-0.38267371) q[2];
sx q[2];
rz(-1.0901701) q[2];
rz(-2.8160461) q[3];
sx q[3];
rz(-1.5109589) q[3];
sx q[3];
rz(2.8946099) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22661041) q[0];
sx q[0];
rz(-2.60422) q[0];
sx q[0];
rz(-1.9052624) q[0];
rz(-2.5044598) q[1];
sx q[1];
rz(-2.5814711) q[1];
sx q[1];
rz(-2.2083652) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7819311) q[0];
sx q[0];
rz(-2.1225516) q[0];
sx q[0];
rz(-2.7130068) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9901598) q[2];
sx q[2];
rz(-1.7438816) q[2];
sx q[2];
rz(1.2746379) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2843202) q[1];
sx q[1];
rz(-0.30306268) q[1];
sx q[1];
rz(-0.28115718) q[1];
x q[2];
rz(-2.0863462) q[3];
sx q[3];
rz(-1.7988867) q[3];
sx q[3];
rz(-2.4683964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8895662) q[2];
sx q[2];
rz(-1.3875952) q[2];
sx q[2];
rz(1.265556) q[2];
rz(1.7172074) q[3];
sx q[3];
rz(-0.5286743) q[3];
sx q[3];
rz(-0.85844794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76765656) q[0];
sx q[0];
rz(-2.7601384) q[0];
sx q[0];
rz(-0.23442991) q[0];
rz(-2.7081721) q[1];
sx q[1];
rz(-1.0973009) q[1];
sx q[1];
rz(-1.1385328) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94110452) q[0];
sx q[0];
rz(-0.88747665) q[0];
sx q[0];
rz(-2.3462458) q[0];
x q[1];
rz(1.5628596) q[2];
sx q[2];
rz(-1.0801549) q[2];
sx q[2];
rz(2.454941) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93764544) q[1];
sx q[1];
rz(-2.0729985) q[1];
sx q[1];
rz(-0.03003386) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3765923) q[3];
sx q[3];
rz(-0.89316165) q[3];
sx q[3];
rz(1.6226131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66702691) q[2];
sx q[2];
rz(-0.67956769) q[2];
sx q[2];
rz(-0.95728528) q[2];
rz(-2.3840617) q[3];
sx q[3];
rz(-1.6672971) q[3];
sx q[3];
rz(2.9981414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677419) q[0];
sx q[0];
rz(-2.0389281) q[0];
sx q[0];
rz(2.962501) q[0];
rz(0.20009072) q[1];
sx q[1];
rz(-2.4724019) q[1];
sx q[1];
rz(-2.3650513) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0901129) q[0];
sx q[0];
rz(-2.5830373) q[0];
sx q[0];
rz(1.2249951) q[0];
rz(-pi) q[1];
rz(-0.29751038) q[2];
sx q[2];
rz(-2.1665467) q[2];
sx q[2];
rz(1.3551301) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4458865) q[1];
sx q[1];
rz(-1.4974471) q[1];
sx q[1];
rz(-2.785745) q[1];
rz(-pi) q[2];
x q[2];
rz(1.253088) q[3];
sx q[3];
rz(-1.5220257) q[3];
sx q[3];
rz(-1.932285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9658003) q[2];
sx q[2];
rz(-1.0719904) q[2];
sx q[2];
rz(1.3624066) q[2];
rz(1.8386748) q[3];
sx q[3];
rz(-1.4785654) q[3];
sx q[3];
rz(2.1030857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68504828) q[0];
sx q[0];
rz(-1.9467204) q[0];
sx q[0];
rz(3.0273279) q[0];
rz(1.201913) q[1];
sx q[1];
rz(-0.54113954) q[1];
sx q[1];
rz(3.0466383) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9499493) q[0];
sx q[0];
rz(-1.7472634) q[0];
sx q[0];
rz(-1.0383738) q[0];
rz(-pi) q[1];
rz(-1.0509148) q[2];
sx q[2];
rz(-0.78453817) q[2];
sx q[2];
rz(0.21518165) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.488738) q[1];
sx q[1];
rz(-1.514637) q[1];
sx q[1];
rz(-0.55027669) q[1];
rz(-pi) q[2];
rz(2.0801274) q[3];
sx q[3];
rz(-1.9373139) q[3];
sx q[3];
rz(-1.4118005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3007043) q[2];
sx q[2];
rz(-1.3450832) q[2];
sx q[2];
rz(-2.4697206) q[2];
rz(0.68615174) q[3];
sx q[3];
rz(-0.66303623) q[3];
sx q[3];
rz(-1.9239976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9753863) q[0];
sx q[0];
rz(-2.9500742) q[0];
sx q[0];
rz(-1.6084877) q[0];
rz(2.8168822) q[1];
sx q[1];
rz(-1.8638116) q[1];
sx q[1];
rz(0.3130354) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9444549) q[0];
sx q[0];
rz(-2.1416683) q[0];
sx q[0];
rz(3.0464493) q[0];
rz(-pi) q[1];
rz(1.3029477) q[2];
sx q[2];
rz(-0.22780475) q[2];
sx q[2];
rz(2.0382263) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5394143) q[1];
sx q[1];
rz(-1.0262607) q[1];
sx q[1];
rz(-0.79411749) q[1];
rz(-pi) q[2];
rz(-1.1129727) q[3];
sx q[3];
rz(-0.99689299) q[3];
sx q[3];
rz(-1.7750268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93840557) q[2];
sx q[2];
rz(-1.9878191) q[2];
sx q[2];
rz(-1.3333092) q[2];
rz(-0.61776727) q[3];
sx q[3];
rz(-2.1963162) q[3];
sx q[3];
rz(1.1355404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4278605) q[0];
sx q[0];
rz(-1.8129616) q[0];
sx q[0];
rz(0.34995361) q[0];
rz(-2.2433544) q[1];
sx q[1];
rz(-2.0571092) q[1];
sx q[1];
rz(2.7877997) q[1];
rz(-2.4785715) q[2];
sx q[2];
rz(-0.86227476) q[2];
sx q[2];
rz(-3.0515565) q[2];
rz(-3.114936) q[3];
sx q[3];
rz(-2.4116357) q[3];
sx q[3];
rz(0.38105376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
