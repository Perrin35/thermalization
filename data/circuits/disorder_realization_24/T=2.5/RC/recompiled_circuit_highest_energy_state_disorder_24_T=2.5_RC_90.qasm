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
rz(-0.51198045) q[0];
sx q[0];
rz(-2.5706302) q[0];
sx q[0];
rz(-2.9538739) q[0];
rz(0.81387782) q[1];
sx q[1];
rz(-1.6698281) q[1];
sx q[1];
rz(-1.8522813) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0388237) q[0];
sx q[0];
rz(-1.7784405) q[0];
sx q[0];
rz(0.24216346) q[0];
x q[1];
rz(-0.046229049) q[2];
sx q[2];
rz(-2.5704489) q[2];
sx q[2];
rz(2.2276304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6179168) q[1];
sx q[1];
rz(-2.6111228) q[1];
sx q[1];
rz(-0.50650774) q[1];
rz(-pi) q[2];
rz(-1.0008775) q[3];
sx q[3];
rz(-2.7288247) q[3];
sx q[3];
rz(1.4980396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2671555) q[2];
sx q[2];
rz(-2.0384553) q[2];
sx q[2];
rz(-0.77679408) q[2];
rz(1.3022425) q[3];
sx q[3];
rz(-1.4812508) q[3];
sx q[3];
rz(1.2031901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60748196) q[0];
sx q[0];
rz(-2.6597839) q[0];
sx q[0];
rz(-0.99824655) q[0];
rz(-0.36969319) q[1];
sx q[1];
rz(-1.0414711) q[1];
sx q[1];
rz(1.1179771) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8449865) q[0];
sx q[0];
rz(-0.70195508) q[0];
sx q[0];
rz(0.92092307) q[0];
x q[1];
rz(0.67310272) q[2];
sx q[2];
rz(-2.3323698) q[2];
sx q[2];
rz(-0.26462091) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3433132) q[1];
sx q[1];
rz(-1.4237397) q[1];
sx q[1];
rz(-1.1215854) q[1];
rz(-0.26882986) q[3];
sx q[3];
rz(-0.95655609) q[3];
sx q[3];
rz(-1.5419095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.86576858) q[2];
sx q[2];
rz(-1.3752702) q[2];
sx q[2];
rz(0.38530525) q[2];
rz(3.1028808) q[3];
sx q[3];
rz(-0.35274115) q[3];
sx q[3];
rz(0.64046162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5019219) q[0];
sx q[0];
rz(-2.6646035) q[0];
sx q[0];
rz(1.3013526) q[0];
rz(-3.0838857) q[1];
sx q[1];
rz(-2.6951908) q[1];
sx q[1];
rz(-1.3267964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037687) q[0];
sx q[0];
rz(-2.38842) q[0];
sx q[0];
rz(-0.84363787) q[0];
rz(1.7102555) q[2];
sx q[2];
rz(-1.4118952) q[2];
sx q[2];
rz(3.0099843) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4918265) q[1];
sx q[1];
rz(-0.61188243) q[1];
sx q[1];
rz(-0.51458451) q[1];
x q[2];
rz(-0.84872021) q[3];
sx q[3];
rz(-1.3009138) q[3];
sx q[3];
rz(0.11972846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79746276) q[2];
sx q[2];
rz(-0.8351438) q[2];
sx q[2];
rz(0.75378913) q[2];
rz(-1.3695184) q[3];
sx q[3];
rz(-1.6794208) q[3];
sx q[3];
rz(0.65521017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21514431) q[0];
sx q[0];
rz(-2.0560052) q[0];
sx q[0];
rz(1.3947067) q[0];
rz(1.9649547) q[1];
sx q[1];
rz(-1.5086915) q[1];
sx q[1];
rz(-2.7746157) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0622371) q[0];
sx q[0];
rz(-1.5141762) q[0];
sx q[0];
rz(1.6376154) q[0];
rz(-pi) q[1];
rz(0.33432482) q[2];
sx q[2];
rz(-1.8211289) q[2];
sx q[2];
rz(1.804753) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1550459) q[1];
sx q[1];
rz(-1.2147702) q[1];
sx q[1];
rz(1.3638391) q[1];
x q[2];
rz(-0.36716299) q[3];
sx q[3];
rz(-1.6510186) q[3];
sx q[3];
rz(1.0368766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7299812) q[2];
sx q[2];
rz(-1.9133762) q[2];
sx q[2];
rz(-1.8709987) q[2];
rz(0.96108428) q[3];
sx q[3];
rz(-1.4189439) q[3];
sx q[3];
rz(0.050475033) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5807895) q[0];
sx q[0];
rz(-0.92629782) q[0];
sx q[0];
rz(1.8286937) q[0];
rz(1.1543697) q[1];
sx q[1];
rz(-1.9998974) q[1];
sx q[1];
rz(-0.55783522) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9168222) q[0];
sx q[0];
rz(-2.0799412) q[0];
sx q[0];
rz(-2.1780685) q[0];
rz(0.16061546) q[2];
sx q[2];
rz(-1.4945507) q[2];
sx q[2];
rz(0.74197021) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.055649672) q[1];
sx q[1];
rz(-0.80971566) q[1];
sx q[1];
rz(1.147406) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6799503) q[3];
sx q[3];
rz(-1.7064693) q[3];
sx q[3];
rz(1.6035348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.64342) q[2];
sx q[2];
rz(-1.2926481) q[2];
sx q[2];
rz(0.84929973) q[2];
rz(0.15527209) q[3];
sx q[3];
rz(-2.1657491) q[3];
sx q[3];
rz(-0.37080216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7840541) q[0];
sx q[0];
rz(-2.5548866) q[0];
sx q[0];
rz(0.46911711) q[0];
rz(0.47438374) q[1];
sx q[1];
rz(-0.7862888) q[1];
sx q[1];
rz(0.75538409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3123388) q[0];
sx q[0];
rz(-2.8702998) q[0];
sx q[0];
rz(1.5543227) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7003661) q[2];
sx q[2];
rz(-0.6775113) q[2];
sx q[2];
rz(-0.73523075) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50690097) q[1];
sx q[1];
rz(-1.4441057) q[1];
sx q[1];
rz(1.7853398) q[1];
rz(-1.1047885) q[3];
sx q[3];
rz(-2.3475523) q[3];
sx q[3];
rz(-0.023762437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0581806) q[2];
sx q[2];
rz(-1.3231134) q[2];
sx q[2];
rz(2.9537436) q[2];
rz(1.8958873) q[3];
sx q[3];
rz(-1.3789504) q[3];
sx q[3];
rz(1.0904306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2096527) q[0];
sx q[0];
rz(-1.8745475) q[0];
sx q[0];
rz(2.7506822) q[0];
rz(-2.2115425) q[1];
sx q[1];
rz(-1.7101733) q[1];
sx q[1];
rz(1.8375058) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6041038) q[0];
sx q[0];
rz(-1.9216282) q[0];
sx q[0];
rz(-3.119191) q[0];
x q[1];
rz(1.485059) q[2];
sx q[2];
rz(-1.9112327) q[2];
sx q[2];
rz(2.1742333) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.48941055) q[1];
sx q[1];
rz(-0.48508185) q[1];
sx q[1];
rz(1.9368383) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29629064) q[3];
sx q[3];
rz(-1.5888831) q[3];
sx q[3];
rz(-2.2294105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2685252) q[2];
sx q[2];
rz(-2.8748942) q[2];
sx q[2];
rz(3.0282057) q[2];
rz(-0.015965613) q[3];
sx q[3];
rz(-1.6458052) q[3];
sx q[3];
rz(-2.9769843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3575344) q[0];
sx q[0];
rz(-1.6679732) q[0];
sx q[0];
rz(-3.0175324) q[0];
rz(3.0198174) q[1];
sx q[1];
rz(-0.75141326) q[1];
sx q[1];
rz(1.442499) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.34262) q[0];
sx q[0];
rz(-2.7443287) q[0];
sx q[0];
rz(1.8333866) q[0];
x q[1];
rz(1.5273483) q[2];
sx q[2];
rz(-1.279131) q[2];
sx q[2];
rz(-1.0857605) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.6119163) q[1];
sx q[1];
rz(-1.2280591) q[1];
sx q[1];
rz(0.14118282) q[1];
rz(-0.7300709) q[3];
sx q[3];
rz(-1.7582809) q[3];
sx q[3];
rz(1.2227675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1759935) q[2];
sx q[2];
rz(-1.7815353) q[2];
sx q[2];
rz(0.74481258) q[2];
rz(2.2143769) q[3];
sx q[3];
rz(-1.5085446) q[3];
sx q[3];
rz(-2.0492699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6532779) q[0];
sx q[0];
rz(-1.4396242) q[0];
sx q[0];
rz(-2.9162245) q[0];
rz(2.0291746) q[1];
sx q[1];
rz(-2.367159) q[1];
sx q[1];
rz(-2.2427799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32466896) q[0];
sx q[0];
rz(-0.82624871) q[0];
sx q[0];
rz(-2.3284495) q[0];
rz(1.1583011) q[2];
sx q[2];
rz(-1.0935244) q[2];
sx q[2];
rz(2.0606747) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.73279335) q[1];
sx q[1];
rz(-2.2299754) q[1];
sx q[1];
rz(0.2562457) q[1];
rz(-pi) q[2];
rz(1.6055272) q[3];
sx q[3];
rz(-0.29303778) q[3];
sx q[3];
rz(-1.7188096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0972458) q[2];
sx q[2];
rz(-1.711742) q[2];
sx q[2];
rz(-2.782235) q[2];
rz(2.8906631) q[3];
sx q[3];
rz(-2.0693306) q[3];
sx q[3];
rz(-0.76771626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74442416) q[0];
sx q[0];
rz(-3.011062) q[0];
sx q[0];
rz(0.8771483) q[0];
rz(-1.5769222) q[1];
sx q[1];
rz(-1.501187) q[1];
sx q[1];
rz(-2.3559779) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1164074) q[0];
sx q[0];
rz(-1.3337564) q[0];
sx q[0];
rz(-2.9171506) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4841275) q[2];
sx q[2];
rz(-0.6415002) q[2];
sx q[2];
rz(2.7335707) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4686376) q[1];
sx q[1];
rz(-1.6412927) q[1];
sx q[1];
rz(2.7967909) q[1];
x q[2];
rz(-1.4203868) q[3];
sx q[3];
rz(-2.0724943) q[3];
sx q[3];
rz(-3.0034163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2239573) q[2];
sx q[2];
rz(-0.79006299) q[2];
sx q[2];
rz(0.7381953) q[2];
rz(-0.92292845) q[3];
sx q[3];
rz(-1.3083369) q[3];
sx q[3];
rz(2.8452528) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34687635) q[0];
sx q[0];
rz(-1.3608169) q[0];
sx q[0];
rz(0.97186744) q[0];
rz(0.90732668) q[1];
sx q[1];
rz(-1.0846039) q[1];
sx q[1];
rz(-0.32122282) q[1];
rz(2.1335294) q[2];
sx q[2];
rz(-0.65958644) q[2];
sx q[2];
rz(1.6749254) q[2];
rz(1.0492269) q[3];
sx q[3];
rz(-1.0913625) q[3];
sx q[3];
rz(2.8186295) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
