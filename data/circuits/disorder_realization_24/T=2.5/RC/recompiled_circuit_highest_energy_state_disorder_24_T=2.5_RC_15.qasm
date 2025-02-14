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
rz(-2.3277148) q[1];
sx q[1];
rz(-1.4717646) q[1];
sx q[1];
rz(-1.2893113) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5587414) q[0];
sx q[0];
rz(-1.8076573) q[0];
sx q[0];
rz(1.7844957) q[0];
rz(3.0953636) q[2];
sx q[2];
rz(-0.57114375) q[2];
sx q[2];
rz(-2.2276304) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.047886176) q[1];
sx q[1];
rz(-2.0290861) q[1];
sx q[1];
rz(1.2935733) q[1];
rz(-1.2175519) q[3];
sx q[3];
rz(-1.3526256) q[3];
sx q[3];
rz(-0.45807236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2671555) q[2];
sx q[2];
rz(-2.0384553) q[2];
sx q[2];
rz(-0.77679408) q[2];
rz(-1.3022425) q[3];
sx q[3];
rz(-1.4812508) q[3];
sx q[3];
rz(1.9384025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60748196) q[0];
sx q[0];
rz(-2.6597839) q[0];
sx q[0];
rz(-0.99824655) q[0];
rz(0.36969319) q[1];
sx q[1];
rz(-2.1001215) q[1];
sx q[1];
rz(1.1179771) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6551483) q[0];
sx q[0];
rz(-2.110743) q[0];
sx q[0];
rz(-0.47294954) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68685617) q[2];
sx q[2];
rz(-2.0389028) q[2];
sx q[2];
rz(1.3324225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4395968) q[1];
sx q[1];
rz(-1.1267822) q[1];
sx q[1];
rz(0.16298144) q[1];
x q[2];
rz(1.2106154) q[3];
sx q[3];
rz(-0.66347117) q[3];
sx q[3];
rz(-2.0455895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86576858) q[2];
sx q[2];
rz(-1.7663225) q[2];
sx q[2];
rz(0.38530525) q[2];
rz(0.038711874) q[3];
sx q[3];
rz(-2.7888515) q[3];
sx q[3];
rz(-2.501131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63967079) q[0];
sx q[0];
rz(-0.47698912) q[0];
sx q[0];
rz(1.3013526) q[0];
rz(3.0838857) q[1];
sx q[1];
rz(-0.4464018) q[1];
sx q[1];
rz(1.8147963) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0327137) q[0];
sx q[0];
rz(-2.0427867) q[0];
sx q[0];
rz(-0.95979877) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4270664) q[2];
sx q[2];
rz(-0.21103141) q[2];
sx q[2];
rz(-2.2843366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51242764) q[1];
sx q[1];
rz(-1.2841793) q[1];
sx q[1];
rz(-0.54835876) q[1];
rz(1.1744099) q[3];
sx q[3];
rz(-0.76226888) q[3];
sx q[3];
rz(-1.9844733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3441299) q[2];
sx q[2];
rz(-0.8351438) q[2];
sx q[2];
rz(-0.75378913) q[2];
rz(1.3695184) q[3];
sx q[3];
rz(-1.6794208) q[3];
sx q[3];
rz(-0.65521017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9264483) q[0];
sx q[0];
rz(-1.0855874) q[0];
sx q[0];
rz(-1.7468859) q[0];
rz(1.1766379) q[1];
sx q[1];
rz(-1.6329012) q[1];
sx q[1];
rz(0.36697695) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94818753) q[0];
sx q[0];
rz(-0.087554878) q[0];
sx q[0];
rz(2.2746536) q[0];
rz(-pi) q[1];
rz(2.4796333) q[2];
sx q[2];
rz(-0.41482224) q[2];
sx q[2];
rz(-0.38554672) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6527892) q[1];
sx q[1];
rz(-1.7646043) q[1];
sx q[1];
rz(-0.36313063) q[1];
x q[2];
rz(0.36716299) q[3];
sx q[3];
rz(-1.490574) q[3];
sx q[3];
rz(1.0368766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41161141) q[2];
sx q[2];
rz(-1.9133762) q[2];
sx q[2];
rz(1.270594) q[2];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56080317) q[0];
sx q[0];
rz(-2.2152948) q[0];
sx q[0];
rz(-1.8286937) q[0];
rz(1.987223) q[1];
sx q[1];
rz(-1.9998974) q[1];
sx q[1];
rz(0.55783522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8759039) q[0];
sx q[0];
rz(-0.77123986) q[0];
sx q[0];
rz(2.3452167) q[0];
rz(2.6959582) q[2];
sx q[2];
rz(-2.9639396) q[2];
sx q[2];
rz(-0.38933094) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5231042) q[1];
sx q[1];
rz(-2.2918211) q[1];
sx q[1];
rz(-2.7343661) q[1];
x q[2];
rz(1.4194896) q[3];
sx q[3];
rz(-2.0278721) q[3];
sx q[3];
rz(-3.0416656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.64342) q[2];
sx q[2];
rz(-1.2926481) q[2];
sx q[2];
rz(-2.2922929) q[2];
rz(-2.9863206) q[3];
sx q[3];
rz(-0.97584358) q[3];
sx q[3];
rz(0.37080216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3575386) q[0];
sx q[0];
rz(-2.5548866) q[0];
sx q[0];
rz(-2.6724755) q[0];
rz(0.47438374) q[1];
sx q[1];
rz(-2.3553039) q[1];
sx q[1];
rz(2.3862086) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7256713) q[0];
sx q[0];
rz(-1.5752107) q[0];
sx q[0];
rz(-1.2995385) q[0];
x q[1];
rz(2.5126476) q[2];
sx q[2];
rz(-1.2997932) q[2];
sx q[2];
rz(2.6586201) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5522105) q[1];
sx q[1];
rz(-0.24866074) q[1];
sx q[1];
rz(2.1099439) q[1];
x q[2];
rz(1.1047885) q[3];
sx q[3];
rz(-0.7940404) q[3];
sx q[3];
rz(3.1178302) q[3];
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
rz(-2.0511621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.93194) q[0];
sx q[0];
rz(-1.8745475) q[0];
sx q[0];
rz(-0.39091045) q[0];
rz(-2.2115425) q[1];
sx q[1];
rz(-1.4314194) q[1];
sx q[1];
rz(1.3040868) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6025897) q[0];
sx q[0];
rz(-2.7900759) q[0];
sx q[0];
rz(1.6319265) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.485059) q[2];
sx q[2];
rz(-1.2303599) q[2];
sx q[2];
rz(-0.96735937) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.75445777) q[1];
sx q[1];
rz(-1.7384733) q[1];
sx q[1];
rz(-1.1134336) q[1];
rz(-pi) q[2];
rz(2.845302) q[3];
sx q[3];
rz(-1.5888831) q[3];
sx q[3];
rz(-2.2294105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8730674) q[2];
sx q[2];
rz(-0.26669845) q[2];
sx q[2];
rz(3.0282057) q[2];
rz(-0.015965613) q[3];
sx q[3];
rz(-1.4957875) q[3];
sx q[3];
rz(2.9769843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3575344) q[0];
sx q[0];
rz(-1.4736195) q[0];
sx q[0];
rz(-3.0175324) q[0];
rz(0.1217753) q[1];
sx q[1];
rz(-2.3901794) q[1];
sx q[1];
rz(1.442499) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6704491) q[0];
sx q[0];
rz(-1.6713977) q[0];
sx q[0];
rz(-1.1858245) q[0];
rz(2.997918) q[2];
sx q[2];
rz(-2.8468003) q[2];
sx q[2];
rz(-0.93570342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.21181078) q[1];
sx q[1];
rz(-0.36961296) q[1];
sx q[1];
rz(-1.1951642) q[1];
rz(-pi) q[2];
rz(0.27713953) q[3];
sx q[3];
rz(-0.74944121) q[3];
sx q[3];
rz(-0.14271862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1759935) q[2];
sx q[2];
rz(-1.3600574) q[2];
sx q[2];
rz(2.3967801) q[2];
rz(2.2143769) q[3];
sx q[3];
rz(-1.6330481) q[3];
sx q[3];
rz(2.0492699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6532779) q[0];
sx q[0];
rz(-1.7019685) q[0];
sx q[0];
rz(0.22536817) q[0];
rz(-1.1124181) q[1];
sx q[1];
rz(-0.77443361) q[1];
sx q[1];
rz(2.2427799) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32466896) q[0];
sx q[0];
rz(-2.3153439) q[0];
sx q[0];
rz(0.81314317) q[0];
rz(-pi) q[1];
rz(-1.9832916) q[2];
sx q[2];
rz(-1.0935244) q[2];
sx q[2];
rz(-1.080918) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73279335) q[1];
sx q[1];
rz(-2.2299754) q[1];
sx q[1];
rz(-0.2562457) q[1];
rz(-pi) q[2];
rz(1.6055272) q[3];
sx q[3];
rz(-2.8485549) q[3];
sx q[3];
rz(1.7188096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0972458) q[2];
sx q[2];
rz(-1.4298507) q[2];
sx q[2];
rz(0.35935768) q[2];
rz(0.25092956) q[3];
sx q[3];
rz(-2.0693306) q[3];
sx q[3];
rz(-2.3738764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3971685) q[0];
sx q[0];
rz(-3.011062) q[0];
sx q[0];
rz(2.2644444) q[0];
rz(-1.5769222) q[1];
sx q[1];
rz(-1.501187) q[1];
sx q[1];
rz(-2.3559779) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6495384) q[0];
sx q[0];
rz(-1.3527332) q[0];
sx q[0];
rz(-1.3278924) q[0];
x q[1];
rz(-0.53391407) q[2];
sx q[2];
rz(-1.1964239) q[2];
sx q[2];
rz(1.716734) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.091508063) q[1];
sx q[1];
rz(-2.7899402) q[1];
sx q[1];
rz(-0.20594724) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7212058) q[3];
sx q[3];
rz(-2.0724943) q[3];
sx q[3];
rz(-3.0034163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.91763532) q[2];
sx q[2];
rz(-0.79006299) q[2];
sx q[2];
rz(-2.4033974) q[2];
rz(2.2186642) q[3];
sx q[3];
rz(-1.3083369) q[3];
sx q[3];
rz(2.8452528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.34687635) q[0];
sx q[0];
rz(-1.3608169) q[0];
sx q[0];
rz(0.97186744) q[0];
rz(-0.90732668) q[1];
sx q[1];
rz(-2.0569888) q[1];
sx q[1];
rz(2.8203698) q[1];
rz(0.99030607) q[2];
sx q[2];
rz(-1.9038426) q[2];
sx q[2];
rz(2.7833084) q[2];
rz(-2.3774556) q[3];
sx q[3];
rz(-0.69307477) q[3];
sx q[3];
rz(0.57144036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
