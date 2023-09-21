OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5157226) q[0];
sx q[0];
rz(-0.54870257) q[0];
sx q[0];
rz(-0.8843511) q[0];
rz(-1.7110775) q[1];
sx q[1];
rz(-0.95354748) q[1];
sx q[1];
rz(1.6391099) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2384773) q[0];
sx q[0];
rz(-1.7503386) q[0];
sx q[0];
rz(1.205501) q[0];
x q[1];
rz(3.1332364) q[2];
sx q[2];
rz(-2.5764321) q[2];
sx q[2];
rz(1.1430119) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1391746) q[1];
sx q[1];
rz(-0.25484172) q[1];
sx q[1];
rz(0.26267085) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9804269) q[3];
sx q[3];
rz(-0.75244609) q[3];
sx q[3];
rz(2.9890271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.78645906) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(-0.65594977) q[2];
rz(1.2077228) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(2.1470127) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2475964) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(-2.7080652) q[0];
rz(2.9128089) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(-0.0072335009) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15046209) q[0];
sx q[0];
rz(-2.7379002) q[0];
sx q[0];
rz(-1.3315721) q[0];
rz(-pi) q[1];
rz(2.8480808) q[2];
sx q[2];
rz(-1.268317) q[2];
sx q[2];
rz(-0.40700618) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.05939535) q[1];
sx q[1];
rz(-2.7121183) q[1];
sx q[1];
rz(-2.5897964) q[1];
x q[2];
rz(2.1528483) q[3];
sx q[3];
rz(-2.3855004) q[3];
sx q[3];
rz(-2.3088629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6146415) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(2.4439404) q[2];
rz(-3.0200322) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(0.30383032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6737297) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(1.3695705) q[0];
rz(1.9000152) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(0.26161584) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3221489) q[0];
sx q[0];
rz(-1.561164) q[0];
sx q[0];
rz(1.5887512) q[0];
x q[1];
rz(0.82654731) q[2];
sx q[2];
rz(-1.4790223) q[2];
sx q[2];
rz(-1.9020136) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9237823) q[1];
sx q[1];
rz(-2.3465956) q[1];
sx q[1];
rz(1.7407106) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6286078) q[3];
sx q[3];
rz(-1.5981711) q[3];
sx q[3];
rz(2.5740636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6802784) q[2];
sx q[2];
rz(-1.5051944) q[2];
sx q[2];
rz(-0.20351163) q[2];
rz(-2.2198548) q[3];
sx q[3];
rz(-1.2676055) q[3];
sx q[3];
rz(0.27954277) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97776425) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(1.5699566) q[0];
rz(-2.1381901) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(1.2483695) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0715863) q[0];
sx q[0];
rz(-1.6674111) q[0];
sx q[0];
rz(-0.065652547) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3893045) q[2];
sx q[2];
rz(-1.3647807) q[2];
sx q[2];
rz(0.69603053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8716988) q[1];
sx q[1];
rz(-1.6819685) q[1];
sx q[1];
rz(0.65472366) q[1];
x q[2];
rz(1.3094041) q[3];
sx q[3];
rz(-2.2769615) q[3];
sx q[3];
rz(-0.9725001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0288329) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(-2.3045585) q[2];
rz(1.933243) q[3];
sx q[3];
rz(-1.2669867) q[3];
sx q[3];
rz(-0.78554955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0060624881) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(2.3663882) q[0];
rz(-0.40183055) q[1];
sx q[1];
rz(-2.1907175) q[1];
sx q[1];
rz(0.90243375) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1989312) q[0];
sx q[0];
rz(-0.64093243) q[0];
sx q[0];
rz(-0.076365691) q[0];
rz(0.47307737) q[2];
sx q[2];
rz(-0.50191754) q[2];
sx q[2];
rz(-0.15854533) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.93859766) q[1];
sx q[1];
rz(-2.1122167) q[1];
sx q[1];
rz(-2.2375537) q[1];
rz(-pi) q[2];
rz(0.0028354672) q[3];
sx q[3];
rz(-1.8875202) q[3];
sx q[3];
rz(2.1383274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.19501413) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(-1.1966594) q[2];
rz(-1.6992016) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(1.4590013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8114132) q[0];
sx q[0];
rz(-0.17689642) q[0];
sx q[0];
rz(0.50317558) q[0];
rz(-1.4563837) q[1];
sx q[1];
rz(-1.0738942) q[1];
sx q[1];
rz(-2.9398289) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1267032) q[0];
sx q[0];
rz(-0.12173437) q[0];
sx q[0];
rz(-0.99862167) q[0];
rz(-2.6642338) q[2];
sx q[2];
rz(-1.9743894) q[2];
sx q[2];
rz(-1.8237643) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.7355431) q[1];
sx q[1];
rz(-2.4821401) q[1];
sx q[1];
rz(-2.4156648) q[1];
rz(-pi) q[2];
rz(0.73268907) q[3];
sx q[3];
rz(-1.9273888) q[3];
sx q[3];
rz(0.69571146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.21489828) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(0.26724896) q[2];
rz(-2.3184508) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(0.87583035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.293752) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(-0.057549495) q[0];
rz(-1.4808222) q[1];
sx q[1];
rz(-1.3596423) q[1];
sx q[1];
rz(-2.1988791) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4276886) q[0];
sx q[0];
rz(-1.9644992) q[0];
sx q[0];
rz(0.64481553) q[0];
rz(0.96037453) q[2];
sx q[2];
rz(-2.865961) q[2];
sx q[2];
rz(0.20197091) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6616933) q[1];
sx q[1];
rz(-0.53240314) q[1];
sx q[1];
rz(-0.9049306) q[1];
x q[2];
rz(2.8816678) q[3];
sx q[3];
rz(-1.0726895) q[3];
sx q[3];
rz(-0.92344027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1192347) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(-0.38267246) q[2];
rz(2.102397) q[3];
sx q[3];
rz(-1.305205) q[3];
sx q[3];
rz(-2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7169749) q[0];
sx q[0];
rz(-0.096352339) q[0];
sx q[0];
rz(2.8714645) q[0];
rz(-2.5121636) q[1];
sx q[1];
rz(-2.4286178) q[1];
sx q[1];
rz(-0.28392917) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89284183) q[0];
sx q[0];
rz(-0.98309702) q[0];
sx q[0];
rz(2.6177004) q[0];
x q[1];
rz(-0.59992744) q[2];
sx q[2];
rz(-1.6919961) q[2];
sx q[2];
rz(1.3231414) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4979907) q[1];
sx q[1];
rz(-1.0712578) q[1];
sx q[1];
rz(-1.040578) q[1];
rz(-pi) q[2];
rz(3.1030032) q[3];
sx q[3];
rz(-0.66017294) q[3];
sx q[3];
rz(0.71036464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.55390629) q[2];
sx q[2];
rz(-2.9710785) q[2];
sx q[2];
rz(-1.930687) q[2];
rz(0.25990137) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07638409) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(-2.912345) q[0];
rz(-0.30300888) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(-1.680826) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54366771) q[0];
sx q[0];
rz(-0.28408465) q[0];
sx q[0];
rz(3.0790867) q[0];
x q[1];
rz(-0.37947189) q[2];
sx q[2];
rz(-1.4037637) q[2];
sx q[2];
rz(2.6118979) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8742121) q[1];
sx q[1];
rz(-2.0955718) q[1];
sx q[1];
rz(-1.5966148) q[1];
rz(1.9577515) q[3];
sx q[3];
rz(-2.317252) q[3];
sx q[3];
rz(-1.4617621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4298657) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(2.4460068) q[2];
rz(-2.7097278) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(2.4263884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5678976) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(2.8046872) q[0];
rz(-0.20740549) q[1];
sx q[1];
rz(-2.1131056) q[1];
sx q[1];
rz(-2.7609603) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99285179) q[0];
sx q[0];
rz(-1.7721575) q[0];
sx q[0];
rz(-3.0701748) q[0];
x q[1];
rz(2.5194174) q[2];
sx q[2];
rz(-0.84465775) q[2];
sx q[2];
rz(-0.13571339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.759728) q[1];
sx q[1];
rz(-1.3761531) q[1];
sx q[1];
rz(1.0641644) q[1];
rz(2.5467039) q[3];
sx q[3];
rz(-2.4747362) q[3];
sx q[3];
rz(2.8951333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1404861) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(-2.005119) q[2];
rz(0.040955695) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(-1.3142746) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3363591) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(2.1144755) q[1];
sx q[1];
rz(-1.2925016) q[1];
sx q[1];
rz(2.1137994) q[1];
rz(2.3958191) q[2];
sx q[2];
rz(-2.8095828) q[2];
sx q[2];
rz(-2.7381894) q[2];
rz(-0.20053486) q[3];
sx q[3];
rz(-1.1391098) q[3];
sx q[3];
rz(-1.1548635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];