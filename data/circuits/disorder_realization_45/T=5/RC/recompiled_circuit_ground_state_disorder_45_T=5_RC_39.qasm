OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1096119) q[0];
sx q[0];
rz(-1.0053827) q[0];
sx q[0];
rz(-2.1729852) q[0];
rz(-2.3321416) q[1];
sx q[1];
rz(-1.6074155) q[1];
sx q[1];
rz(0.029835116) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3650509) q[0];
sx q[0];
rz(-1.33857) q[0];
sx q[0];
rz(-2.5864629) q[0];
x q[1];
rz(-1.6462684) q[2];
sx q[2];
rz(-0.75368893) q[2];
sx q[2];
rz(1.8729295) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0152143) q[1];
sx q[1];
rz(-0.83288899) q[1];
sx q[1];
rz(2.6070541) q[1];
x q[2];
rz(-2.5621595) q[3];
sx q[3];
rz(-2.979485) q[3];
sx q[3];
rz(2.1539089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6385103) q[2];
sx q[2];
rz(-2.5677887) q[2];
sx q[2];
rz(-0.32162515) q[2];
rz(-2.5173748) q[3];
sx q[3];
rz(-1.7523242) q[3];
sx q[3];
rz(1.483884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44840789) q[0];
sx q[0];
rz(-0.25633651) q[0];
sx q[0];
rz(-2.2053027) q[0];
rz(1.869465) q[1];
sx q[1];
rz(-1.5969758) q[1];
sx q[1];
rz(1.3880656) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0569071) q[0];
sx q[0];
rz(-1.7900247) q[0];
sx q[0];
rz(0.058554336) q[0];
rz(-0.53728062) q[2];
sx q[2];
rz(-0.85308077) q[2];
sx q[2];
rz(-3.0514015) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.1759049) q[1];
sx q[1];
rz(-0.85077326) q[1];
sx q[1];
rz(-0.73890025) q[1];
x q[2];
rz(1.8348946) q[3];
sx q[3];
rz(-1.3624884) q[3];
sx q[3];
rz(-1.1549319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0973413) q[2];
sx q[2];
rz(-0.96434957) q[2];
sx q[2];
rz(-2.468289) q[2];
rz(1.8246957) q[3];
sx q[3];
rz(-1.293106) q[3];
sx q[3];
rz(0.84002686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3299385) q[0];
sx q[0];
rz(-1.543777) q[0];
sx q[0];
rz(0.04976186) q[0];
rz(0.20637575) q[1];
sx q[1];
rz(-1.3866837) q[1];
sx q[1];
rz(1.4899563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7508036) q[0];
sx q[0];
rz(-0.62701462) q[0];
sx q[0];
rz(0.20558447) q[0];
x q[1];
rz(-1.1300722) q[2];
sx q[2];
rz(-1.0606546) q[2];
sx q[2];
rz(-2.5464692) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4088285) q[1];
sx q[1];
rz(-1.2868501) q[1];
sx q[1];
rz(-0.27393053) q[1];
rz(0.95155119) q[3];
sx q[3];
rz(-2.5971322) q[3];
sx q[3];
rz(1.8306554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9516307) q[2];
sx q[2];
rz(-2.5962574) q[2];
sx q[2];
rz(3.087431) q[2];
rz(1.1545898) q[3];
sx q[3];
rz(-1.2879939) q[3];
sx q[3];
rz(1.1156999) q[3];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40256777) q[0];
sx q[0];
rz(-1.716528) q[0];
sx q[0];
rz(2.2836852) q[0];
rz(0.67445451) q[1];
sx q[1];
rz(-2.1916788) q[1];
sx q[1];
rz(1.3309006) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9337263) q[0];
sx q[0];
rz(-0.27369341) q[0];
sx q[0];
rz(2.6143608) q[0];
rz(-pi) q[1];
rz(-1.2670008) q[2];
sx q[2];
rz(-1.7335606) q[2];
sx q[2];
rz(-0.13653423) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7965296) q[1];
sx q[1];
rz(-1.3994675) q[1];
sx q[1];
rz(-1.0061128) q[1];
rz(-1.5783806) q[3];
sx q[3];
rz(-2.7223153) q[3];
sx q[3];
rz(-1.5602001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3746419) q[2];
sx q[2];
rz(-2.9810814) q[2];
sx q[2];
rz(-0.10425076) q[2];
rz(0.57779622) q[3];
sx q[3];
rz(-1.0775074) q[3];
sx q[3];
rz(-0.92208636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1330971) q[0];
sx q[0];
rz(-1.6660026) q[0];
sx q[0];
rz(2.1971028) q[0];
rz(1.7780444) q[1];
sx q[1];
rz(-1.8133769) q[1];
sx q[1];
rz(1.4469226) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1773374) q[0];
sx q[0];
rz(-0.87133555) q[0];
sx q[0];
rz(0.2486458) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7124921) q[2];
sx q[2];
rz(-2.8887914) q[2];
sx q[2];
rz(1.6293874) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8268298) q[1];
sx q[1];
rz(-0.75195014) q[1];
sx q[1];
rz(0.58677267) q[1];
rz(-0.99926104) q[3];
sx q[3];
rz(-2.0700698) q[3];
sx q[3];
rz(1.5972135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.038534433) q[2];
sx q[2];
rz(-2.6699622) q[2];
sx q[2];
rz(-0.72059694) q[2];
rz(3.0125812) q[3];
sx q[3];
rz(-1.7804264) q[3];
sx q[3];
rz(-1.2989929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1593889) q[0];
sx q[0];
rz(-2.144618) q[0];
sx q[0];
rz(2.5250315) q[0];
rz(-2.189134) q[1];
sx q[1];
rz(-1.9472232) q[1];
sx q[1];
rz(1.8985101) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0117802) q[0];
sx q[0];
rz(-1.3115378) q[0];
sx q[0];
rz(-3.0002322) q[0];
rz(-2.0427734) q[2];
sx q[2];
rz(-2.3772559) q[2];
sx q[2];
rz(2.7175222) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3042708) q[1];
sx q[1];
rz(-0.72266403) q[1];
sx q[1];
rz(-0.091462084) q[1];
rz(-pi) q[2];
rz(0.9493306) q[3];
sx q[3];
rz(-1.5104745) q[3];
sx q[3];
rz(0.5746791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9055966) q[2];
sx q[2];
rz(-2.8541028) q[2];
sx q[2];
rz(3.0671885) q[2];
rz(-2.060804) q[3];
sx q[3];
rz(-1.6483043) q[3];
sx q[3];
rz(1.0409482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3022795) q[0];
sx q[0];
rz(-1.9957207) q[0];
sx q[0];
rz(-3.0795414) q[0];
rz(1.7265559) q[1];
sx q[1];
rz(-2.3686385) q[1];
sx q[1];
rz(-0.089990377) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0768698) q[0];
sx q[0];
rz(-1.8241811) q[0];
sx q[0];
rz(-0.86034448) q[0];
rz(-2.2587218) q[2];
sx q[2];
rz(-1.6992879) q[2];
sx q[2];
rz(1.524432) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6288694) q[1];
sx q[1];
rz(-0.97755331) q[1];
sx q[1];
rz(-2.7166461) q[1];
x q[2];
rz(2.2136647) q[3];
sx q[3];
rz(-1.0786329) q[3];
sx q[3];
rz(-1.6097691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.0061079582) q[2];
sx q[2];
rz(-1.8718448) q[2];
sx q[2];
rz(-0.6832068) q[2];
rz(-1.3225383) q[3];
sx q[3];
rz(-0.95219487) q[3];
sx q[3];
rz(-1.7542084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4676062) q[0];
sx q[0];
rz(-1.645393) q[0];
sx q[0];
rz(0.4678539) q[0];
rz(-0.6894919) q[1];
sx q[1];
rz(-2.4263224) q[1];
sx q[1];
rz(2.668344) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5734542) q[0];
sx q[0];
rz(-1.3949419) q[0];
sx q[0];
rz(-3.0699845) q[0];
rz(-2.6976349) q[2];
sx q[2];
rz(-1.7164162) q[2];
sx q[2];
rz(-0.99010003) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6618408) q[1];
sx q[1];
rz(-1.5532232) q[1];
sx q[1];
rz(0.015546202) q[1];
x q[2];
rz(-0.28041081) q[3];
sx q[3];
rz(-1.8160928) q[3];
sx q[3];
rz(-2.6287358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1271992) q[2];
sx q[2];
rz(-1.8202123) q[2];
sx q[2];
rz(-0.42352208) q[2];
rz(-1.315377) q[3];
sx q[3];
rz(-2.6500621) q[3];
sx q[3];
rz(1.2929085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92726436) q[0];
sx q[0];
rz(-0.94896999) q[0];
sx q[0];
rz(1.3342791) q[0];
rz(-1.3182053) q[1];
sx q[1];
rz(-0.91000906) q[1];
sx q[1];
rz(-0.013890161) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2023808) q[0];
sx q[0];
rz(-1.9343543) q[0];
sx q[0];
rz(-0.027287539) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7282762) q[2];
sx q[2];
rz(-1.9183927) q[2];
sx q[2];
rz(1.9633788) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.397032) q[1];
sx q[1];
rz(-2.5066527) q[1];
sx q[1];
rz(1.145361) q[1];
rz(0.12356476) q[3];
sx q[3];
rz(-0.84513226) q[3];
sx q[3];
rz(-0.35620505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5615329) q[2];
sx q[2];
rz(-1.9999028) q[2];
sx q[2];
rz(-2.0884936) q[2];
rz(-0.92755353) q[3];
sx q[3];
rz(-1.2404697) q[3];
sx q[3];
rz(-0.48228669) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7841566) q[0];
sx q[0];
rz(-2.744839) q[0];
sx q[0];
rz(-2.5411153) q[0];
rz(-0.26957574) q[1];
sx q[1];
rz(-0.69200626) q[1];
sx q[1];
rz(1.762766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3665837) q[0];
sx q[0];
rz(-1.3008769) q[0];
sx q[0];
rz(-1.2054514) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7711033) q[2];
sx q[2];
rz(-1.2844119) q[2];
sx q[2];
rz(0.011025393) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3311171) q[1];
sx q[1];
rz(-1.1441191) q[1];
sx q[1];
rz(1.7716549) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61609488) q[3];
sx q[3];
rz(-1.2047676) q[3];
sx q[3];
rz(0.6835365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1158585) q[2];
sx q[2];
rz(-2.8074042) q[2];
sx q[2];
rz(2.1791229) q[2];
rz(0.95364237) q[3];
sx q[3];
rz(-1.0735268) q[3];
sx q[3];
rz(-1.7474705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7643323) q[0];
sx q[0];
rz(-2.2567516) q[0];
sx q[0];
rz(1.2010038) q[0];
rz(1.7878905) q[1];
sx q[1];
rz(-1.1693015) q[1];
sx q[1];
rz(2.404626) q[1];
rz(2.216966) q[2];
sx q[2];
rz(-2.1344816) q[2];
sx q[2];
rz(-0.62350694) q[2];
rz(0.50696452) q[3];
sx q[3];
rz(-1.6929814) q[3];
sx q[3];
rz(-0.89043989) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
