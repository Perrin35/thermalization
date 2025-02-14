OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7656443) q[0];
sx q[0];
rz(-0.41416895) q[0];
sx q[0];
rz(-0.8015269) q[0];
rz(-0.0030567788) q[1];
sx q[1];
rz(-0.86441511) q[1];
sx q[1];
rz(-0.094223082) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.191358) q[0];
sx q[0];
rz(-1.9468029) q[0];
sx q[0];
rz(-2.3467499) q[0];
rz(-pi) q[1];
rz(-1.4613737) q[2];
sx q[2];
rz(-2.8673807) q[2];
sx q[2];
rz(0.69137979) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7677418) q[1];
sx q[1];
rz(-2.5818425) q[1];
sx q[1];
rz(2.1327726) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5685969) q[3];
sx q[3];
rz(-1.723395) q[3];
sx q[3];
rz(-0.27657498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6039383) q[2];
sx q[2];
rz(-0.64017576) q[2];
sx q[2];
rz(1.5113277) q[2];
rz(-0.18621914) q[3];
sx q[3];
rz(-0.43034601) q[3];
sx q[3];
rz(-2.490624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52361012) q[0];
sx q[0];
rz(-1.1815434) q[0];
sx q[0];
rz(3.004916) q[0];
rz(0.094820529) q[1];
sx q[1];
rz(-0.46351981) q[1];
sx q[1];
rz(-2.3725841) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13327293) q[0];
sx q[0];
rz(-1.1964487) q[0];
sx q[0];
rz(1.4248614) q[0];
x q[1];
rz(0.42119512) q[2];
sx q[2];
rz(-2.4026818) q[2];
sx q[2];
rz(0.74904672) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31615801) q[1];
sx q[1];
rz(-2.3806751) q[1];
sx q[1];
rz(-0.23657152) q[1];
rz(-1.8212998) q[3];
sx q[3];
rz(-0.13938306) q[3];
sx q[3];
rz(-1.2047307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7424221) q[2];
sx q[2];
rz(-0.48226446) q[2];
sx q[2];
rz(1.7102309) q[2];
rz(-2.6509905) q[3];
sx q[3];
rz(-0.49121818) q[3];
sx q[3];
rz(-0.77559364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2630149) q[0];
sx q[0];
rz(-0.37913015) q[0];
sx q[0];
rz(-1.0947134) q[0];
rz(-1.8393983) q[1];
sx q[1];
rz(-2.1523988) q[1];
sx q[1];
rz(3.1375695) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1510295) q[0];
sx q[0];
rz(-0.73712611) q[0];
sx q[0];
rz(1.766979) q[0];
rz(2.3962069) q[2];
sx q[2];
rz(-0.89085397) q[2];
sx q[2];
rz(-2.5542575) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6371256) q[1];
sx q[1];
rz(-1.6011304) q[1];
sx q[1];
rz(2.8843969) q[1];
rz(-0.71378848) q[3];
sx q[3];
rz(-1.3828438) q[3];
sx q[3];
rz(-2.1655708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5912938) q[2];
sx q[2];
rz(-0.58612263) q[2];
sx q[2];
rz(-0.44814056) q[2];
rz(-0.48629931) q[3];
sx q[3];
rz(-2.2794162) q[3];
sx q[3];
rz(2.015131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.22911856) q[0];
sx q[0];
rz(-1.4294701) q[0];
sx q[0];
rz(2.4718156) q[0];
rz(1.9122596) q[1];
sx q[1];
rz(-0.45893097) q[1];
sx q[1];
rz(-2.5965447) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6908643) q[0];
sx q[0];
rz(-2.0929027) q[0];
sx q[0];
rz(-2.1073116) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5169536) q[2];
sx q[2];
rz(-2.1339941) q[2];
sx q[2];
rz(0.27499946) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4008444) q[1];
sx q[1];
rz(-1.5103589) q[1];
sx q[1];
rz(1.7857741) q[1];
rz(-pi) q[2];
rz(-2.3432076) q[3];
sx q[3];
rz(-1.9163864) q[3];
sx q[3];
rz(1.2452728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5547319) q[2];
sx q[2];
rz(-1.6394337) q[2];
sx q[2];
rz(-2.9736605) q[2];
rz(1.933291) q[3];
sx q[3];
rz(-2.8390563) q[3];
sx q[3];
rz(-1.4471853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.55030441) q[0];
sx q[0];
rz(-1.8809603) q[0];
sx q[0];
rz(3.1377129) q[0];
rz(2.8344391) q[1];
sx q[1];
rz(-2.3852564) q[1];
sx q[1];
rz(-0.53877962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2505441) q[0];
sx q[0];
rz(-1.7122503) q[0];
sx q[0];
rz(0.23141872) q[0];
x q[1];
rz(-1.7096552) q[2];
sx q[2];
rz(-1.9966239) q[2];
sx q[2];
rz(-0.52639293) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51018029) q[1];
sx q[1];
rz(-2.4459527) q[1];
sx q[1];
rz(-2.5744777) q[1];
rz(-pi) q[2];
x q[2];
rz(2.410048) q[3];
sx q[3];
rz(-1.2243825) q[3];
sx q[3];
rz(-0.15958318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6946081) q[2];
sx q[2];
rz(-1.0993404) q[2];
sx q[2];
rz(1.851932) q[2];
rz(-1.5223711) q[3];
sx q[3];
rz(-0.74519849) q[3];
sx q[3];
rz(-1.6600018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30037844) q[0];
sx q[0];
rz(-2.2024246) q[0];
sx q[0];
rz(0.10264957) q[0];
rz(2.4094021) q[1];
sx q[1];
rz(-2.3468572) q[1];
sx q[1];
rz(0.57714677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6711376) q[0];
sx q[0];
rz(-1.5744493) q[0];
sx q[0];
rz(-1.5994344) q[0];
rz(-pi) q[1];
rz(2.5546165) q[2];
sx q[2];
rz(-1.7086626) q[2];
sx q[2];
rz(1.1816813) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53458339) q[1];
sx q[1];
rz(-1.680169) q[1];
sx q[1];
rz(1.6499728) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34586819) q[3];
sx q[3];
rz(-0.78639275) q[3];
sx q[3];
rz(2.6707471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.010043667) q[2];
sx q[2];
rz(-2.4250344) q[2];
sx q[2];
rz(-1.5232167) q[2];
rz(-0.067954436) q[3];
sx q[3];
rz(-2.8165635) q[3];
sx q[3];
rz(2.3006191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1286569) q[0];
sx q[0];
rz(-1.034863) q[0];
sx q[0];
rz(2.3701684) q[0];
rz(0.5386638) q[1];
sx q[1];
rz(-1.3464876) q[1];
sx q[1];
rz(-2.0753986) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5504042) q[0];
sx q[0];
rz(-0.26445358) q[0];
sx q[0];
rz(1.868684) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0838924) q[2];
sx q[2];
rz(-1.7745681) q[2];
sx q[2];
rz(-2.9625826) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33246189) q[1];
sx q[1];
rz(-1.5839424) q[1];
sx q[1];
rz(2.9244208) q[1];
x q[2];
rz(-1.3049502) q[3];
sx q[3];
rz(-1.2472594) q[3];
sx q[3];
rz(2.8697447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.44925877) q[2];
sx q[2];
rz(-2.3350495) q[2];
sx q[2];
rz(0.44245693) q[2];
rz(2.9509406) q[3];
sx q[3];
rz(-2.4176044) q[3];
sx q[3];
rz(-0.75907069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.2451179) q[0];
sx q[0];
rz(-2.9560095) q[0];
sx q[0];
rz(1.3099439) q[0];
rz(-0.39778057) q[1];
sx q[1];
rz(-2.3186627) q[1];
sx q[1];
rz(1.3936874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8373972) q[0];
sx q[0];
rz(-2.1050354) q[0];
sx q[0];
rz(-0.5128575) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91266914) q[2];
sx q[2];
rz(-2.7670112) q[2];
sx q[2];
rz(-2.7233469) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.950395) q[1];
sx q[1];
rz(-1.1069595) q[1];
sx q[1];
rz(1.2359124) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10618629) q[3];
sx q[3];
rz(-1.7683889) q[3];
sx q[3];
rz(1.1923517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.028367793) q[2];
sx q[2];
rz(-0.71030474) q[2];
sx q[2];
rz(3.0681211) q[2];
rz(2.4469589) q[3];
sx q[3];
rz(-1.2725384) q[3];
sx q[3];
rz(1.0276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6756814) q[0];
sx q[0];
rz(-2.9407192) q[0];
sx q[0];
rz(2.1821816) q[0];
rz(2.361946) q[1];
sx q[1];
rz(-2.1387073) q[1];
sx q[1];
rz(-0.27228212) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0060385614) q[0];
sx q[0];
rz(-1.5840496) q[0];
sx q[0];
rz(3.1363961) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.831976) q[2];
sx q[2];
rz(-2.2520492) q[2];
sx q[2];
rz(-0.75971425) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.694937) q[1];
sx q[1];
rz(-2.4866101) q[1];
sx q[1];
rz(-0.86465624) q[1];
rz(-pi) q[2];
rz(0.63377728) q[3];
sx q[3];
rz(-2.2867775) q[3];
sx q[3];
rz(0.26308003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4043364) q[2];
sx q[2];
rz(-1.5055089) q[2];
sx q[2];
rz(2.9340202) q[2];
rz(-0.10457822) q[3];
sx q[3];
rz(-0.50555491) q[3];
sx q[3];
rz(-1.526621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.680147) q[0];
sx q[0];
rz(-1.2274281) q[0];
sx q[0];
rz(1.0567868) q[0];
rz(2.5959004) q[1];
sx q[1];
rz(-0.97815198) q[1];
sx q[1];
rz(2.812885) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0595039) q[0];
sx q[0];
rz(-0.57841136) q[0];
sx q[0];
rz(0.93579458) q[0];
x q[1];
rz(-0.64735554) q[2];
sx q[2];
rz(-0.090687625) q[2];
sx q[2];
rz(-0.91435223) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4663965) q[1];
sx q[1];
rz(-1.5541142) q[1];
sx q[1];
rz(1.9128602) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77158934) q[3];
sx q[3];
rz(-1.3719012) q[3];
sx q[3];
rz(-1.0071013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1591961) q[2];
sx q[2];
rz(-0.13267645) q[2];
sx q[2];
rz(-1.1245493) q[2];
rz(-3.0647965) q[3];
sx q[3];
rz(-2.7087637) q[3];
sx q[3];
rz(-0.92397773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15976739) q[0];
sx q[0];
rz(-1.8851017) q[0];
sx q[0];
rz(-1.5633352) q[0];
rz(0.81623296) q[1];
sx q[1];
rz(-1.4030133) q[1];
sx q[1];
rz(-1.0907008) q[1];
rz(-0.74538415) q[2];
sx q[2];
rz(-1.1911662) q[2];
sx q[2];
rz(0.62698812) q[2];
rz(0.0088469395) q[3];
sx q[3];
rz(-0.4053517) q[3];
sx q[3];
rz(2.3978618) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
