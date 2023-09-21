OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3333617) q[0];
sx q[0];
rz(-0.97167492) q[0];
sx q[0];
rz(1.6605641) q[0];
rz(3.0193168) q[1];
sx q[1];
rz(-3.0552157) q[1];
sx q[1];
rz(3.1187305) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6883144) q[0];
sx q[0];
rz(-2.9524321) q[0];
sx q[0];
rz(-2.2384089) q[0];
x q[1];
rz(0.68732287) q[2];
sx q[2];
rz(-2.5097373) q[2];
sx q[2];
rz(-0.99422115) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3373128) q[1];
sx q[1];
rz(-0.37706456) q[1];
sx q[1];
rz(-1.1537329) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1894737) q[3];
sx q[3];
rz(-1.7784356) q[3];
sx q[3];
rz(1.729897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7227398) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(2.375405) q[2];
rz(-2.9700759) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(0.086659327) q[0];
rz(-0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-0.08509732) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3413977) q[0];
sx q[0];
rz(-1.5872103) q[0];
sx q[0];
rz(2.5192758) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30479635) q[2];
sx q[2];
rz(-2.1509503) q[2];
sx q[2];
rz(3.0457029) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47463402) q[1];
sx q[1];
rz(-0.14161319) q[1];
sx q[1];
rz(-0.9160362) q[1];
x q[2];
rz(1.1005136) q[3];
sx q[3];
rz(-1.8575875) q[3];
sx q[3];
rz(-1.4834529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(1.6764199) q[2];
rz(-1.9654467) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(-2.8951077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572606) q[0];
sx q[0];
rz(-0.13467877) q[0];
sx q[0];
rz(0.9056257) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.7571626) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15554409) q[0];
sx q[0];
rz(-1.4882003) q[0];
sx q[0];
rz(0.47604576) q[0];
rz(-pi) q[1];
rz(0.39321123) q[2];
sx q[2];
rz(-1.097743) q[2];
sx q[2];
rz(-0.87393239) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3123734) q[1];
sx q[1];
rz(-2.3260498) q[1];
sx q[1];
rz(-1.6014598) q[1];
rz(-pi) q[2];
rz(2.609385) q[3];
sx q[3];
rz(-2.0720707) q[3];
sx q[3];
rz(1.7621079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.24421346) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(1.3282233) q[2];
rz(-1.7437079) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.06552799) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-2.4705825) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(0.20733325) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57732338) q[0];
sx q[0];
rz(-0.042357001) q[0];
sx q[0];
rz(1.3576515) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3992594) q[2];
sx q[2];
rz(-2.9651387) q[2];
sx q[2];
rz(1.2982969) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85992766) q[1];
sx q[1];
rz(-1.7237701) q[1];
sx q[1];
rz(1.330919) q[1];
rz(-pi) q[2];
rz(1.5563577) q[3];
sx q[3];
rz(-1.0591905) q[3];
sx q[3];
rz(-2.5479941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.31071445) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(-2.3332398) q[2];
rz(-2.3338142) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(-0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59072524) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(0.89548683) q[0];
rz(2.8764309) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(-0.53363824) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69913188) q[0];
sx q[0];
rz(-1.8457544) q[0];
sx q[0];
rz(1.2752227) q[0];
rz(0.76782121) q[2];
sx q[2];
rz(-1.4066833) q[2];
sx q[2];
rz(-2.5196911) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6753908) q[1];
sx q[1];
rz(-1.7168105) q[1];
sx q[1];
rz(0.22202613) q[1];
rz(-pi) q[2];
x q[2];
rz(2.93612) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(-0.56838244) q[2];
rz(-1.9249632) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(-2.7980428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.50826532) q[0];
sx q[0];
rz(-2.2631622) q[0];
sx q[0];
rz(1.1761965) q[0];
rz(-1.5856702) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(-1.0046545) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3096136) q[0];
sx q[0];
rz(-1.2011315) q[0];
sx q[0];
rz(-1.2020822) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19095687) q[2];
sx q[2];
rz(-0.40340427) q[2];
sx q[2];
rz(-2.0463338) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.55014729) q[1];
sx q[1];
rz(-1.0284371) q[1];
sx q[1];
rz(-1.28246) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9846356) q[3];
sx q[3];
rz(-1.5817747) q[3];
sx q[3];
rz(-1.4423086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4040318) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(-2.8435977) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649243) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(-0.65814322) q[0];
rz(1.8572042) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(2.4694494) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9400529) q[0];
sx q[0];
rz(-1.4048028) q[0];
sx q[0];
rz(-3.1340909) q[0];
rz(2.3357046) q[2];
sx q[2];
rz(-0.58008367) q[2];
sx q[2];
rz(0.16575955) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0828515) q[1];
sx q[1];
rz(-0.62742678) q[1];
sx q[1];
rz(-0.65545603) q[1];
rz(-pi) q[2];
rz(2.9716306) q[3];
sx q[3];
rz(-2.4331577) q[3];
sx q[3];
rz(1.8426614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.60823524) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(-0.39880025) q[2];
rz(-2.3170025) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7010715) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(-2.9833802) q[0];
rz(2.054706) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(0.80668443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0772484) q[0];
sx q[0];
rz(-1.6883388) q[0];
sx q[0];
rz(2.006152) q[0];
rz(-pi) q[1];
rz(-2.9358747) q[2];
sx q[2];
rz(-1.6763902) q[2];
sx q[2];
rz(-0.92668698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.17864171) q[1];
sx q[1];
rz(-2.2370403) q[1];
sx q[1];
rz(-1.9402177) q[1];
x q[2];
rz(1.4679457) q[3];
sx q[3];
rz(-2.8392483) q[3];
sx q[3];
rz(2.9000988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.57758254) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(0.99964833) q[2];
rz(3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0209811) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(-0.28655562) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(2.8009169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0176312) q[0];
sx q[0];
rz(-2.258856) q[0];
sx q[0];
rz(2.7244363) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0057783) q[2];
sx q[2];
rz(-1.3823595) q[2];
sx q[2];
rz(2.4609158) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.162519) q[1];
sx q[1];
rz(-0.58774978) q[1];
sx q[1];
rz(-1.9025004) q[1];
rz(-pi) q[2];
rz(1.4021923) q[3];
sx q[3];
rz(-1.3753034) q[3];
sx q[3];
rz(-2.33778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.75491536) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(2.6436451) q[2];
rz(-0.30161101) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.4122445) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(2.8630032) q[0];
rz(0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(3.0648807) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0522703) q[0];
sx q[0];
rz(-0.33903402) q[0];
sx q[0];
rz(-1.1917172) q[0];
rz(-pi) q[1];
rz(-2.9701091) q[2];
sx q[2];
rz(-1.0043317) q[2];
sx q[2];
rz(2.7330074) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.50347603) q[1];
sx q[1];
rz(-1.5498811) q[1];
sx q[1];
rz(-2.7760844) q[1];
rz(-pi) q[2];
rz(-1.65927) q[3];
sx q[3];
rz(-1.7469329) q[3];
sx q[3];
rz(-2.4000771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.77438337) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(2.7101743) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(-0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(-2.5554399) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(-2.8557599) q[2];
sx q[2];
rz(-1.1163859) q[2];
sx q[2];
rz(-3.1039539) q[2];
rz(2.7145731) q[3];
sx q[3];
rz(-1.2357161) q[3];
sx q[3];
rz(-1.2824035) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];