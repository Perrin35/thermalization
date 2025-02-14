OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47867632) q[0];
sx q[0];
rz(-1.588568) q[0];
sx q[0];
rz(1.8832062) q[0];
rz(1.2619184) q[1];
sx q[1];
rz(-2.6231397) q[1];
sx q[1];
rz(3.0814896) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4661048) q[0];
sx q[0];
rz(-1.6610613) q[0];
sx q[0];
rz(-1.6016927) q[0];
rz(1.7243027) q[2];
sx q[2];
rz(-1.0833086) q[2];
sx q[2];
rz(1.9755028) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2025238) q[1];
sx q[1];
rz(-1.9323404) q[1];
sx q[1];
rz(-0.099379813) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93688776) q[3];
sx q[3];
rz(-2.4408615) q[3];
sx q[3];
rz(2.7752293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1897159) q[2];
sx q[2];
rz(-2.3592301) q[2];
sx q[2];
rz(0.18048364) q[2];
rz(-0.30432025) q[3];
sx q[3];
rz(-2.2173827) q[3];
sx q[3];
rz(-2.9144104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3825398) q[0];
sx q[0];
rz(-0.57342425) q[0];
sx q[0];
rz(0.10391129) q[0];
rz(2.3567764) q[1];
sx q[1];
rz(-1.8321313) q[1];
sx q[1];
rz(0.58194247) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4020445) q[0];
sx q[0];
rz(-0.53240896) q[0];
sx q[0];
rz(3.0976712) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58217073) q[2];
sx q[2];
rz(-2.5522531) q[2];
sx q[2];
rz(-1.2784399) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.90408976) q[1];
sx q[1];
rz(-2.6640035) q[1];
sx q[1];
rz(1.5144004) q[1];
x q[2];
rz(1.9690912) q[3];
sx q[3];
rz(-1.952926) q[3];
sx q[3];
rz(-2.7492461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.71501422) q[2];
sx q[2];
rz(-0.38807401) q[2];
sx q[2];
rz(2.1132054) q[2];
rz(-0.55108023) q[3];
sx q[3];
rz(-1.9469399) q[3];
sx q[3];
rz(-2.6330131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55117637) q[0];
sx q[0];
rz(-1.6352147) q[0];
sx q[0];
rz(-1.2580385) q[0];
rz(-2.538077) q[1];
sx q[1];
rz(-1.3110833) q[1];
sx q[1];
rz(-2.9579732) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98417898) q[0];
sx q[0];
rz(-1.8350198) q[0];
sx q[0];
rz(-1.1452894) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48086353) q[2];
sx q[2];
rz(-1.4521952) q[2];
sx q[2];
rz(3.1080217) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10079846) q[1];
sx q[1];
rz(-1.14193) q[1];
sx q[1];
rz(-2.3593192) q[1];
x q[2];
rz(1.5294936) q[3];
sx q[3];
rz(-1.777711) q[3];
sx q[3];
rz(0.69118354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52643481) q[2];
sx q[2];
rz(-1.7650471) q[2];
sx q[2];
rz(-2.5395425) q[2];
rz(-0.43271068) q[3];
sx q[3];
rz(-1.5475464) q[3];
sx q[3];
rz(1.4759147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7567619) q[0];
sx q[0];
rz(-2.6669406) q[0];
sx q[0];
rz(-2.5153644) q[0];
rz(-1.2681883) q[1];
sx q[1];
rz(-1.405193) q[1];
sx q[1];
rz(2.7558806) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.884366) q[0];
sx q[0];
rz(-2.2472839) q[0];
sx q[0];
rz(-1.0067329) q[0];
rz(-pi) q[1];
rz(2.5587444) q[2];
sx q[2];
rz(-1.1839317) q[2];
sx q[2];
rz(-0.033535784) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.048679439) q[1];
sx q[1];
rz(-1.5041222) q[1];
sx q[1];
rz(2.9408022) q[1];
rz(-pi) q[2];
rz(1.7712426) q[3];
sx q[3];
rz(-2.5894937) q[3];
sx q[3];
rz(-2.828966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4422153) q[2];
sx q[2];
rz(-1.135004) q[2];
sx q[2];
rz(2.8507612) q[2];
rz(1.1902827) q[3];
sx q[3];
rz(-2.5255327) q[3];
sx q[3];
rz(-1.0440089) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46665835) q[0];
sx q[0];
rz(-0.062677296) q[0];
sx q[0];
rz(1.4516996) q[0];
rz(2.2845204) q[1];
sx q[1];
rz(-1.8430201) q[1];
sx q[1];
rz(-0.42795408) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1728522) q[0];
sx q[0];
rz(-1.4442872) q[0];
sx q[0];
rz(1.1304024) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0796135) q[2];
sx q[2];
rz(-2.2428838) q[2];
sx q[2];
rz(-0.88269688) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9792069) q[1];
sx q[1];
rz(-0.89079327) q[1];
sx q[1];
rz(1.1080145) q[1];
x q[2];
rz(-0.6537207) q[3];
sx q[3];
rz(-0.80721426) q[3];
sx q[3];
rz(-1.2098055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2032623) q[2];
sx q[2];
rz(-1.0469971) q[2];
sx q[2];
rz(-2.9948575) q[2];
rz(-2.9444368) q[3];
sx q[3];
rz(-1.4898172) q[3];
sx q[3];
rz(-2.3728235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7305304) q[0];
sx q[0];
rz(-2.125232) q[0];
sx q[0];
rz(-2.9732669) q[0];
rz(-0.39237818) q[1];
sx q[1];
rz(-2.7057251) q[1];
sx q[1];
rz(1.6900774) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0453059) q[0];
sx q[0];
rz(-1.1340965) q[0];
sx q[0];
rz(-0.023604579) q[0];
rz(-pi) q[1];
rz(-1.8658357) q[2];
sx q[2];
rz(-0.74271281) q[2];
sx q[2];
rz(-2.4969522) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5785006) q[1];
sx q[1];
rz(-0.76116409) q[1];
sx q[1];
rz(-0.20462402) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1216954) q[3];
sx q[3];
rz(-0.83415745) q[3];
sx q[3];
rz(2.4569805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3671941) q[2];
sx q[2];
rz(-0.44782475) q[2];
sx q[2];
rz(-1.9642584) q[2];
rz(-0.92159739) q[3];
sx q[3];
rz(-0.84601837) q[3];
sx q[3];
rz(2.6109076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0717936) q[0];
sx q[0];
rz(-1.256218) q[0];
sx q[0];
rz(-1.0569093) q[0];
rz(-2.9131556) q[1];
sx q[1];
rz(-0.7799131) q[1];
sx q[1];
rz(2.8905919) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91755501) q[0];
sx q[0];
rz(-0.38781181) q[0];
sx q[0];
rz(2.4075131) q[0];
rz(-2.2784581) q[2];
sx q[2];
rz(-2.302711) q[2];
sx q[2];
rz(-1.9307435) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31902203) q[1];
sx q[1];
rz(-0.21253702) q[1];
sx q[1];
rz(-0.6937278) q[1];
rz(-pi) q[2];
rz(-1.1993221) q[3];
sx q[3];
rz(-1.087093) q[3];
sx q[3];
rz(-0.35292816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2061657) q[2];
sx q[2];
rz(-2.3212104) q[2];
sx q[2];
rz(2.7247735) q[2];
rz(2.5721512) q[3];
sx q[3];
rz(-2.0408401) q[3];
sx q[3];
rz(-2.4429564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.57132974) q[0];
sx q[0];
rz(-1.8710192) q[0];
sx q[0];
rz(-0.50462333) q[0];
rz(0.2568256) q[1];
sx q[1];
rz(-0.80373126) q[1];
sx q[1];
rz(-1.9907192) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68783142) q[0];
sx q[0];
rz(-2.6021299) q[0];
sx q[0];
rz(-1.9786754) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1033642) q[2];
sx q[2];
rz(-2.4515477) q[2];
sx q[2];
rz(-2.0375117) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0432277) q[1];
sx q[1];
rz(-0.72611626) q[1];
sx q[1];
rz(0.67732201) q[1];
x q[2];
rz(2.7959149) q[3];
sx q[3];
rz(-0.44000235) q[3];
sx q[3];
rz(-2.1328762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5095832) q[2];
sx q[2];
rz(-0.42319599) q[2];
sx q[2];
rz(-0.43935856) q[2];
rz(-2.0158694) q[3];
sx q[3];
rz(-1.537354) q[3];
sx q[3];
rz(1.2996947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87089649) q[0];
sx q[0];
rz(-2.3211711) q[0];
sx q[0];
rz(-2.6743555) q[0];
rz(-1.0386508) q[1];
sx q[1];
rz(-0.50913441) q[1];
sx q[1];
rz(-1.4899303) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2307546) q[0];
sx q[0];
rz(-1.0808792) q[0];
sx q[0];
rz(1.3719489) q[0];
rz(-pi) q[1];
rz(0.57180228) q[2];
sx q[2];
rz(-2.9677113) q[2];
sx q[2];
rz(1.7553365) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.4626325) q[1];
sx q[1];
rz(-0.81607807) q[1];
sx q[1];
rz(-0.036463375) q[1];
x q[2];
rz(-2.0155679) q[3];
sx q[3];
rz(-2.091134) q[3];
sx q[3];
rz(-0.57131469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6963639) q[2];
sx q[2];
rz(-0.35564056) q[2];
sx q[2];
rz(-1.5871555) q[2];
rz(2.1571531) q[3];
sx q[3];
rz(-0.90960228) q[3];
sx q[3];
rz(1.4136774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5928891) q[0];
sx q[0];
rz(-2.2991572) q[0];
sx q[0];
rz(2.5133251) q[0];
rz(0.20333044) q[1];
sx q[1];
rz(-1.1140099) q[1];
sx q[1];
rz(-0.59439739) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3056476) q[0];
sx q[0];
rz(-1.3653269) q[0];
sx q[0];
rz(0.97272689) q[0];
rz(-pi) q[1];
rz(1.3514694) q[2];
sx q[2];
rz(-2.0989053) q[2];
sx q[2];
rz(-2.2161432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0815142) q[1];
sx q[1];
rz(-1.8412672) q[1];
sx q[1];
rz(0.31281506) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10108642) q[3];
sx q[3];
rz(-0.6154664) q[3];
sx q[3];
rz(-0.044199889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2079042) q[2];
sx q[2];
rz(-0.7242569) q[2];
sx q[2];
rz(-2.9653463) q[2];
rz(1.8267953) q[3];
sx q[3];
rz(-2.3254471) q[3];
sx q[3];
rz(0.76752457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15277302) q[0];
sx q[0];
rz(-1.7012699) q[0];
sx q[0];
rz(-1.758601) q[0];
rz(-0.11463166) q[1];
sx q[1];
rz(-0.75420598) q[1];
sx q[1];
rz(2.4892714) q[1];
rz(0.15936511) q[2];
sx q[2];
rz(-0.2401611) q[2];
sx q[2];
rz(-2.9693749) q[2];
rz(0.38269855) q[3];
sx q[3];
rz(-1.1094339) q[3];
sx q[3];
rz(-1.4239428) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
