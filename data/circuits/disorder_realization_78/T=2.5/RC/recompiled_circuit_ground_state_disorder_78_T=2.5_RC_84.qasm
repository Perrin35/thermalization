OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.652777) q[0];
sx q[0];
rz(4.1069694) q[0];
sx q[0];
rz(10.341636) q[0];
rz(-1.7011473) q[1];
sx q[1];
rz(-2.1135766) q[1];
sx q[1];
rz(-2.4457959) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3549933) q[0];
sx q[0];
rz(-2.7114641) q[0];
sx q[0];
rz(-1.7010078) q[0];
rz(0.90593018) q[2];
sx q[2];
rz(-3.0810789) q[2];
sx q[2];
rz(-1.1031594) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3506608) q[1];
sx q[1];
rz(-0.41484851) q[1];
sx q[1];
rz(2.7069339) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0761818) q[3];
sx q[3];
rz(-1.2418223) q[3];
sx q[3];
rz(-1.9291725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9838788) q[2];
sx q[2];
rz(-0.82008755) q[2];
sx q[2];
rz(-0.90041089) q[2];
rz(-2.3484717) q[3];
sx q[3];
rz(-1.6499062) q[3];
sx q[3];
rz(0.69860631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5291491) q[0];
sx q[0];
rz(-1.9731033) q[0];
sx q[0];
rz(0.68552619) q[0];
rz(-0.43288747) q[1];
sx q[1];
rz(-1.2665749) q[1];
sx q[1];
rz(-1.8276851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0746611) q[0];
sx q[0];
rz(-1.0767796) q[0];
sx q[0];
rz(-0.95290174) q[0];
rz(-pi) q[1];
rz(1.0772583) q[2];
sx q[2];
rz(-0.80776513) q[2];
sx q[2];
rz(1.1459717) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6095143) q[1];
sx q[1];
rz(-2.0348961) q[1];
sx q[1];
rz(2.1183234) q[1];
rz(-pi) q[2];
rz(0.39137381) q[3];
sx q[3];
rz(-2.5487566) q[3];
sx q[3];
rz(-2.9212869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96585298) q[2];
sx q[2];
rz(-2.7312835) q[2];
sx q[2];
rz(-2.8941594) q[2];
rz(-0.11582173) q[3];
sx q[3];
rz(-1.7146866) q[3];
sx q[3];
rz(-0.57870948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9847357) q[0];
sx q[0];
rz(-2.7655089) q[0];
sx q[0];
rz(1.0648741) q[0];
rz(-0.6330601) q[1];
sx q[1];
rz(-1.9828911) q[1];
sx q[1];
rz(-3.0050468) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9158123) q[0];
sx q[0];
rz(-1.5951711) q[0];
sx q[0];
rz(-2.4052785) q[0];
rz(2.2011287) q[2];
sx q[2];
rz(-2.5879938) q[2];
sx q[2];
rz(2.3950837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3202793) q[1];
sx q[1];
rz(-1.1007231) q[1];
sx q[1];
rz(2.4072985) q[1];
rz(2.7259215) q[3];
sx q[3];
rz(-1.9686254) q[3];
sx q[3];
rz(1.3597387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.69918767) q[2];
sx q[2];
rz(-1.2523315) q[2];
sx q[2];
rz(-1.5999954) q[2];
rz(1.0934746) q[3];
sx q[3];
rz(-1.6396061) q[3];
sx q[3];
rz(2.7961965) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71444756) q[0];
sx q[0];
rz(-0.26145014) q[0];
sx q[0];
rz(2.5647822) q[0];
rz(0.72751865) q[1];
sx q[1];
rz(-2.0501761) q[1];
sx q[1];
rz(-2.2732546) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24791323) q[0];
sx q[0];
rz(-0.94847877) q[0];
sx q[0];
rz(-1.076368) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98855726) q[2];
sx q[2];
rz(-1.120627) q[2];
sx q[2];
rz(-0.52272138) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35981813) q[1];
sx q[1];
rz(-1.9137772) q[1];
sx q[1];
rz(0.69056679) q[1];
x q[2];
rz(2.6336977) q[3];
sx q[3];
rz(-2.1334071) q[3];
sx q[3];
rz(0.25319448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.19520983) q[2];
sx q[2];
rz(-1.3888487) q[2];
sx q[2];
rz(0.69551224) q[2];
rz(-2.9269311) q[3];
sx q[3];
rz(-1.5011468) q[3];
sx q[3];
rz(-2.3282839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5401841) q[0];
sx q[0];
rz(-0.61487991) q[0];
sx q[0];
rz(2.0414798) q[0];
rz(2.9310138) q[1];
sx q[1];
rz(-2.5081878) q[1];
sx q[1];
rz(1.4932154) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8029047) q[0];
sx q[0];
rz(-0.7006104) q[0];
sx q[0];
rz(-2.3466097) q[0];
rz(-3.0576747) q[2];
sx q[2];
rz(-1.7482107) q[2];
sx q[2];
rz(2.5179297) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0303585) q[1];
sx q[1];
rz(-0.514563) q[1];
sx q[1];
rz(-1.2505934) q[1];
rz(-pi) q[2];
rz(-2.4744002) q[3];
sx q[3];
rz(-1.3722739) q[3];
sx q[3];
rz(-0.39451724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7275927) q[2];
sx q[2];
rz(-1.4294581) q[2];
sx q[2];
rz(2.2445938) q[2];
rz(-1.2796848) q[3];
sx q[3];
rz(-1.0106267) q[3];
sx q[3];
rz(0.047209386) q[3];
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
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996662) q[0];
sx q[0];
rz(-2.5157295) q[0];
sx q[0];
rz(-0.33011398) q[0];
rz(-0.64078981) q[1];
sx q[1];
rz(-1.7658486) q[1];
sx q[1];
rz(0.97553387) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2084239) q[0];
sx q[0];
rz(-2.5481249) q[0];
sx q[0];
rz(3.0400601) q[0];
x q[1];
rz(-1.2519205) q[2];
sx q[2];
rz(-1.826573) q[2];
sx q[2];
rz(2.6915611) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2157877) q[1];
sx q[1];
rz(-0.5277718) q[1];
sx q[1];
rz(2.5592519) q[1];
rz(-0.18838547) q[3];
sx q[3];
rz(-1.9864161) q[3];
sx q[3];
rz(-1.3481026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8581802) q[2];
sx q[2];
rz(-2.8874669) q[2];
sx q[2];
rz(0.88968366) q[2];
rz(0.59600082) q[3];
sx q[3];
rz(-2.0723876) q[3];
sx q[3];
rz(3.0349351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0058873) q[0];
sx q[0];
rz(-1.1361253) q[0];
sx q[0];
rz(-1.1329875) q[0];
rz(-2.5021482) q[1];
sx q[1];
rz(-2.0638128) q[1];
sx q[1];
rz(2.1741672) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3302878) q[0];
sx q[0];
rz(-2.0795224) q[0];
sx q[0];
rz(2.9697493) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.899951) q[2];
sx q[2];
rz(-1.2094922) q[2];
sx q[2];
rz(-0.45223927) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.514545) q[1];
sx q[1];
rz(-2.8882898) q[1];
sx q[1];
rz(-2.6652418) q[1];
rz(-pi) q[2];
rz(-0.65770517) q[3];
sx q[3];
rz(-1.8175565) q[3];
sx q[3];
rz(1.9384428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3966668) q[2];
sx q[2];
rz(-1.0090088) q[2];
sx q[2];
rz(-2.1972556) q[2];
rz(2.5189597) q[3];
sx q[3];
rz(-0.98832744) q[3];
sx q[3];
rz(2.3545789) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866163) q[0];
sx q[0];
rz(-2.3353307) q[0];
sx q[0];
rz(-3.1373613) q[0];
rz(-1.4684234) q[1];
sx q[1];
rz(-2.3783042) q[1];
sx q[1];
rz(-0.17446336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3136351) q[0];
sx q[0];
rz(-2.9330755) q[0];
sx q[0];
rz(2.0329518) q[0];
x q[1];
rz(2.9315591) q[2];
sx q[2];
rz(-2.6843417) q[2];
sx q[2];
rz(-2.5041472) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4694388) q[1];
sx q[1];
rz(-1.6368333) q[1];
sx q[1];
rz(-0.87182893) q[1];
rz(1.8352083) q[3];
sx q[3];
rz(-2.8972662) q[3];
sx q[3];
rz(-0.78247386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8909495) q[2];
sx q[2];
rz(-3.0477016) q[2];
sx q[2];
rz(-0.54035652) q[2];
rz(3.12449) q[3];
sx q[3];
rz(-2.159724) q[3];
sx q[3];
rz(0.66177773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.0509725) q[0];
sx q[0];
rz(-0.80247387) q[0];
sx q[0];
rz(0.78659868) q[0];
rz(1.1057009) q[1];
sx q[1];
rz(-1.6603755) q[1];
sx q[1];
rz(-0.22122637) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3174789) q[0];
sx q[0];
rz(-2.2750018) q[0];
sx q[0];
rz(-0.368035) q[0];
rz(-pi) q[1];
rz(-1.3972261) q[2];
sx q[2];
rz(-1.6643057) q[2];
sx q[2];
rz(1.5154523) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7936642) q[1];
sx q[1];
rz(-1.778144) q[1];
sx q[1];
rz(-0.4936895) q[1];
rz(1.33696) q[3];
sx q[3];
rz(-0.89409308) q[3];
sx q[3];
rz(-1.2016344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1456387) q[2];
sx q[2];
rz(-0.24253878) q[2];
sx q[2];
rz(1.185574) q[2];
rz(1.1061741) q[3];
sx q[3];
rz(-1.0217051) q[3];
sx q[3];
rz(-2.1574028) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36158654) q[0];
sx q[0];
rz(-1.3573283) q[0];
sx q[0];
rz(2.4318045) q[0];
rz(-2.2478726) q[1];
sx q[1];
rz(-0.62785405) q[1];
sx q[1];
rz(-0.46700221) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119962) q[0];
sx q[0];
rz(-1.3043367) q[0];
sx q[0];
rz(-0.14273739) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94552083) q[2];
sx q[2];
rz(-1.7904803) q[2];
sx q[2];
rz(-1.4131119) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4255387) q[1];
sx q[1];
rz(-1.3322833) q[1];
sx q[1];
rz(-0.038360049) q[1];
rz(-2.4670842) q[3];
sx q[3];
rz(-2.149636) q[3];
sx q[3];
rz(1.9306435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3422157) q[2];
sx q[2];
rz(-0.90505427) q[2];
sx q[2];
rz(-0.13599914) q[2];
rz(-2.811725) q[3];
sx q[3];
rz(-1.0672652) q[3];
sx q[3];
rz(-0.29712591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3099404) q[0];
sx q[0];
rz(-2.0257873) q[0];
sx q[0];
rz(-2.5153487) q[0];
rz(-2.2533439) q[1];
sx q[1];
rz(-1.2445969) q[1];
sx q[1];
rz(0.1319763) q[1];
rz(-2.6425003) q[2];
sx q[2];
rz(-0.44514984) q[2];
sx q[2];
rz(0.48322074) q[2];
rz(1.1430525) q[3];
sx q[3];
rz(-2.2528867) q[3];
sx q[3];
rz(-2.5611655) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
